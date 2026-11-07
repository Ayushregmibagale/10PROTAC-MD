"""Explicit E2 geometric descriptors; no kinetic or degradation probabilities."""
import numpy as np
from protac_qc import ValidationError


def coordinates(value,name):
    array=np.asarray(value,dtype=float)
    if array.ndim!=2 or array.shape[1]!=3 or not len(array) or not np.isfinite(array).all():
        raise ValidationError('Finite nonempty Nx3 coordinates required: '+name)
    return array


def fit_reference_points(reference,mobile,points):
    """Proper least-squares rigid transform, preserving point-frame convention."""
    reference=coordinates(reference,'reference');mobile=coordinates(mobile,'mobile')
    points=coordinates(points,'points')
    if reference.shape!=mobile.shape or len(reference)<3:
        raise ValidationError('Matched anchor coordinates with at least three atoms required')
    a=reference-reference.mean(0);b=mobile-mobile.mean(0)
    if np.linalg.matrix_rank(a)<2 or np.linalg.matrix_rank(b)<2:
        raise ValidationError('Anchor geometry is degenerate')
    left,_,right=np.linalg.svd(a.T@b)
    rotation=right.T@np.diag([1.,1.,np.sign(np.linalg.det(right.T@left.T))])@left.T
    fitted=a@rotation.T+mobile.mean(0)
    mapped=(points-reference.mean(0))@rotation.T+mobile.mean(0)
    return mapped,float(np.sqrt(np.mean(np.sum((fitted-mobile)**2,axis=1)))),rotation


def validate_time_alignment(metric_times,trajectory_times,tolerance_ns):
    metric=np.asarray(metric_times,dtype=float);trajectory=np.asarray(trajectory_times,dtype=float)
    if not np.isfinite(tolerance_ns) or tolerance_ns<0:
        raise ValidationError('Finite nonnegative timestamp tolerance required')
    if metric.ndim!=1 or trajectory.shape!=metric.shape or not len(metric):
        raise ValidationError('Metric and trajectory frame counts differ; never truncate')
    if not np.isfinite(metric).all() or not np.isfinite(trajectory).all() or np.any(np.diff(metric)<0) or np.any(np.diff(trajectory)<0):
        raise ValidationError('Metric/trajectory timestamps are nonfinite or decreasing')
    error=float(np.max(np.abs(metric-trajectory)))
    if error>tolerance_ns:raise ValidationError('Metric and trajectory timestamps differ; no positional pairing')
    return error


def geometric_reachability(points,lysines,thresholds,box=None):
    from MDAnalysis.lib.distances import distance_array
    points=coordinates(points,'E2 points');lysines=coordinates(lysines,'lysine NZ')
    cutoffs=np.asarray(thresholds,dtype=float)
    if cutoffs.ndim!=1 or not len(cutoffs) or not np.isfinite(cutoffs).all() or np.any(cutoffs<=0) or len(set(cutoffs))!=len(cutoffs):
        raise ValidationError('Unique finite positive reach thresholds required')
    distances=distance_array(points,lysines,box=box).min(axis=1)
    return float(distances.min()),[float(np.mean(distances<cutoff)) for cutoff in cutoffs]


class WholeProteinGuard:
    """Skip make_whole only when every bond is already shorter than 3 A.

    With every periodic cell height above 6 A, these bond displacements are
    already inside the minimum-image cell. Otherwise use the original unwrap.
    No reference atom/fragment is recentered by this operation.
    """
    def __init__(self,protein):
        from MDAnalysis.transformations import unwrap
        self.protein=protein;self.universe=protein.universe
        self.bonds=np.asarray(protein.bonds.indices,dtype=int)
        if not len(self.bonds) or not np.isin(self.bonds,protein.indices).all():
            raise ValidationError('Complete protein-only bond mapping required')
        self.fallback=unwrap(protein)

    def __call__(self,ts):
        from MDAnalysis.lib.mdamath import triclinic_vectors
        box=np.asarray(ts.dimensions,dtype=float)
        if box.shape!=(6,) or not np.isfinite(box).all():raise ValidationError('Finite periodic box required')
        vectors=triclinic_vectors(box,dtype=np.float64)
        volume=abs(float(np.linalg.det(vectors)))
        areas=np.array([np.linalg.norm(np.cross(vectors[(i+1)%3],vectors[(i+2)%3])) for i in range(3)])
        if volume<=0 or np.any(areas<=0):raise ValidationError('Nondegenerate periodic cell required')
        positions=self.universe.atoms.positions
        deltas=positions[self.bonds[:,0]]-positions[self.bonds[:,1]]
        if not np.isfinite(deltas).all():raise ValidationError('Nonfinite protein bond coordinates')
        if np.min(volume/areas)>6. and np.max(np.sum(deltas*deltas,axis=1))<9.:
            return False
        self.fallback(ts)
        return True


class WholeProteinTree:
    """Reconstruct bonded fragments using MIC edges and a fixed Euler tree.

    Root atoms remain fixed. Every bonded edge is checked after reconstruction;
    inconsistent cyclic geometry falls back to the original make_whole operation.
    This is a numerical implementation change requiring reference comparisons.
    """
    def __init__(self,protein):
        from MDAnalysis.transformations import unwrap
        self.protein=protein;self.fallback=unwrap(protein)
        ids=protein.indices;lookup=np.full(len(protein.universe.atoms),-1,dtype=int);lookup[ids]=np.arange(len(ids))
        bonds=lookup[protein.bonds.indices]
        if not len(bonds) or np.any(bonds<0):raise ValidationError('Complete protein bond groups required')
        self.bonds=bonds;adj=[[] for _ in ids]
        for i,j in bonds:adj[i].append(j);adj[j].append(i)
        parents=np.full(len(ids),-2,dtype=int);root_of=np.zeros(len(ids),dtype=int);entries=np.full(len(ids),-1,dtype=int);root_start={};children=[];parent_nodes=[];events=[];signs=[];edge_of={}
        for root in range(len(ids)):
            if parents[root]!=-2:continue
            parents[root]=-1;root_of[root]=root;root_start[root]=len(events)-1;stack=[(root,iter(adj[root]))]
            while stack:
                node,neighbors=stack[-1]
                try:child=next(neighbors)
                except StopIteration:
                    stack.pop()
                    if parents[node]>=0:events.append(edge_of[node]);signs.append(-1.)
                    continue
                if parents[child]!=-2:continue
                parents[child]=node;root_of[child]=root;edge_of[child]=len(children);children.append(child);parent_nodes.append(node);entries[child]=len(events);events.append(edge_of[child]);signs.append(1.);stack.append((child,iter(adj[child])))
        self.children=np.asarray(children,dtype=int);self.parents=np.asarray(parent_nodes,dtype=int);self.events=np.asarray(events,dtype=int);self.signs=np.asarray(signs);self.entries=entries;self.root_of=root_of;self.root_start=root_start

    def __call__(self,ts):
        from MDAnalysis.lib.distances import minimize_vectors
        box=np.asarray(ts.dimensions,dtype=float);raw=coordinates(self.protein.positions,'protein')
        if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0):raise ValidationError('Valid periodic cell required')
        deltas=minimize_vectors(raw[self.children]-raw[self.parents],box)
        prefix=np.cumsum(deltas[self.events]*self.signs[:,None],axis=0)
        whole=raw.copy()
        for root,start in self.root_start.items():
            children=self.children[self.root_of[self.children]==root]
            baseline=prefix[start] if start>=0 else 0.
            whole[children]=raw[root]+prefix[self.entries[children]]-baseline
        expected=minimize_vectors(raw[self.bonds[:,1]]-raw[self.bonds[:,0]],box)
        error=np.max(np.linalg.norm(whole[self.bonds[:,1]]-whole[self.bonds[:,0]]-expected,axis=1))
        if not np.isfinite(error) or error>.001:
            self.fallback(ts);return True
        self.protein.positions=whole
        return False
