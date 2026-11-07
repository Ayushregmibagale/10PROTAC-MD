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
