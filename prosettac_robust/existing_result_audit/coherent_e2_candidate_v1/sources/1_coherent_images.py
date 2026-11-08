"""Candidate whole-component imaging conventions for retained E2 analysis.

Translations are lattice vectors only. This defines a geometric representation,
not a productive-complex test, kinetic model or automatic training release.
"""
import numpy as np
from protac_qc import ValidationError
from e2_geometry import coordinates,fit_reference_points


def cell_vectors(box):
    from MDAnalysis.lib.mdamath import triclinic_vectors
    box=np.asarray(box,dtype=float)
    if box.shape!=(6,) or not np.isfinite(box).all() or np.any(box[:3]<=0) or np.any((box[3:]<=0)|(box[3:]>=180)):
        raise ValidationError('Finite physical triclinic box required')
    cosine=np.cos(np.deg2rad(box[3:]))
    gram_factor=1+2*np.prod(cosine)-np.sum(cosine*cosine)
    if gram_factor<=1e-12:
        raise ValidationError('Degenerate or numerically singular periodic cell')
    vectors=triclinic_vectors(box,dtype=np.float64)
    if np.linalg.det(vectors)<=np.prod(box[:3])*1e-12:
        raise ValidationError('Nondegenerate right-handed periodic cell required')
    return vectors


def lattice_shift(delta,box):
    """Translation moving a displacement to its nearest periodic image."""
    from MDAnalysis.lib.distances import minimize_vectors
    delta=np.asarray(delta,dtype=float);vectors=cell_vectors(box)
    if delta.shape!=(3,) or not np.isfinite(delta).all():raise ValidationError('Finite image displacement required')
    reduced=minimize_vectors(delta[None],np.asarray(box,dtype=float))[0]
    coefficients=np.linalg.solve(vectors.T,reduced-delta);integers=np.rint(coefficients)
    if not np.allclose(coefficients,integers,atol=5e-5,rtol=0):raise ValidationError('Computed image shift is not a lattice translation')
    return integers@vectors,integers.astype(np.int64)


def image_to_reference(mobile,reference,box):
    mobile=coordinates(mobile,'image mobile');reference=coordinates(reference,'image reference')
    if mobile.shape!=reference.shape:raise ValidationError('Matched reference atoms required for image selection')
    shift,integers=lattice_shift(mobile.mean(0)-reference.mean(0),box)
    return mobile+shift,shift,integers


def image_by_closest_pair(anchor,mobile,box):
    from MDAnalysis.lib.distances import distance_array
    anchor=coordinates(anchor,'contact anchor');mobile=coordinates(mobile,'contact mobile')
    distances=distance_array(anchor,mobile,box=np.asarray(box,dtype=float));i,j=np.unravel_index(np.argmin(distances),distances.shape)
    shift,integers=lattice_shift(mobile[j]-anchor[i],box)
    return mobile+shift,shift,integers,(int(i),int(j))


def assemble(anchor,vhl,poi,reference_anchor,reference_vhl,reference_poi,box):
    """Return reference-anchored VHL and closest-periodic-pair POI images.

    A separate reference-pose POI image supplies an agreement diagnostic.
    The virtual E2 prior is intentionally absent from the selection rules.
    """
    cell_vectors(box)
    vhl=coordinates(vhl,'VHL');poi=coordinates(poi,'POI')
    predicted,rmsd,rotation=fit_reference_points(reference_anchor,anchor,np.vstack([reference_vhl,reference_poi]))
    predicted_vhl=predicted[:len(reference_vhl)];predicted_poi=predicted[len(reference_vhl):]
    lifted_vhl,shift_vhl,ivhl=image_to_reference(vhl,predicted_vhl,box)
    lifted_poi,shift_poi,ipoi,pair=image_by_closest_pair(lifted_vhl,poi,box)
    reference_poi_image,ref_shift,iref=image_to_reference(poi,predicted_poi,box)
    agreement=bool(np.linalg.norm(shift_poi-ref_shift)<.001)
    return {'vhl':lifted_vhl,'poi':lifted_poi,'reference_poi':reference_poi_image,'vhl_shift_A':shift_vhl,'poi_shift_A':shift_poi,'reference_poi_shift_A':ref_shift,'vhl_lattice_shift':ivhl,'poi_lattice_shift':ipoi,'reference_poi_lattice_shift':iref,'POI_image_methods_agree':agreement,'closest_pair_CA_indices':pair,'ELOC_fit_RMSD_A':rmsd,'VHL_reference_fit_RMSD_A':float(np.sqrt(np.mean(np.sum((lifted_vhl-predicted_vhl)**2,axis=1)))),'rotation_determinant':float(np.linalg.det(rotation)),'automatic_training_release':False}
