"""Sequence-based identification of observed protein C-alpha residues.

Coverage is of the observed construct, not canonical full-length completeness.
"""
from Bio.Align import PairwiseAligner,substitution_matrices
from protac_qc import ValidationError
ALIGNER=PairwiseAligner();ALIGNER.mode='local';ALIGNER.substitution_matrix=substitution_matrices.load('BLOSUM62');ALIGNER.open_gap_score=-10;ALIGNER.extend_gap_score=-.5


def identify(sequence,entries):
    matches=[]
    for uid,target in entries.items():
        alignment=ALIGNER.align(target,sequence)[0]
        pairs=[(i,j) for (a,b),(c,d) in zip(*alignment.aligned) for i,j in zip(range(a,b),range(c,d))]
        if not pairs:continue
        identity=sum(target[i]==sequence[j] for i,j in pairs)/len(pairs);coverage=len(pairs)/len(sequence)
        if identity>=.90 and coverage>=.90:matches.append({'uniprot':uid,'identity':identity,'observed_CA_coverage':coverage,'canonical_first':min(i for i,j in pairs)+1,'canonical_last':max(i for i,j in pairs)+1})
    if len(matches)!=1:raise ValidationError('Protein role not uniquely identified with >=90% identity and observed-sequence coverage')
    return matches[0]
