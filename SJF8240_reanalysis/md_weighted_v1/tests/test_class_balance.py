from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
import torch,pytest
from common import balanced_class_weights,balanced_bce
LABELS={'a':0,'b':0,'c':1,'d':1,'e':1,'f':1,'g':1,'h':1}
@pytest.mark.parametrize('held,ratio',[('a',6.0),('c',2.5)])
def test_training_classes_have_equal_total_weight(held,ratio):
 train=[p for p in LABELS if p!=held];w=balanced_class_weights(LABELS,train)
 assert w[0]/w[1]==pytest.approx(ratio)
 assert sum(w[LABELS[p]] for p in train if LABELS[p]==0)==pytest.approx(sum(w[LABELS[p]] for p in train if LABELS[p]==1))
def test_held_out_label_cannot_affect_weights():
 train=[p for p in LABELS if p!='a'];changed=dict(LABELS,a=1)
 assert balanced_class_weights(LABELS,train)==balanced_class_weights(changed,train)
def test_negative_gradient_is_upweighted_not_positive():
 w=balanced_class_weights(LABELS,list(LABELS));a=torch.tensor([0.],requires_grad=True);b=torch.tensor([0.],requires_grad=True)
 balanced_bce(a,0,w).backward();balanced_bce(b,1,w).backward()
 assert abs(a.grad.item())/abs(b.grad.item())==pytest.approx(3.)
