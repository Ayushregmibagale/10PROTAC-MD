# Retained GNN model provenance

Read-only, CPU `torch.load(weights_only=True)` inspection found two original MGDB pretraining checkpoints, `pretrain_best.pt` and `pretrain_final.pt`, and their training log. Both deserialize successfully and all tensor values are finite. Original hashes were checked; no model was instantiated or trained.

These weights belong to InterfaceGNN with 24 node features and one edge feature, referenced by the separate 04f adapter experiment. They are not the main 04g E2 dual-encoder fold states (nine static node features, three dynamic node features and four E2 scalars). Finding these weights therefore corrects any blanket statement that no trained weights exist, while leaving exact replay of the historical seventeen-row E2 evaluation unsupported by the inspected artifacts.

Evidence, tensor metadata, source hashes and limitations are in `existing_result_audit/retained_gnn_weights_v1`. No new predictions, biological training or conventional ML model was added.
