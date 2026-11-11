# Option1: class balance

User requested emphasis on the two non-degraded targets. Verified PyTorch2.10 BCEWithLogitsLoss documentation: pos_weight multiplies the positive term; pos_weight3 would emphasize the six degraded targets. Implemented balanced per-target weights N_train/(2*N_class), computed independently inside each seven-target LOPO training fold. The negative/positive ratios are6 or2.5 depending on the held-out class; no held-out label enters the computation. All old graphs/labels/model architecture and300-epoch seed42 protocol are preserved. Only loss weighting changes.

This fix changes optimization contributions, not sample size. It cannot guarantee learning the negative class from two proteins. No conventional predictive model, new docking or MD is run. E2 features are reviewed separately for Option2.
