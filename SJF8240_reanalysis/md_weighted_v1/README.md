# Fold-balanced MD-only GNN

Option1: same eight MD targets and model as corrected_md_only_v1, with fold-local balanced BCE weights. Degrader=1; ABL1/AXL=0. pos_weight=3 would upweight degraders, so use explicit weights N_train/(2*N_class). Negative/positive weight ratios are6 when a negative is held out,2.5 when a positive is held out. Compute from training targets only. Original300epochs, seed42, original MD-only class, no tuning. ELOB control and true VHL. No new docking or MD.
