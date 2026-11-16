# Six-compound Ernie technical-pilot preparation

All6technical-control input directories are staged with supplier molecule identity and source-file hashes. INPUT_AUDIT.json records all-file integrity and independent head/protein steric diagnostics including explicit MS143 phosphoThr308atom selection. These are technical controls, not full kinome coverage or new biological labels. SJF_delta/MET is only a parent-head technical control.

Current status: HOLD, not submitted. Reference qualification is running locally. Full input state/sterics qualification and MS143 end-to-end PTM-aware post-QC remain required. The draft six-task configuration uses4cores,24GiB RSS,24GiB cumulative output,50GiB minimum free,4h cap pertask. Those are caps, not runtime predictions or allocated resources.

A private executable draft is at /tmp/six_protac_ernie_pilot_draft_v1. Its controller refuses execution at the HOLD gate before creating runs or STATUS; DRAFT_GATE_CHECK.json records this negative control. Do not edit the gate to bypass qualification; build a new qualified frozen release after all checks pass.

REMOTE_ENVIRONMENT_CHECK.json confirms read-only Ernie connectivity,96visibleCPUs and about1.49TBfree at the timechecked. No scheduler commands appeared on its noninteractive PATH; this does not establish site allocation/scheduling policy or permission to occupy allCPUs. No job was submitted. Final concrete submission goes to the user for approval when qualified.
