# Active six-compound static campaign

Latest user instruction, 15 November 2026: continue the project, update GitHub and prioritize pRosettaC on Ernie. All six compounds were purchased from Tocris/Bio-Techne. MS143 (8081) is now identity-verified in the registry.

Explicit constraint: WE WILL NOT DO MD FOR THESE SIX COMPOUNDS. Use pRosettaC ensembles, static structural analysis and later an ensemble GNN. Existing historical MD is separate supporting material only. Do not resume the old MD pilot or require MD equilibration as a docking release gate.

Next: reconcile stopped native-recovery outputs, prepare MS143 bound-head mappings, verify Ernie access/runtime/storage and release an isolated repaired docking pilot. Preserve scientific QC, raw failure counts, chemical identity and legacy outputs. No failed finite search is a biological negative. Save regular source/report checkpoints to GitHub.

Read PLAN.md. The old SJF8240 supplier-correct pilot remains historical/paused; its docking outputs can be inventoried/reused without restarting its MD stages.

## Checkpoint: 15 November 2026, static reference pilot

Read RESEARCH_CONTRACT.md before changing protocol or releasing data. User explicitly prioritizes scientific quality; no MD for these compounds. Old completed small-budget references have now been assessed: MZ1 misses the frozen coarse gate; GSK215 has one passing retained pose but top1fails. Full campaign release remains pending.

MS143 source/heads are under MS143/ and references/. Its modified AKT head is modeled, not experimental; TPO308is retained in the parent protein. Do not silently submit it as a validated ternary pose.

An isolated bounded reference/resource pilot was launched on Ernie at /home/s1638308/ernie_static_reference_pilot_v1, supervisor PID313406. Inspect remote STATUS.json and receipts before any retry. Four-core affinity,24GiB summed process RSS cap,10GiB output cap,50GiBminimum free disk; native integration cap30min, each of two global runs cap8h. Global budget50×10×20,seed42. STOP sentinel terminates the current task. This pilot is not an authorization to release failed geometry or start kinome production. Licensed dependencies are reused privately from the remote protac environment.

Commands: ssh ernie 'cat /home/s1638308/ernie_static_reference_pilot_v1/STATUS.json'; ssh ernie 'touch /home/s1638308/ernie_static_reference_pilot_v1/STOP'. Monitor/update receipts and GitHub while actively working; no claim of continuous human/agent monitoring after the turn.

Latest result: native-frame Ernie control completed with0/3feasible linkers at20conformers; this is NOT a successful native-recovery certificate. Global GSK215 benchmark is running; MZ1 follows if execution succeeds. This is a budget diagnostic, not a scientific production release. Initial MS143 lowest-internal-energy candidate fails the independent vdW screen; strict_review_v2retains one passing seed4candidate (maximum overlap0.576Åversus deposited parent0.483Å). Do not use the initial root candidate as input. All failed trials remain recorded.

Latest user instruction: keep working autonomously; ask only when ready to submit the Ernie campaign. No additional remote submission beyond the already-running pilot queue without that approval. Current editable checkout: /mnt/c/Users/Ayush Regmi Bagale/10PROTAC_static_campaign_work . The original Linux working tree is preserved, including its staged checkpoint and unrelated changes. GitHub updates come from this isolated checkout.
