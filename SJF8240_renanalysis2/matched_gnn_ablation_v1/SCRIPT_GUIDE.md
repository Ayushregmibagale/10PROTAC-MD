# Ordered scripts

Use /home/ayushregmibagale/miniconda3/bin/python3. Dependencies:PyTorch,PyG,NumPy,pandas,SciPy,scikit-learn,matplotlib,python-docx. The reviewed native-PyG guard disables incompatible optional compiled modules before PyG import. Do not import PyG first or run Python with optimization flags.

| Order | Entry point | Purpose |
|---|---|---|
|01|01_preflight.py|Verify16 matchedgraphinputs,whole-target splits,currentuserlabels and architecture equivalence.|
|02|02_train_matched_pairs.py --workers 2|Train24fixed300epochfolds:12pairs×2variants. Existingfolds require matchingcode/input/checkpointhashes;orphanweights cause explicitfailure rather than overwrite. --limit isdiagnostic only.|
|03|03_verify_cpu_replay.py|Separately recreateinitialstates andloadall24checkpoints onCPU;replay48scores andcheckheld-out exclusions.|
|04|04_review_experimental_evidence.py|Read-only LF/PS review for8targets;keep source versions/contrasts distinct;no relabeling or statisticalrefit.|
|05|05_compare_ablation.py|Compare all12 matchedpairs andboth negatives atfixed0.5;export allscores,margins,per-target summaries andchangedpairs.|
|06|06_write_findings.py|Write finalWord/Markdown report andPNG/PDF comparison.|
|07|07_check_package.py|Verify allweights,sourcereceipts,CPUreplays andsavedoutputs.|
|08|08_record_runtime.py|Save software versions andnative-backend configuration.|

Allsteps havealreadycompleted. A new method,seed,loss orinput selection requires a newprotocol/version. Never edit frozenmanifests toaccept stale results. The primary pairedAUROC isnot the pooledROC ofprobabilities produced bydifferentmodels. Per-targetmeans/min/max aredescriptive anddonotestablish calibratedclassprobabilities.

The original29features andeveryedge areunchanged. Addedfivecolumns arefixednative-reference lysine descriptors. The controlzero-pads thosecolumns;bothmodels haveidentical34inputarchitecture andpairedinitialstate hashes. Eachheld-out pair leaves5positives+1negative fortraining. Eachprotein isheldout withallitsposes/windows;trainingdoesnotuseheld-outtargetlosses ornormalization.

No lossweighting,thresholdcalibration orbest-checkpoint selection wasadded. The300epochfinalcheckpoint isalwaysused. Allfoldweights andimmediatereplayreceipts areinresults/<variant>/. Independent verification isvalidation/independent_cpu_verification.json. The source experiment,historical weights anduserlabels areunchanged. Experimental sourcefiles andtheold reviewed modeldefinition remainread-onlydependencies forprovenance checks;copiedgraphs/model/runtime needed fortraining areinthisfolder.

TheRPE1LF `_new` version doesnotindependentlyestablish everybinarylabel:EPHA2 isinconclusive there. Thehistorical computationalchemistry differsfromsupplier-matchedSJF8240. These limitations remain explicit;theablation isnot independentbiological validation orpermissionforproductionkinomepredictions.
