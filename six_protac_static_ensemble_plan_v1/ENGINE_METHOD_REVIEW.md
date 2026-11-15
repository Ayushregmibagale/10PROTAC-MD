# Engine-method qualification before production

The currently tested engine is an adapted pRosettaC workflow, not an exact reproduction of the published implementation. The inherited adapter uses a seeded random global-orientation search in place of PatchDock, PyRosetta local docking/repacking, and fixed-head ETKDG/UFF linker construction. These differences must be stated in any methods section and qualified by structural benchmarks.

Original pRosettaC uses PatchDock for global shape-complementarity docking and RosettaDock local refinement: https://github.com/LondonLab/PRosettaC and https://pmc.ncbi.nlm.nih.gov/articles/PMC7592117/ . Both heads start in supported binary binding modes. ATP-pocket constraints are distinct from prescribing the VHL–kinase orientation.

A bounded search found no working PatchDock binary in the inspected local paths or within five levels of /home/s1638308 on Ernie. The old local PatchDock.tgz is HTML, not an installed program. The current official download page requests identity/affiliation/email and consent to store information: https://bio3d.cs.huji.ac.il/webserver/patchdock/ . No form was submitted and no personal details were supplied. This does not prove absence from all shared installations.

Production choices remain open: qualify the explicitly described adapter, or benchmark an actual PatchDock/Rosetta alternative separately. Do not silently substitute an untested search, rename random sampling PatchDock, or change the running pilot snapshot. Adaptations need demonstrated recovery rather than inheriting the pRosettaC publication's performance.

The old seed42 reference anchor ranges include the native geometry: GSK215 native distance2.706Å is within1.345–5.345Å; MZ1 native8.350Å is within7.610–11.610Å. Direct native-geometry linker construction passes three seeds for both references. These controls do not establish protein-clash freedom or global docking recovery.

The Ernie pilot has completed and its raw data are being retrieved. Structural recovery, source-calibrated steric review and resource accounting remain required. Additional Ernie submissions require user approval when ready. No MD.
