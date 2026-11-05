# 01 — Scientific Rationale

## 1. What are PROTACs?

PROTACs (PROteolysis TArgeting Chimeras) are bifunctional small molecules that degrade a
protein of interest (POI) by hijacking the cell's ubiquitin-proteasome system. A PROTAC
molecule has three parts:

1. **POI warhead** — binds the target protein's active site (usually a known inhibitor)
2. **Linker** — connects the two warheads; flexibility/length determines the PROTAC's
   conformational space
3. **E3 ligase warhead** — recruits an E3 ubiquitin ligase (here: VHL)

When both warheads are simultaneously bound, a **ternary complex** (E3·PROTAC·POI) forms.
The E3 ligase then ubiquitinates a lysine on the POI surface, flagging it for proteasomal
degradation. The PROTAC is catalytic — once the POI is degraded it is released and can
recruit another POI molecule.

## 2. Why some proteins are not degraded despite tight PROTAC binding

A critical and counter-intuitive finding in the field is that **high binding affinity of the
PROTAC to the POI is necessary but not sufficient for degradation**. Non-degraders can show
nanomolar POI binding (Kd) while producing no measurable degradation (DC50 = undefined).

The main mechanistic explanations are:

1. **Hook effect / ternary complex geometry**: The ternary complex cannot form a
   productive geometry — the POI and E3 ligase surfaces are incompatible for simultaneous
   binding, or the complex dissociates faster than ubiquitination can occur.
2. **No accessible lysine**: The POI has no surface lysine within the ubiquitination reach
   of the E2 enzyme (~35–55 Å from the VHL-associated E2 active site cysteine). Ubiquitin
   is transferred from E2 to a POI lysine; if none is in range the reaction cannot occur.
3. **Cooperativity (alpha)**: Ternary complex cooperativity (α > 1 = cooperative, α < 1 =
   anticooperative) affects the equilibrium population of productive ternary complexes.
4. **Structural occlusion**: The POI surface near a potential ubiquitination lysine may be
   occluded in the ternary complex, preventing E2 access.

## 3. SJF8240: the PROTAC studied here

SJF8240 is a **pan-kinase PROTAC** built around:

- **VHL warhead**: A (2S,4R)-4-hydroxyproline-based VHL ligand (same chemotype as VH032,
  widely used in VHL-recruiting PROTACs). This warhead binds in the VHL β-sheet pocket
  with sub-micromolar affinity.
- **Kinase warhead**: The foretinib/XL880/GSK1363089 pharmacophore — a potent multi-kinase
  inhibitor that targets the ATP-binding site of many receptor tyrosine kinases (MET, AXL,
  DDR1/2, EphB2, VEGFR2) and some other kinases.
- **Linker**: 12-atom flexible PEG/alkyl chain connecting the two warheads.

SJF8240 has been experimentally profiled against a panel of kinases. For some kinases it
causes robust degradation; for others with equally tight binding it causes none. This makes
it an ideal system for studying the structural determinants of degradation.

### Experimental classification used in this study

| Class | Criteria | POIs |
|-------|----------|------|
| High-affinity degraders | Kd < 100 nM to PROTAC, DC50 measured | MET, DDR2, RIPK2, EPHB2, MAPK14 |
| High-affinity non-degraders | Kd < 100 nM to PROTAC, no degradation | AXL, SLK, ABL1, EPHA2, MAP4K5 |
| Low-affinity degraders | Kd > 100 nM, degradation observed | TNIK, PIP4K2C, CDK17 |

The high-affinity paired panel (5 degraders + 5 non-degraders, all with tight PROTAC
binding) is the most scientifically powerful comparison because affinity cannot explain the
difference — only ternary complex geometry can.

## 4. Why pRosettaC?

pRosettaC (Zaidman et al., *Nature Chemical Biology*, 2020) is the leading published
computational method for PROTAC ternary complex modelling. It:

1. **Generates diverse ternary complex geometries** by running rigid-body docking
   (PatchDock) of the E3·PROTAC half-complex against the POI, then applying local Rosetta
   relaxation to each docking pose.
2. **Builds a physically reasonable linker** by running conformer generation with
   pharmacophoric constraints anchoring the two warheads to their crystallographically
   determined binding poses.
3. **Clusters the resulting ensemble** to identify the dominant ternary complex geometries.
4. **Has been experimentally validated** — in the original paper it predicted productive
   geometries for known degraders and found no productive geometries for known non-degraders.

pRosettaC is appropriate here because:
- We have co-crystal structures (or very close homologues) for all POIs
- The warheads' binding poses are well-characterised
- We have a clear degrader/non-degrader ground truth to train against
- The method produces interpretable structural features (interface contacts, E2–Lys
  distances) that can be used in a simple ML classifier

## 5. Hypothesis and expected discriminating features

We hypothesise that the following features will differ between degraders and non-degraders:

1. **Interface contacts (5Å cutoff)**: Productive ternary complexes require a stable
   POI–E3 interface. More contacts = more stable complex = higher probability of
   ubiquitination before complex dissociation.
2. **E2–lysine distance**: The active cysteine of the E2 enzyme (UBE2D3/UbcH5c in the
   VHL pathway) must be within ~35–55 Å of an accessible POI lysine to catalyse
   ubiquitin transfer. POIs with no lysine in this range cannot be ubiquitinated.
3. **Number of E2-optimal lysines**: More lysines in the optimal distance range increases
   the probability that at least one will be ubiquitinated.
4. **Polar hydrogen bond pairs at interface**: Specific H-bonds stabilise the ternary
   complex and may be needed for the complex to persist long enough for ubiquitination.
5. **Ternary complex score (REU)**: Lower Rosetta energy = more favourable complex.

## 6. Expected output and downstream use

After all 13 POI runs complete, we will have:

1. **Structural ensemble** for each POI: ~6–20 cluster representatives per POI
2. **Feature table** (`ml_features_best.csv`): one row per POI, ~50 structural/geometric
   features
3. **Classification model**: Leave-one-out logistic regression + decision tree trained on
   the 10 high-affinity POIs (5+5), with the 3 low-affinity degraders as an independent
   test set
4. **Mechanistic insight**: Which structural features (and therefore which aspects of
   ternary complex geometry) explain degradation outcomes

This work is designed to provide:
- A retrospective validation of pRosettaC on a well-characterised PROTAC
- A set of interpretable structural rules for predicting degradation
- A foundation for prospective PROTAC design optimisation
