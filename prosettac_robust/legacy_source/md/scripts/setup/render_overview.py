# Visual inspection renders for the 8 READY ternary complexes + MET cluster comparison
# Run: pymol -cq md/scripts/setup/render_overview.py
# Colour scheme: VHL = grey70, ElonginB = skyblue, ElonginC = palecyan,
#                POI = orange (degrader) / teal (non-degrader), PROTAC (PTC) = magenta sticks

from pymol import cmd

base = "/home/ayushregmibagale/PROTAC degradation Prediction/pRosettaC/runs"
outdir = "/home/ayushregmibagale/PROTAC degradation Prediction/md/figures"

pois = {
    "MET":    ("high_affinity_degraders",     "cluster1", "combined_48_15_0001.pdb",  "degrader"),
    "DDR2":   ("high_affinity_degraders",     "cluster1", "combined_38_41_0001.pdb",  "degrader"),
    "RIPK2":  ("high_affinity_degraders",     "cluster1", "combined_99_30_0001.pdb",  "degrader"),
    "EPHB2":  ("high_affinity_degraders",     "cluster1", "combined_38_48_0001.pdb",  "degrader"),
    "MAPK14": ("high_affinity_degraders",     "cluster1", "combined_67_50_0001.pdb",  "degrader"),
    "AXL":    ("high_affinity_no_degradation","cluster1", "combined_10_45_0001.pdb",  "non-degrader"),
    "ABL1":   ("high_affinity_no_degradation","cluster1", "combined_177_38_0001.pdb", "non-degrader"),
    "EPHA2":  ("high_affinity_no_degradation","cluster1", "combined_134_42_0001.pdb", "non-degrader"),
}

cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)

# ---------------------------------------------------------------------------
# 1) One overview render per POI (cluster1 ternary complex)
# ---------------------------------------------------------------------------
for poi, (cat, cluster, pdb, label) in pois.items():
    path = f"{base}/{cat}/{poi}/Results/{cluster}/{pdb}"
    cmd.load(path, poi)
    cmd.hide("everything", poi)
    cmd.show("cartoon", f"{poi} and polymer")
    cmd.color("grey70", f"{poi} and chain A and resi 1-103")     # VHL
    cmd.color("skyblue", f"{poi} and chain B")                    # ElonginB
    cmd.color("palecyan", f"{poi} and chain C")                   # ElonginC
    poi_color = "orange" if label == "degrader" else "teal"
    cmd.color(poi_color, f"{poi} and chain A and resi 104-9999")  # POI
    cmd.show("sticks", f"{poi} and resn PTC")
    cmd.color("magenta", f"{poi} and resn PTC")
    cmd.orient(poi)
    cmd.zoom(poi, buffer=3)
    cmd.png(f"{outdir}/{poi}_cluster1.png", width=600, height=600, dpi=150, ray=1)
    cmd.delete(poi)

# ---------------------------------------------------------------------------
# 2) MET cluster1 vs cluster2 vs cluster6 comparison, aligned on ElonginB+C
#    (the "fixed" E3 side) to show how POI+PROTAC pose differs between clusters
# ---------------------------------------------------------------------------
met_clusters = {
    "cluster1": ("combined_48_15_0001.pdb", "cyan"),
    "cluster2": ("combined_21_38_0001.pdb", "magenta"),
    "cluster6": ("combined_16_4_0001.pdb",  "yellow"),
}

objs = []
for cl, (pdb, color) in met_clusters.items():
    path = f"{base}/high_affinity_degraders/MET/Results/{cl}/{pdb}"
    objname = f"MET_{cl}"
    cmd.load(path, objname)
    objs.append((objname, color))

ref = objs[0][0]
for objname, _ in objs[1:]:
    cmd.align(f"{objname} and chain B+C", f"{ref} and chain B+C")

cmd.hide("everything")
# shared VHL + ElonginB/C reference frame from cluster1, in grey
ref_sel = f"{ref} and (chain B+C or (chain A and resi 1-103))"
cmd.show("cartoon", ref_sel)
cmd.color("grey80", ref_sel)

for objname, color in objs:
    poi_sel = f"{objname} and chain A and resi 104-9999"
    cmd.show("cartoon", poi_sel)
    cmd.color(color, poi_sel)
    cmd.show("sticks", f"{objname} and resn PTC")
    cmd.color(color, f"{objname} and resn PTC")

cmd.orient(ref)
cmd.zoom("all", buffer=5)
cmd.png(f"{outdir}/MET_cluster_comparison.png", width=900, height=900, dpi=150, ray=1)

print("Done. Renders written to md/figures/")
