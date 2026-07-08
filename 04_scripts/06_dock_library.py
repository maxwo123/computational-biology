"""
Steps 15-18: high-throughput docking of the full library against all three
receptors (HIV-1 protease [primary], Cathepsin D + BACE1 [human off-targets]),
and extraction of the best binding affinity (mode 1) from each Vina run.

Each ligand x receptor pose file is cached; re-running skips completed dockings.

Run:  .venv/bin/python 04_scripts/06_dock_library.py
"""
import csv, os, subprocess, time

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
DOCK = os.path.join(BASE, "03_docking")
LIBDIR = os.path.join(BASE, "02_ligands", "library")
POSES = os.path.join(DOCK, "poses")
os.makedirs(POSES, exist_ok=True)

RECEPTORS = ["hiv_protease", "cathepsin_d", "bace1"]
library = list(csv.DictReader(open(os.path.join(BASE, "02_ligands", "library.csv"))))

def best_affinity(vina_stdout):
    for line in vina_stdout.splitlines():
        s = line.split()
        if len(s) >= 4 and s[0] == "1":
            return float(s[1])
    return None

results = []
t0 = time.time()
total = len(RECEPTORS) * len(library)
done = 0
for rec in RECEPTORS:
    conf = os.path.join(DOCK, f"{rec}_vina.conf")
    for lig in library:
        lid = lig["ligand_id"]
        ligpdbqt = os.path.join(LIBDIR, f"{lid}.pdbqt")
        out = os.path.join(POSES, f"{rec}__{lid}.pdbqt")
        log = os.path.join(POSES, f"{rec}__{lid}.log")
        done += 1
        if os.path.exists(out) and os.path.exists(log):
            aff = best_affinity(open(log).read())
            results.append({"ligand_id": lid, "name": lig["name"], "status": lig["status"],
                            "receptor": rec, "affinity": aff})
            print(f"[{done}/{total}] cached  {rec:13s} {lid:14s} {aff}")
            continue
        v = subprocess.run(["vina", "--config", conf, "--ligand", ligpdbqt,
                            "--out", out, "--seed", "1"], capture_output=True, text=True)
        open(log, "w").write(v.stdout)
        aff = best_affinity(v.stdout)
        results.append({"ligand_id": lid, "name": lig["name"], "status": lig["status"],
                        "receptor": rec, "affinity": aff})
        el = time.time() - t0
        print(f"[{done}/{total}] {rec:13s} {lid:14s} {aff} kcal/mol  ({el:.0f}s elapsed)")

out_csv = os.path.join(DOCK, "docking_results.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand_id", "name", "status", "receptor", "affinity"])
    w.writeheader(); w.writerows(results)
print(f"\nDocked {total} pairs in {time.time()-t0:.0f}s -> {out_csv}")
