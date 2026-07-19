"""
Build gcp/manifest.csv = every (receptor x prepared-candidate) docking pair.

Enumerates the candidate .pdbqt files actually present in 02_ligands/candidates/
(so it only lists ligands that were successfully prepared), crossed with the
three receptors. Run locally BEFORE staging inputs to the bucket.

Uses only the standard library, so it runs under a bare python3 (e.g. in Cloud
Shell) as well as under the project venv.

Run:  python3 gcp/make_manifest.py
"""
import csv, glob, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CANDIR = os.path.join(BASE, "02_ligands", "candidates")
OUT = os.path.join(BASE, "gcp", "manifest.csv")
RECEPTORS = ["hiv_protease", "cathepsin_d", "bace1"]

ligands = sorted(os.path.splitext(os.path.basename(p))[0]
                 for p in glob.glob(os.path.join(CANDIR, "*.pdbqt")))
if not ligands:
    raise SystemExit(f"No candidate .pdbqt files in {CANDIR} — run step 12 first.")

pairs = [{"receptor": rec, "ligand_id": lid} for rec in RECEPTORS for lid in ligands]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "ligand_id"])
    w.writeheader(); w.writerows(pairs)

print(f"{len(ligands)} ligands x {len(RECEPTORS)} receptors = {len(pairs)} pairs -> {OUT}")
