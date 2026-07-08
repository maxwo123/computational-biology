"""
Steps 10 + 11 for receptors: protonate at physiological pH (PROPKA pKa assignment,
which sets the catalytic aspartate dyad ionization) and convert to PDBQT for Vina.

Pipeline per receptor:  *_clean.pdb --pdb2pqr--> *.pqr --meeko--> *.pdbqt (+ Vina box)

Run with the venv python:
    .venv/bin/python 04_scripts/02_receptors_to_pdbqt.py
"""
import csv, os, subprocess, sys

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
PROT = os.path.join(BASE, "01_proteins")
VENV = os.path.join(BASE, ".venv", "bin")

with open(os.path.join(BASE, "03_docking", "grid_boxes.csv")) as f:
    boxes = list(csv.DictReader(f))

for b in boxes:
    tag = b["tag"]
    clean = os.path.join(PROT, f"{tag}_clean.pdb")
    pqr   = os.path.join(PROT, f"{tag}.pqr")
    print(f"\n=== {tag} ===")

    # Step 10: pH-7 protonation + PROPKA pKa (protonates catalytic Asp dyad correctly).
    # Emit BOTH a PQR (charges/pKa record) and a protonated PDB (input for typing).
    pdb_h = os.path.join(PROT, f"{tag}_H.pdb")
    subprocess.run([
        os.path.join(VENV, "pdb2pqr"),
        "--ff=AMBER", "--with-ph=7.0",
        "--titration-state-method=propka", "--keep-chain",
        "--pdb-output", pdb_h,
        clean, pqr,
    ], check=True, capture_output=True, text=True)
    print(f"  protonated -> {os.path.basename(pdb_h)}")

    # Step 11: rigid receptor PDBQT via Open Babel (-xr) with Gasteiger charges.
    # NOTE: read the protonated *PDB*, not the PQR — reading PQR made Open Babel emit
    # nonsensical charges (|q|~100) that overflow the fixed-width column and corrupt
    # the atom-type field, breaking Vina's parser. Uniform route across all three
    # receptors and, unlike Meeko's template engine, handles BACE1's disulfide bonds.
    pdbqt = os.path.join(PROT, f"{tag}.pdbqt")
    r = subprocess.run([
        os.path.join(VENV, "obabel"), pdb_h, "-O", pdbqt, "-xr", "--partialcharge", "gasteiger",
    ], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(pdbqt):
        print(f"  OBABEL FAILED:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        sys.exit(1)
    n_atoms = sum(1 for l in open(pdbqt) if l.startswith(("ATOM", "HETATM")))
    print(f"  PDBQT -> {os.path.basename(pdbqt)} ({n_atoms} atoms)")

print("\nAll receptors prepared.")
