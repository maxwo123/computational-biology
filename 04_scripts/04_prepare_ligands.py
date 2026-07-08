"""
Steps 7, 8, 11 for the screening library.

Per compound:
  Step 7  protonation  : dominant microspecies at pH 7.4 (Open Babel -p)
  Step 8  3D + minimize: RDKit ETKDGv3 embedding + MMFF94 energy minimization
  Step 11 PDBQT        : Meeko assigns AutoDock atom types + rotatable-bond torsion tree

Outputs one <ligand_id>.sdf (minimized 3D) and <ligand_id>.pdbqt per compound in
02_ligands/library/, plus a prep_report.csv.

Run:  .venv/bin/python 04_scripts/04_prepare_ligands.py
"""
import csv, os, subprocess, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
VENV = os.path.join(BASE, ".venv", "bin")
LIBDIR = os.path.join(BASE, "02_ligands", "library")
os.makedirs(LIBDIR, exist_ok=True)

def protonate_ph74(smiles):
    """Return the dominant protonation state at pH 7.4 (Open Babel)."""
    r = subprocess.run([os.path.join(VENV, "obabel"), f"-:{smiles}", "-osmi", "-p", "7.4"],
                       capture_output=True, text=True)
    out = r.stdout.strip().split("\t")[0].strip()
    return out or smiles

def embed_minimize(smiles):
    """SMILES -> 3D RDKit mol, MMFF94-minimized. Returns mol or None."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, "unparseable SMILES"
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3(); params.randomSeed = 0xf00d
    if AllChem.EmbedMolecule(m, params) != 0:
        # retry with random coordinates for awkward macrocycles
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(m, params) != 0:
            return None, "embedding failed"
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception as e:
        return None, f"MMFF failed: {e}"
    return m, "ok"

rows = list(csv.DictReader(open(os.path.join(BASE, "02_ligands", "library.csv"))))
report = []
for r in rows:
    lid, smi = r["ligand_id"], r["smiles"]
    prot = protonate_ph74(smi)
    mol, status = embed_minimize(prot)
    if mol is None:
        print(f"  FAIL  {lid:16s} {status}")
        report.append({"ligand_id": lid, "status": status, "protonated_smiles": prot, "n_torsions": ""})
        continue
    sdf = os.path.join(LIBDIR, f"{lid}.sdf")
    Chem.MolToMolFile(mol, sdf)
    pdbqt = os.path.join(LIBDIR, f"{lid}.pdbqt")
    mk = subprocess.run([os.path.join(VENV, "python"), os.path.join(VENV, "mk_prepare_ligand.py"),
                         "-i", sdf, "-o", pdbqt], capture_output=True, text=True)
    if not os.path.exists(pdbqt):
        print(f"  FAIL  {lid:16s} meeko: {mk.stdout[-200:]} {mk.stderr[-200:]}")
        report.append({"ligand_id": lid, "status": "meeko_failed", "protonated_smiles": prot, "n_torsions": ""})
        continue
    ntor = sum(1 for l in open(pdbqt) if l.startswith("BRANCH"))
    print(f"  OK    {lid:16s} torsions={ntor}")
    report.append({"ligand_id": lid, "status": "ok", "protonated_smiles": prot, "n_torsions": ntor})

with open(os.path.join(LIBDIR, "prep_report.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand_id", "status", "n_torsions", "protonated_smiles"])
    w.writeheader(); w.writerows(report)

ok = sum(1 for x in report if x["status"] == "ok")
print(f"\nPrepared {ok}/{len(rows)} ligands -> {LIBDIR}")
