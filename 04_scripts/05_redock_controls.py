"""
Steps 13-14: native re-docking controls + RMSD validation.

Redock each receptor's own co-crystallized inhibitor back into its pocket, then
measure heavy-atom RMSD between the top Vina pose and the crystallographic pose.
RMSD < 2.0 A is the accepted threshold that the box + protocol reproduce reality.

Run:  .venv/bin/python 04_scripts/05_redock_controls.py
"""
import csv, os, subprocess

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
VENV = os.path.join(BASE, ".venv", "bin")
OBABEL = os.path.join(VENV, "obabel")
_OB = os.path.join(BASE, ".venv/lib/python3.9/site-packages/openbabel")
OBRMS = os.path.join(_OB, "bin", "obrms")
# the raw obrms binary (unlike the obabel wrapper) needs plugin dirs pointed out
OBRMS_ENV = dict(os.environ,
                 BABEL_LIBDIR=os.path.join(_OB, "lib", "openbabel", "3.1.0"),
                 BABEL_DATADIR=os.path.join(_OB, "share", "openbabel", "3.1.0"))
DOCK = os.path.join(BASE, "03_docking")
CTRL = os.path.join(DOCK, "redock_controls")
os.makedirs(CTRL, exist_ok=True)

# (receptor_tag, native_name, crystal_ligand_pdb, ligand_pdbqt_to_dock_or_None)
# For ritonavir we reuse the clean library PDBQT; Vina re-samples the pose anyway.
CONTROLS = [
    ("hiv_protease", "ritonavir",
     os.path.join(BASE, "02_ligands", "1hxw_ligand.pdb"),
     os.path.join(BASE, "02_ligands", "library", "ritonavir.pdbqt")),
    ("cathepsin_d", "pepstatin",
     os.path.join(BASE, "02_ligands", "cathepsin_d_native_pepstatin.pdb"), None),
    ("bace1", "om99-2",
     os.path.join(BASE, "02_ligands", "bace1_native_om99-2.pdb"), None),
]

def prep_ligand_from_pdb(pdb, out_pdbqt):
    """Protonate a crystal ligand at pH 7.4 and write a Vina ligand PDBQT (Open Babel)."""
    subprocess.run([OBABEL, pdb, "-O", out_pdbqt, "-p", "7.4", "--partialcharge", "gasteiger"],
                   capture_output=True, text=True)
    return os.path.exists(out_pdbqt)

results = []
for tag, name, crystal_pdb, lig_pdbqt in CONTROLS:
    print(f"\n=== redock {name} into {tag} ===")
    if lig_pdbqt is None:
        lig_pdbqt = os.path.join(CTRL, f"{name}.pdbqt")
        if not prep_ligand_from_pdb(crystal_pdb, lig_pdbqt):
            print(f"  ligand prep failed for {name}")
            results.append({"receptor": tag, "native": name, "top_affinity": "", "rmsd": "prep_failed"})
            continue

    out = os.path.join(CTRL, f"{tag}_{name}_redock.pdbqt")
    v = subprocess.run(["vina", "--config", os.path.join(DOCK, f"{tag}_vina.conf"),
                        "--ligand", lig_pdbqt, "--out", out, "--seed", "1"],
                       capture_output=True, text=True)
    # parse best affinity (mode 1)
    aff = ""
    for line in v.stdout.splitlines():
        s = line.split()
        if len(s) >= 4 and s[0] == "1":
            aff = s[1]; break
    if not os.path.exists(out):
        print(f"  vina failed: {v.stdout[-300:]}{v.stderr[-300:]}")
        results.append({"receptor": tag, "native": name, "top_affinity": aff, "rmsd": "dock_failed"})
        continue

    # extract top pose (first model) and reference, then RMSD
    pose = os.path.join(CTRL, f"{tag}_{name}_pose1.sdf")
    ref = os.path.join(CTRL, f"{tag}_{name}_crystal.sdf")
    subprocess.run([OBABEL, out, "-O", pose, "-f", "1", "-l", "1"], capture_output=True, text=True)
    subprocess.run([OBABEL, crystal_pdb, "-O", ref], capture_output=True, text=True)
    rms = subprocess.run([OBRMS, ref, pose], capture_output=True, text=True, env=OBRMS_ENV)
    rmsd = ""
    if rms.stdout.strip():
        # obrms prints: "RMSD <file1>:<file2> <value>"
        try:
            rmsd = round(float(rms.stdout.strip().split()[-1]), 3)
        except ValueError:
            rmsd = rms.stdout.strip()
    flag = ""
    if isinstance(rmsd, float):
        flag = "  PASS (<2A)" if rmsd < 2.0 else "  (>2A)"
    print(f"  top affinity = {aff} kcal/mol   RMSD = {rmsd}{flag}")
    results.append({"receptor": tag, "native": name, "top_affinity": aff, "rmsd": rmsd})

with open(os.path.join(DOCK, "redock_validation.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["receptor", "native", "top_affinity", "rmsd"])
    w.writeheader(); w.writerows(results)
print(f"\nWrote {DOCK}/redock_validation.csv")
