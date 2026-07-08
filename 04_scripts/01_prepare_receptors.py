"""
Step 9 (receptor cleaning) + Step 12 (grid box centers) for all three receptors.

Run with PyMOL's python:
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq -d "run 04_scripts/01_prepare_receptors.py"

For each receptor we:
  - load the raw PDB
  - record the co-crystallized inhibitor's centroid = the binding-pocket center
    (this is the docking grid-box center used later by AutoDock Vina)
  - strip waters, the inhibitor, ions, and (for Cathepsin D) N-linked glycan sugars
  - keep only the biologically relevant protein chains, and save a clean receptor PDB
"""
from pymol import cmd
import os, csv

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
PROT = os.path.join(BASE, "01_proteins")
LIGD = os.path.join(BASE, "02_ligands")

# (tag, raw_structure_path, protein_chains_to_keep, inhibitor_selection, native_ligand_name)
# HIV-1 protease is an obligate HOMODIMER: its active site sits at the A/B interface,
# so BOTH chains A and B must be kept.
RECEPTORS = [
    ("hiv_protease", os.path.join(LIGD, "1hxw.cif"), "chain A+B", "resn RIT", "ritonavir"),
    ("cathepsin_d",  os.path.join(PROT, "1lyb.pdb"), "chain A+B", "chain I",  "pepstatin"),
    ("bace1",        os.path.join(PROT, "1fkn.pdb"), "chain A",   "chain C",  "om99-2"),
]

rows = []
for tag, raw_path, keep_chains, inhib_sel, lig_name in RECEPTORS:
    cmd.delete("all")
    cmd.load(raw_path, tag)

    # --- pocket center from the native inhibitor ---
    inhib = f"({tag}) and ({inhib_sel})"
    n_inhib = cmd.count_atoms(inhib)
    com = cmd.centerofmass(inhib)
    (minc, maxc) = cmd.get_extent(inhib)
    size = [maxc[i] - minc[i] for i in range(3)]
    print(f"[{tag}] inhibitor '{inhib_sel}' atoms={n_inhib}  center=({com[0]:.2f},{com[1]:.2f},{com[2]:.2f})")

    # save the off-target native ligand for redocking controls (HIV ritonavir already saved earlier)
    if tag != "hiv_protease":
        cmd.save(os.path.join(LIGD, f"{tag}_native_{lig_name}.pdb"), inhib)

    # --- clean receptor: protein of the kept chains only, no solvent/hetero ---
    clean_sel = f"({tag}) and ({keep_chains}) and polymer.protein and not hetatm"
    # (polymer.protein already excludes waters/sugars/inhibitor; the explicit filters are belt-and-suspenders)
    n_clean = cmd.count_atoms(clean_sel)
    out_pdb = os.path.join(PROT, f"{tag}_clean.pdb")
    cmd.create(f"{tag}_clean", clean_sel)

    # Remove insertion codes by renumbering residues sequentially. BACE1 (memapsin 2)
    # uses insertion codes (e.g. "43P") that downstream PQR parsing can't handle.
    # Docking boxes are defined by coordinates, not residue IDs, so this is harmless.
    _rn = {"n": 0, "last": None}
    def _renumber(chain, resi, segi):
        key = (segi, chain, resi)
        if key != _rn["last"]:
            _rn["n"] += 1
            _rn["last"] = key
        return str(_rn["n"])
    cmd.alter(f"{tag}_clean", "resi = _renumber(chain, resi, segi)", space={"_renumber": _renumber})
    cmd.alter(f"{tag}_clean", 'ins = ""')
    cmd.sort(f"{tag}_clean")
    cmd.save(out_pdb, f"{tag}_clean")
    print(f"[{tag}] clean receptor atoms={n_clean} -> {os.path.basename(out_pdb)}")

    # pad the ligand extent to a cubic box (min 22.5 A/side) for Vina search space
    box = [max(22.5, round(s + 10.0, 1)) for s in size]
    rows.append({
        "tag": tag, "native_ligand": lig_name,
        "center_x": round(com[0], 3), "center_y": round(com[1], 3), "center_z": round(com[2], 3),
        "size_x": box[0], "size_y": box[1], "size_z": box[2],
    })

# write consolidated grid-box config
out_csv = os.path.join(BASE, "03_docking", "grid_boxes.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"\nWrote {out_csv}")
for r in rows:
    print(r)
