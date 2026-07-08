"""
PyMOL visualization helper for the docking results.

The project's 7:2:2026.pse session only holds the Phase-1 crystal structures.
The DOCKING results live in separate files, so this script loads them into a
nice scene for you.

HOW TO USE (in the PyMOL GUI):
    1. Open PyMOL
    2. File > Run Script...  ->  select this file
       (or in the PyMOL command line:  run 04_scripts/08_pymol_view.py )
    3. It loads a demo scene: darunavir docked in HIV-1 protease.
    4. To view any other pair, type in the PyMOL command line, e.g.:
           show_pose 'cathepsin_d', 'darunavir'
           show_pose 'bace1', 'tipranavir'
       To overlay a redocked native on its crystal pose (validation):
           show_redock 'hiv_protease', 'ritonavir'

Receptors:  hiv_protease | cathepsin_d | bace1
Ligands  :  any ligand_id in 02_ligands/library.csv (saquinavir, ritonavir, ...)
"""
import os
from pymol import cmd

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"

def _style_receptor(obj):
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.color("gray80", obj)
    cmd.set("cartoon_transparency", 0.3, obj)

def show_pose(receptor="hiv_protease", ligand="darunavir", pose=1):
    """Load a receptor + the top docked pose of a ligand, and frame the pocket."""
    cmd.delete("all")
    rec = os.path.join(BASE, "01_proteins", f"{receptor}.pdbqt")
    pos = os.path.join(BASE, "03_docking", "poses", f"{receptor}__{ligand}.pdbqt")
    if not os.path.exists(pos):
        print(f"!! no pose file for {receptor} + {ligand}: {pos}")
        return
    cmd.load(rec, "receptor")
    cmd.load(pos, "ligand")
    cmd.split_states("ligand")          # poses are multi-model; keep the best (state/model `pose`)
    cmd.delete("ligand")
    best = f"ligand_{pose:04d}"
    cmd.set_name(best, "pose")
    cmd.delete("ligand_*")
    _style_receptor("receptor")
    cmd.show("sticks", "pose")
    cmd.color("yellow", "pose and elem C")
    # highlight pocket residues within 5 A of the ligand
    cmd.select("pocket", f"receptor within 5 of pose")
    cmd.show("sticks", "pocket")
    cmd.color("cyan", "pocket and elem C")
    cmd.set("cartoon_transparency", 0.5, "receptor")
    cmd.orient("pose")
    cmd.zoom("pose", 8)
    print(f"Showing {ligand} docked in {receptor} (pose {pose}). Pocket residues within 5A shown in cyan.")

def show_redock(receptor="hiv_protease", native="ritonavir"):
    """Overlay a redocked native ligand (green) on its crystal pose (magenta)."""
    cmd.delete("all")
    rec = os.path.join(BASE, "01_proteins", f"{receptor}.pdbqt")
    redock = os.path.join(BASE, "03_docking", "redock_controls", f"{receptor}_{native}_pose1.sdf")
    crystal = os.path.join(BASE, "03_docking", "redock_controls", f"{receptor}_{native}_crystal.sdf")
    cmd.load(rec, "receptor"); _style_receptor("receptor")
    if os.path.exists(crystal):
        cmd.load(crystal, "crystal"); cmd.show("sticks", "crystal"); cmd.color("magenta", "crystal and elem C")
    if os.path.exists(redock):
        cmd.load(redock, "redocked"); cmd.show("sticks", "redocked"); cmd.color("green", "redocked and elem C")
    cmd.orient("crystal or redocked"); cmd.zoom("crystal or redocked", 6)
    print(f"Crystal (magenta) vs redocked (green) for {native} in {receptor}.")

# register as PyMOL commands so you can type them in the command line
cmd.extend("show_pose", show_pose)
cmd.extend("show_redock", show_redock)

# load a default demo scene on run
show_pose("hiv_protease", "darunavir")
print("\nType e.g.  show_pose 'bace1', 'tipranavir'   or   show_redock 'hiv_protease', 'ritonavir'")
