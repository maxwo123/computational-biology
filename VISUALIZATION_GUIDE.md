# Laboratory Manual — Visualizing the Docking Results

**Module goal.** By the end of this module you will be able to open your protein
structures and docked drug poses in several different visualization programs,
understand *what you are actually looking at* on the screen, and read a docked
pose the way a medicinal chemist does — identifying the binding pocket, the drug,
and the specific amino-acid contacts that explain **selectivity**.

This guide assumes **no prior experience** with molecular graphics software. Every
technical term is defined the first time it appears, and every instruction explains
*why* you are doing it, not just *what* to click.

---

## Part 0 — Essential concepts before you open anything

Molecular visualization is just a way of drawing the 3-D coordinates of atoms.
Before touching the software, you need to understand the three *kinds* of files
in this project, because opening the wrong one is the single most common source of
confusion.

### 0.1 The three file types you will visualize

| File type | Extension | What it is, biochemically | Example in this project |
|-----------|-----------|---------------------------|-------------------------|
| **Coordinate file** | `.pdb`, `.cif` | A list of every atom in a molecule and its (x, y, z) position in space, measured in Ångströms (Å; 1 Å = 10⁻¹⁰ m, roughly the width of one atom). This is the raw 3-D structure of a protein or a drug. | `01_proteins/1lyb.pdb` (Cathepsin D crystal structure) |
| **Docking input file** | `.pdbqt` | The same coordinates, but re-formatted for the docking program AutoDock Vina. The extra letters **Q** and **T** stand for the partial atomic **Q**charge and the **T**orsion (rotatable-bond) information Vina needs to flex the drug. Think of it as a "coordinate file with docking metadata." | `01_proteins/hiv_protease.pdbqt` (the receptor); `02_ligands/library/darunavir.pdbqt` (a drug) |
| **Session file** | `.pse` | A saved *snapshot of the PyMOL program itself* — which molecules were loaded, how they were colored, the camera angle. It is **not** raw data; it only remembers whatever was loaded into PyMOL at the moment you pressed Save. | `7:2:2026.pse` |

> **⚠️ Critical point about the session file.**
> Opening `7:2:2026.pse` shows only the structures that were loaded during **Phase 1**
> of the project: the three HIV-1 protease crystal structures (1HXW, 2IEN, 3QOZ) and
> the pocket-center markers. It does **not** contain the off-target enzymes, the
> prepared drug library, or any of the 60 docked poses — those were generated *later*
> by automated scripts and were never saved back into the session. **To see the
> docking results, you must load the individual `.pdbqt` files directly** (Parts 1–2
> below). The `.pse` is your Phase-1 lab notebook page, not the finished experiment.

### 0.2 Vocabulary you will meet on screen

- **Receptor** — the target protein the drug binds to (here: HIV-1 protease, or one of the two human off-target enzymes). In docking it is held **rigid**.
- **Ligand** — the small molecule (a drug) that binds. In docking it is **flexible** — Vina rotates its bonds to find the best fit.
- **Binding pocket (active site)** — the concave cavity on the protein surface where the ligand sits. For all three enzymes here, this is the catalytic cleft containing the aspartic-acid (Asp) dyad that would normally cleave a peptide bond.
- **Pose** — one specific 3-D arrangement of the ligand inside the pocket that Vina proposes. Vina returns up to 9 poses per drug, ranked best-first by predicted binding energy (ΔG, in kcal/mol; more negative = tighter, more favorable binding).
- **State / model** — when a single `.pdbqt` pose file contains all 9 poses stacked together, PyMOL calls each one a "state" and Mol* calls it a "model." **State 1 = the best (lowest-energy) pose** — this is the one you almost always want.
- **Representation** — the drawing style: *cartoon* (ribbon showing protein fold), *sticks* (bonds drawn as rods — good for drugs and pocket residues), *surface* (solid molecular skin — good for seeing pocket shape), *spheres* (atoms as balls).

---

## Part 1 — PyMOL (primary tool)

**What PyMOL is.** A molecular graphics program for displaying and rendering 3-D
structures. It is the standard tool in structural biology for making publication
figures. On this computer it is installed at
`/Applications/PyMOL.app` — launch it like any other Mac application.

There are two ways to work: an **automated loader** we provide (fastest), and the
**manual method** (teaches you the underlying steps). Do the automated method first
to see a result, then repeat manually to learn the mechanics.

### 1A. Automated method — the provided loader script

We wrote a helper script, `04_scripts/08_pymol_view.py`, that loads a receptor and
a docked drug pose together and automatically highlights the binding pocket. This
saves you ~15 manual clicks per pose.

**Procedure:**

1. **Launch PyMOL.** You will see a black 3-D viewport and, below it, a
   command-line box (the small white text field labeled `PyMOL>`). PyMOL is driven
   both by mouse and by typed commands; we will use both.
2. From the top menu choose **`File > Run Script…`** and select
   `04_scripts/08_pymol_view.py`.
   *(Equivalently, type `run 04_scripts/08_pymol_view.py` in the command line.)*
3. The script immediately builds a **demonstration scene**: the drug **darunavir**
   docked inside **HIV-1 protease**. What you are seeing:
   - **Gray ribbon (cartoon)** = the protein backbone fold of the receptor, made semi-transparent so you can see inside the pocket.
   - **Yellow sticks** = darunavir, in its best docked pose.
   - **Cyan sticks** = the amino-acid side chains of the receptor that lie within 5 Å of the drug — i.e., the residues physically lining the binding pocket and making contact with the drug. *These cyan residues are the heart of the selectivity question: they are the protein "hand" gripping the drug.*
4. **To view any other drug–receptor combination,** type a command like the
   following into the `PyMOL>` box and press Enter:
   ```
   show_pose 'cathepsin_d', 'darunavir'
   show_pose 'bace1', 'tipranavir'
   show_pose 'hiv_protease', 'saquinavir'
   ```
   - The first argument is the **receptor**: `hiv_protease`, `cathepsin_d`, or `bace1`.
   - The second is the **drug** (any `ligand_id` listed in `02_ligands/library.csv`, e.g. `ritonavir`, `indinavir`, `lopinavir`, `dmp323`, `mozenavir`, …).
   - *Biochemical purpose:* by loading the **same drug** into the viral target and then into a human off-target, you can visually compare how well it fits each pocket — the structural basis of selectivity.
5. **To visualize the method-validation experiment** (the re-docking control from
   Steps 13–14), type:
   ```
   show_redock 'hiv_protease', 'ritonavir'
   ```
   This overlays two copies of ritonavir in the same pocket:
   - **Magenta sticks** = the *experimental* position of ritonavir, taken directly from the X-ray crystal structure (ground truth).
   - **Green sticks** = the position our docking protocol *predicted*.
   Because these two nearly overlap (they differ by only 1.93 Å — less than the
   width of a single chemical bond), you are looking at visual proof that the
   docking protocol reproduces reality. This is the structural-biology equivalent
   of a positive control in a wet-lab assay.

### 1B. Manual method — loading files yourself

Do this once to understand what the script automates.

1. **Load the receptor.** `File > Open…` → navigate to
   `01_proteins/` → choose `hiv_protease.pdbqt`. A tangle of lines appears; this is
   the protein drawn in "lines" representation.
2. **Make it readable.** In the right-hand object panel, find the row named
   `hiv_protease`. Click the **`H`** (hide) button → `everything`, then the **`S`**
   (show) button → `cartoon`. You now see the ribbon fold instead of every atom —
   far easier to interpret. *Structurally, the cartoon traces the path of the
   protein's peptide backbone, with helices drawn as coils and β-sheets as arrows.*
3. **Load the docked drug.** `File > Open…` →
   `03_docking/poses/` → choose `hiv_protease__darunavir.pdbqt`. (Note the **double
   underscore** in the filename: `receptor__ligand`.)
4. **Show the drug as sticks.** Click the pose's row in the object panel →
   `S > sticks`. The drug is now drawn as rods you can see clearly.
5. **Step through the 9 poses.** Look at the bottom-right of the window for the
   playback controls (◀ ▶). Each frame is one Vina pose; **frame/state 1 is the
   best-scoring pose.** Advancing the frames shows you the alternative binding
   modes Vina considered.
6. **Reveal the pocket.** In the command line, type:
   ```
   select pocket, hiv_protease within 5 of hiv_protease__darunavir
   show sticks, pocket
   ```
   The first line *selects* every receptor atom within 5 Å of the drug (the
   contact residues); the second draws them. You have now manually reconstructed
   what the loader script does automatically.

### 1C. Opening the original Phase-1 session

To revisit your early work, simply double-click `7:2:2026.pse` (or `File > Open`).
Remember from Part 0 that this shows only the three HIV protease crystals and the
pocket markers — it is a Phase-1 snapshot, not the docking results.

### 1D. Mouse controls (all viewers are similar)

- **Left-drag** = rotate the molecule.
- **Scroll wheel** (or right-drag up/down) = zoom.
- **Middle-drag** = pan (slide the view sideways).
- **`orient`** typed in the command line = auto-frame the current selection.

### 1E. Saving a figure

`File > Export Image As > PNG…` → choose **ray-traced** for a smooth, publication-
quality render. Save figures into `05_results/` to keep your lab notebook tidy.

---

## Part 2 — Mol* (web browser, no installation)

**What Mol\* is.** A free molecular viewer that runs entirely inside a web browser —
the same engine that powers the official Protein Data Bank website. Use it when you
want a quick look without opening PyMOL, or on a computer where PyMOL is not
installed.

**Procedure:**

1. In any web browser, go to **https://molstar.org/viewer**.
2. Click the **folder / "Open Files"** icon in the left toolbar.
3. Drag in a receptor file, e.g. `01_proteins/hiv_protease.pdbqt`. The protein
   renders automatically.
4. Drag in a pose file, e.g. `03_docking/poses/hiv_protease__darunavir.pdbqt`.
   The drug appears inside the pocket.
5. Rotate (left-drag), zoom (scroll), and right-click any residue for a menu of
   options such as "focus" or "label."

> **If Mol\* refuses a `.pdbqt` file:** some web viewers prefer the simpler `.pdb`
> or `.sdf` formats. Ask the instructor (or re-run the export step) to convert any
> pose to `.pdb`; those load in every viewer.

---

## Part 3 — py3Dmol (interactive viewer embedded in a web page)

**What py3Dmol is.** A small Python library that writes a *self-contained* HTML file
containing an interactive 3-D viewer. The advantage: once generated, the HTML file
opens in any browser with **no software at all** — ideal for sharing a rotatable
pose with an instructor or embedding in a report.

This tool is **optional** and is not yet set up. If you want it, it requires a
one-time installation into the project's Python environment:

```
.venv/bin/python -m pip install py3Dmol
```

A short generator script would then produce one `.html` file per pose in
`05_results/`. *(The project maintainer can add `04_scripts/09_make_html_viewers.py`
on request.)*

---

## Part 4 — PLIP: 2-D protein–ligand interaction diagrams

**What PLIP is, and why it matters most for this project.** The *Protein–Ligand
Interaction Profiler* analyzes a docked complex and produces a flat, labeled diagram
naming **every specific chemical interaction** between the drug and the pocket:
hydrogen bonds, hydrophobic contacts, salt bridges, π-stacking, and so on.

This is arguably the **most important visualization for the selectivity analysis**
(project Step 21), because selectivity is ultimately decided by *which* pocket
residues a drug can and cannot reach. A 2-D interaction map answers the question
"*why* does darunavir bind HIV-1 protease more tightly than BACE1?" in a way a 3-D
picture cannot, by listing the exact residue contacts — including whether the drug
engages the **S1 subpocket** highlighted in the project brief.

**Two ways to use it:**

- **Web server (no install):** go to **https://plip-tool.biotec.tu-dresden.de**,
  upload a single file that contains *both* the receptor and the docked drug merged
  together, and download the annotated diagram.
- **Local install:** PLIP can be added to the project's Python environment
  (`.venv/bin/python -m pip install plip`).

Either way, PLIP needs the receptor and the chosen pose combined into **one** `.pdb`
"complex" file first. That merge is one short command the maintainer can add as
`04_scripts/10_make_complex.py` — ask if you would like PLIP wired up.

---

## Part 5 — RDKit: a 2-D grid of the whole drug library

**What this shows.** Everything above is 3-D. But sometimes you simply want to see
the **flat chemical structures** of all 20 inhibitors at once — the way they are
drawn in an organic chemistry textbook (Chem 51 style) — to compare their chemical
scaffolds. RDKit (the cheminformatics library already installed for this project)
can render all 20 molecules from their SMILES strings into a single labeled image,
colored by clinical status (approved / discontinued / investigational).

This produces `05_results/library_grid.png`. It is one short script
(`04_scripts/11_library_grid.py`) the maintainer can add on request.

---

## Part 6 — The quantitative result figures (already generated)

These are ordinary images — open them with **Preview** or any image viewer by
double-clicking. They live in `05_results/`:

| File | What it shows | How to read it |
|------|---------------|----------------|
| `fig1_affinity_heatmap.png` | Predicted binding energy (ΔG) of every drug against every receptor | Darker = tighter binding. Compare a drug's HIV column to its off-target columns. |
| `fig2_selectivity.png` | Selectivity score (ΔΔG) per drug, colored by clinical status | Bars to the **right** = selective for HIV protease (safer); bars to the **left** = binds a human off-target as well or better. |
| `fig3_clusters_pca.png` | Compounds grouped by their overall 3-receptor binding "fingerprint" | Drugs near each other behave similarly across all three enzymes. |
| `fig4_selectivity_landscape.png` | HIV affinity vs. best off-target affinity | Points **below the diagonal** are HIV-selective. |

---

## Part 7 — What to actually *look for* (interpreting a pose)

Loading a pose is only step one; reading it is the real skill. When you have a drug
in a pocket on screen, work through this checklist:

1. **Is the drug actually inside the cleft?** A correct pose sits buried in the
   concave active site, not perched on the protein's outer surface. Switch the
   receptor to `surface` representation (PyMOL: `S > surface`) to confirm the drug
   is enclosed.
2. **Find the catalytic aspartate dyad.** All three enzymes cleave peptides using
   two Asp residues. A genuine protease inhibitor almost always makes contact with
   this dyad. In PyMOL: `select cat, receptor and resn ASP and pocket`.
3. **Count the hydrogen bonds.** These are the specific, directional "grip" points.
   More and better-aligned H-bonds generally mean tighter, more specific binding.
   (PLIP, Part 4, quantifies these precisely.)
4. **Compare viral vs. human pockets.** Load the same drug into `hiv_protease` and
   then `cathepsin_d`. If the drug fills the HIV pocket snugly but leaves gaps or
   clashes in the human pocket, that mismatch *is* the structural explanation for
   selectivity — the central finding of this project.
5. **Connect back to the numbers.** Cross-reference what you see with
   `05_results/selectivity_table.csv`: does a drug you *see* fitting HIV better also
   have a positive ΔΔG selectivity score? Visual and quantitative evidence should
   agree.

---

## Quick reference — file locations

```
Receptors (rigid targets):     01_proteins/{hiv_protease,cathepsin_d,bace1}.pdbqt
Prepared drugs (flexible):     02_ligands/library/<ligand_id>.pdbqt
Docked poses (the results):    03_docking/poses/<receptor>__<ligand_id>.pdbqt
Re-docking validation:         03_docking/redock_controls/
Result charts:                 05_results/*.png
PyMOL loader script:           04_scripts/08_pymol_view.py
Phase-1 session snapshot:      7:2:2026.pse
```

Valid receptor names: `hiv_protease`, `cathepsin_d`, `bace1`.
Valid `ligand_id` values: see the first column of `02_ligands/library.csv`.
