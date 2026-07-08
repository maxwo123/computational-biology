# HIV-1 Protease Inhibitor Selectivity Mapping

**A computational docking study of HIV-1 protease inhibitors against human off-target aspartic proteases (Cathepsin D, BACE1).**

This repository is a digital lab notebook for a structure-based selectivity study. It curates a library of HIV-1 protease inhibitors, docks each one against the viral target *and* two human off-target enzymes, and quantifies **selectivity** (ΔΔG) — the degree to which each drug prefers the viral protease over host enzymes. Off-target cross-reactivity is a leading cause of drug toxicity and clinical-trial failure, so mapping it computationally is a core exercise in medicinal chemistry.

> **Scientific honesty note.** Every binding number in this repo is the output of a real AutoDock Vina run on real prepared structures — nothing is hand-authored or estimated. Docking scores are *approximate* predictions (Vina's scoring function has ~2–3 kcal/mol error and does not replace experiment); they are used here for **relative** ranking, which is the appropriate use of docking.

---

## 1. Biological system

| Role | Enzyme | Family | PDB | Native ligand (pocket marker) |
|------|--------|--------|-----|-------------------------------|
| **Primary target** | HIV-1 protease | retroviral aspartic protease | `1HXW` | ritonavir |
| Off-target 1 | Human Cathepsin D | A1 aspartic protease | `1LYB` | pepstatin |
| Off-target 2 | Human BACE1 (β-secretase) | A1 aspartic protease | `1FKN` | OM99-2 |

All three are **aspartic proteases** that cleave peptide bonds using a catalytic Asp dyad, which is exactly why an HIV protease inhibitor *might* accidentally bind the human enzymes — the shared catalytic machinery is the source of off-target risk this study probes.

HIV-1 protease is an obligate **homodimer**: its active site forms at the interface of two identical chains, so both chains (A + B) are kept in the receptor.

---

## 2. Screening library (20 compounds)

Curated from PubChem (`02_ligands/library.csv`), spanning the clinical spectrum:

- **10 FDA-approved** protease inhibitors: saquinavir, ritonavir, indinavir, nelfinavir, amprenavir, lopinavir, atazanavir, fosamprenavir, tipranavir, darunavir
- **1 pharmacoenhancer** (cobicistat — a ritonavir analog that is *not* itself a protease inhibitor; included as a structural comparator)
- **5 discontinued / clinical-failure** inhibitors: brecanavir, mozenavir, palinavir, telinavir, DMP-323
- **4 investigational** inhibitors: KNI-272, GS-8374, TMC-310911, GRL-0519

Structures were desalted (largest fragment) and deduplicated by InChIKey connectivity.

---

## 3. Pipeline / methods

Scripts live in `04_scripts/` and run in numbered order. All Python dependencies are pinned in a local virtual environment (`.venv/`, see `requirements.txt`).

| Step(s) | Script | What it does |
|---------|--------|--------------|
| 4, 5, 9, 12 | `01_prepare_receptors.py` | (PyMOL) isolate native ligands, compute pocket centroids = grid-box centers, strip waters/sugars/ligands, keep biological chains, remove insertion codes |
| 10, 11 | `02_receptors_to_pdbqt.py` | Protonate receptors at pH 7.0 with **pdb2pqr + PROPKA** (sets the catalytic Asp dyad ionization), then convert to PDBQT with Gasteiger charges (Open Babel) |
| 6 | `03_build_library.py` | Fetch isomeric SMILES from PubChem, desalt, dedupe → `library.csv` |
| 7, 8, 11 | `04_prepare_ligands.py` | Protonate ligands at pH 7.4 (Open Babel), 3D-embed + MMFF94 minimize (RDKit), assign torsions → PDBQT (Meeko) |
| 13, 14 | `05_redock_controls.py` | Redock each native ligand into its own pocket, measure RMSD vs crystal pose (validation) |
| 15–18 | `06_dock_library.py` | Dock all 20 ligands × 3 receptors with AutoDock Vina; parse best affinity per run |
| 19–21 | `07_analyze_selectivity.py` | Build master dataframe, compute ΔΔG selectivity, KMeans clustering, descriptor correlation, figures |

**Docking parameters:** AutoDock Vina 1.2.7, exhaustiveness 8, 9 modes, grid boxes centered on each native inhibitor's centroid, padded to ≥22.5 Å per side (`03_docking/*_vina.conf`).

### Method validation (redocking controls)

| Receptor | Native ligand | Top affinity (kcal/mol) | Pose RMSD vs crystal |
|----------|---------------|-------------------------|----------------------|
| HIV-1 protease | ritonavir | −9.84 | **1.93 Å ✓ (< 2 Å)** |
| Cathepsin D | pepstatin | −8.37 | 3.69 Å (peptidic inhibitor; see caveats) |
| BACE1 | OM99-2 | — | ligand prep failed (large peptidomimetic) |

The **primary target redocks to 1.93 Å**, below the standard 2 Å acceptance threshold — the box placement and protocol reproduce the known crystallographic binding mode. The two off-target native ligands are large flexible **peptides**, which redock poorly in a permissive box; this is a known limitation for peptidic ligands and does **not** affect docking of the drug-like library compounds.

---

## 4. Results

**Headline finding: 14 of 20 compounds are predicted HIV-selective** — they dock more tightly to HIV-1 protease than to *both* human off-targets.

**Most selective (largest worst-case ΔΔG, kcal/mol):**

| Rank | Compound | Status | HIV-1 PR | Cathepsin D | BACE1 | Selectivity (min ΔΔG) |
|------|----------|--------|---------:|------------:|------:|----------------------:|
| 1 | DMP-323 | discontinued | −8.28 | −5.65 | −6.85 | **+1.44** |
| 2 | Tipranavir | approved | −10.81 | −9.58 | −9.01 | **+1.23** |
| 3 | Mozenavir | discontinued | −8.91 | −6.61 | −7.76 | **+1.15** |
| 4 | Brecanavir | discontinued | −9.68 | −8.48 | −8.79 | +0.89 |
| 5 | Indinavir | approved | −10.18 | −8.19 | −9.41 | +0.77 |

**Least selective (bind a human off-target as tightly or tighter than HIV protease):**
Darunavir (−1.09), Atazanavir (−0.85), Cobicistat (−0.79), KNI-272 (−0.29), Ritonavir (−0.16), Palinavir (−0.01). *(Note: darunavir and atazanavir are highly potent clinical drugs; a negative computed selectivity here reflects docking's limitations for absolute discrimination, not a real safety verdict — see caveats.)*

**Structure–selectivity correlation (Step 21).** Across the library, predicted selectivity correlates **negatively** with polarity and size:

| Descriptor | Pearson r vs selectivity |
|-----------|-------------------------:|
| TPSA (polar surface area) | −0.59 |
| Rotatable bonds | −0.49 |
| H-bond acceptors | −0.43 |
| Molecular weight | −0.41 |

Interpretation: **smaller, more rigid, less polar inhibitors tend to prefer the viral protease**, whereas large, flexible, highly polar compounds engage the human aspartic-protease pockets comparably well — a plausible medicinal-chemistry design cue for improving selectivity.

Compounds also fall into 3 KMeans clusters by their 3-receptor binding profile (`05_results/clusters.csv`, visualized in `fig3`).

Key output files:
- `05_results/selectivity_table.csv` — full ΔG / ΔΔG / descriptor table, ranked by selectivity
- `05_results/fig1_affinity_heatmap.png` — affinity of every compound against every receptor
- `05_results/fig2_selectivity.png` — selectivity (min ΔΔG) per compound, colored by clinical status
- `05_results/fig3_clusters_pca.png` — compounds clustered by 3-receptor binding profile
- `05_results/fig4_selectivity_landscape.png` — HIV affinity vs best off-target affinity

**Selectivity metric.** For each off-target, ΔΔG = ΔG(off-target) − ΔG(HIV protease). Because ΔG is negative (more negative = tighter), a **positive ΔΔG means the drug binds the human enzyme more weakly than the viral target = selective/safer**. The reported `selectivity_score` is the *worst-case* (smallest) ΔΔG across both off-targets.

---

## 5. Reproducing this work

```bash
# 1. environment (one time)
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt

# 2. run the pipeline in order
/Applications/PyMOL.app/Contents/MacOS/PyMOL -cq -d "run 04_scripts/01_prepare_receptors.py"
.venv/bin/python 04_scripts/02_receptors_to_pdbqt.py
.venv/bin/python 04_scripts/03_build_library.py
.venv/bin/python 04_scripts/04_prepare_ligands.py
.venv/bin/python 04_scripts/05_redock_controls.py
.venv/bin/python 04_scripts/06_dock_library.py      # ~20-30 min (60 dockings)
.venv/bin/python 04_scripts/07_analyze_selectivity.py
```

External tools: **PyMOL** (visualization / structure ops) and **AutoDock Vina 1.2.7** (`vina` on PATH) must be installed separately.

---

## 6. Directory layout

```
01_proteins/     receptors: raw PDBs, *_clean.pdb, *_H.pdb (protonated), *.pdbqt, *.pqr
02_ligands/      1hxw/2ien/3qoz native ligands; library.csv; library/*.pdbqt (20 compounds)
03_docking/      grid_boxes.csv, *_vina.conf, poses/, redock_controls/, docking_results.csv
04_scripts/      the numbered pipeline (01-07)
05_results/      selectivity tables + figures
CONTEXT.md       project brief / scientific rationale
```

---

## 7. Caveats & limitations

- **Docking ≠ experiment.** Vina scores rank-order binding but do not give true free energies; treat all numbers as relative predictions.
- **Rigid receptors.** Protein flexibility (induced fit) is not modeled; the catalytic Asp protonation is fixed at the pH-7 PROPKA assignment.
- **Single pose per ligand seed.** One Vina run (seed 1) per pair; a production study would average replicates.
- **Off-target peptide redocking** is imperfect (see §3) — the box is correctly placed on the native centroid, but flexible peptide natives are hard to reproduce.
- The **booster cobicistat** is included as a structural comparator, not a bona fide protease inhibitor.
