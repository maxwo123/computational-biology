"""
Steps 19-21: selectivity analysis + clustering.

  Step 19  build a master dataframe of all docking scores (ligand x receptor)
  Step 20  selectivity metric  DDG = DG_offtarget - DG_primary   (kcal/mol)
           (DG is negative; a POSITIVE DDG => weaker off-target binding => selective)
  Step 21  cluster compounds by their 3-receptor binding profile and correlate
           selectivity with molecular descriptors

Outputs go to 05_results/ (CSV tables + PNG figures).

Run:  .venv/bin/python 04_scripts/07_analyze_selectivity.py
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
RES = os.path.join(BASE, "05_results")
os.makedirs(RES, exist_ok=True)

PRIMARY = "hiv_protease"
OFFTARGETS = ["cathepsin_d", "bace1"]

# ---- Step 19: master dataframe ----
dock = pd.read_csv(os.path.join(BASE, "03_docking", "docking_results.csv"))
lib = pd.read_csv(os.path.join(BASE, "02_ligands", "library.csv"))

wide = dock.pivot_table(index="ligand_id", columns="receptor", values="affinity")
wide = wide.merge(lib.set_index("ligand_id")[["name", "status", "smiles", "mw"]],
                  left_index=True, right_index=True)

# ---- Step 20: selectivity metrics ----
for off in OFFTARGETS:
    wide[f"ddg_{off}"] = wide[off] - wide[PRIMARY]            # >0 => selective for HIV
# worst-case selectivity = smallest gap to any off-target (the limiting off-target)
wide["selectivity_score"] = wide[[f"ddg_{o}" for o in OFFTARGETS]].min(axis=1)
wide["selective_for_hiv"] = wide["selectivity_score"] > 0

# ---- molecular descriptors for Step 21 correlation ----
def descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return pd.Series({"MW": np.nan, "logP": np.nan, "TPSA": np.nan,
                          "HBD": np.nan, "HBA": np.nan, "RotB": np.nan})
    return pd.Series({
        "MW": Descriptors.MolWt(m), "logP": Crippen.MolLogP(m),
        "TPSA": Descriptors.TPSA(m), "HBD": Descriptors.NumHDonors(m),
        "HBA": Descriptors.NumHAcceptors(m), "RotB": Descriptors.NumRotatableBonds(m),
    })
desc = wide["smiles"].apply(descriptors)
wide = pd.concat([wide, desc], axis=1)

cols = ["name", "status", PRIMARY] + OFFTARGETS + \
       [f"ddg_{o}" for o in OFFTARGETS] + ["selectivity_score", "selective_for_hiv",
        "MW", "logP", "TPSA", "HBD", "HBA", "RotB"]
table = wide[cols].sort_values("selectivity_score", ascending=False)
table.to_csv(os.path.join(RES, "selectivity_table.csv"))
print("=== Selectivity ranking (most HIV-selective first) ===")
print(table[["name", "status", PRIMARY] + OFFTARGETS +
            ["selectivity_score", "selective_for_hiv"]].round(2).to_string())

# ---- Step 21a: KMeans clustering on 3-receptor affinity profile ----
profile = wide[[PRIMARY] + OFFTARGETS].values
Xs = StandardScaler().fit_transform(profile)
k = 3
km = KMeans(n_clusters=k, random_state=0, n_init=10).fit(Xs)
wide["cluster"] = km.labels_
table2 = wide[["name", "status", PRIMARY] + OFFTARGETS + ["selectivity_score", "cluster"]]
table2.to_csv(os.path.join(RES, "clusters.csv"))

# ---- Step 21b: correlation of selectivity with descriptors ----
corr = wide[["selectivity_score", "MW", "logP", "TPSA", "HBD", "HBA", "RotB"]].corr()["selectivity_score"]
corr.to_csv(os.path.join(RES, "selectivity_descriptor_correlation.csv"))
print("\n=== Pearson r: selectivity vs descriptors ===")
print(corr.round(3).to_string())

# ================= FIGURES =================
plt.rcParams.update({"figure.dpi": 130, "font.size": 9})

# Fig 1: affinity heatmap
fig, ax = plt.subplots(figsize=(6, 8))
hm = wide.sort_values("selectivity_score", ascending=False)
mat = hm[[PRIMARY] + OFFTARGETS].values
im = ax.imshow(mat, cmap="viridis_r", aspect="auto")
ax.set_xticks(range(3)); ax.set_xticklabels(["HIV-1\nprotease", "Cathepsin D", "BACE1"])
ax.set_yticks(range(len(hm))); ax.set_yticklabels(hm["name"], fontsize=7)
for i in range(len(hm)):
    for j in range(3):
        ax.text(j, i, f"{mat[i,j]:.1f}", ha="center", va="center", color="w", fontsize=6)
plt.colorbar(im, label="Vina affinity (kcal/mol)")
ax.set_title("Docking affinity: library x receptor")
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig1_affinity_heatmap.png")); plt.close()

# Fig 2: selectivity bar plot
fig, ax = plt.subplots(figsize=(7, 6))
sd = wide.sort_values("selectivity_score")
colors = {"approved": "#2166ac", "discontinued": "#b2182b",
          "investigational": "#762a83", "booster": "#999999"}
ax.barh(sd["name"], sd["selectivity_score"],
        color=[colors.get(s, "#333") for s in sd["status"]])
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("Selectivity  min(ΔΔG)  (kcal/mol)   →  more HIV-selective")
ax.set_title("HIV-1 protease selectivity vs human off-targets")
handles = [plt.Rectangle((0,0),1,1,color=c) for c in colors.values()]
ax.legend(handles, colors.keys(), fontsize=8, loc="lower right")
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig2_selectivity.png")); plt.close()

# Fig 3: PCA of binding profiles, colored by cluster
pca = PCA(n_components=2).fit(Xs); pcs = pca.transform(Xs)
fig, ax = plt.subplots(figsize=(6.5, 5.5))
sc = ax.scatter(pcs[:,0], pcs[:,1], c=wide["cluster"], cmap="Set1", s=60, edgecolor="k")
for i, nm in enumerate(wide["name"]):
    ax.annotate(nm, (pcs[i,0], pcs[i,1]), fontsize=6, xytext=(3,3), textcoords="offset points")
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.0f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.0f}%)")
ax.set_title("Clustering of compounds by 3-receptor binding profile")
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig3_clusters_pca.png")); plt.close()

# Fig 4: HIV vs off-target scatter (selectivity landscape)
fig, ax = plt.subplots(figsize=(6.5, 5.5))
ax.scatter(wide[PRIMARY], wide[OFFTARGETS].min(axis=1),
           c=[colors.get(s, "#333") for s in wide["status"]], s=55, edgecolor="k")
lims = [min(wide[PRIMARY].min(), wide[OFFTARGETS].min().min())-0.5,
        max(wide[PRIMARY].max(), wide[OFFTARGETS].max().max())+0.5]
ax.plot(lims, lims, "k--", lw=0.8, label="equal affinity")
for i, nm in enumerate(wide["name"]):
    ax.annotate(nm, (wide[PRIMARY].iloc[i], wide[OFFTARGETS].min(axis=1).iloc[i]),
                fontsize=6, xytext=(3,3), textcoords="offset points")
ax.set_xlabel("HIV-1 protease affinity (kcal/mol)")
ax.set_ylabel("Best off-target affinity (kcal/mol)")
ax.set_title("Points below the diagonal = HIV-selective")
plt.tight_layout(); plt.savefig(os.path.join(RES, "fig4_selectivity_landscape.png")); plt.close()

n_sel = int(wide["selective_for_hiv"].sum())
print(f"\n{n_sel}/{len(wide)} compounds are HIV-selective (bind HIV-1 protease "
      f"tighter than both off-targets).")
print(f"Figures + tables written to {RES}/")
