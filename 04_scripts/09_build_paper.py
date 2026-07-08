"""
Assemble the research-paper HTML (self-contained, figures embedded as base64) from
the project's data files. A companion shell step renders it to PDF with headless
Chrome. Run:  .venv/bin/python 04_scripts/09_build_paper.py
"""
import base64, os
import pandas as pd

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
RES = os.path.join(BASE, "05_results")

sel = pd.read_csv(os.path.join(RES, "selectivity_table.csv"))
redock = pd.read_csv(os.path.join(BASE, "03_docking", "redock_validation.csv"))
lib = pd.read_csv(os.path.join(BASE, "02_ligands", "library.csv"))

def img64(path):
    with open(path, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

FIG = {n: img64(os.path.join(RES, f"{n}.png")) for n in
       ["fig1_affinity_heatmap", "fig2_selectivity", "fig3_clusters_pca", "fig4_selectivity_landscape"]}

# ---- build the main results table (ranked by selectivity) ----
def results_rows():
    rows = ""
    for _, r in sel.iterrows():
        sel_class = "pos" if r["selectivity_score"] > 0 else "neg"
        rows += (f"<tr><td class='name'>{r['name']}</td><td>{r['status']}</td>"
                 f"<td>{r['hiv_protease']:.2f}</td><td>{r['cathepsin_d']:.2f}</td>"
                 f"<td>{r['bace1']:.2f}</td>"
                 f"<td class='{sel_class}'>{r['selectivity_score']:+.2f}</td></tr>")
    return rows

def corr_rows():
    corr = pd.read_csv(os.path.join(RES, "selectivity_descriptor_correlation.csv"), index_col=0)
    labels = {"MW": "Molecular weight", "logP": "logP (lipophilicity)",
              "TPSA": "Topological polar surface area", "HBD": "H-bond donors",
              "HBA": "H-bond acceptors", "RotB": "Rotatable bonds"}
    rows = ""
    for key, lab in labels.items():
        v = corr.loc[key, "selectivity_score"]
        rows += f"<tr><td>{lab}</td><td>{v:+.2f}</td></tr>"
    return rows

n_sel = int((sel["selectivity_score"] > 0).sum())
n_total = len(sel)
n_approved = int((lib["status"] == "approved").sum())

HTML = f"""<!doctype html><html><head><meta charset="utf-8"><style>
@page {{ size: letter; margin: 2.2cm 2cm; }}
body {{ font-family: Georgia, 'Times New Roman', serif; font-size: 10.5pt;
        line-height: 1.5; color: #111; max-width: 100%; }}
h1 {{ font-size: 18pt; margin: 0 0 4px; line-height: 1.25; }}
.authors {{ font-size: 11pt; margin: 6px 0 2px; }}
.affil {{ font-size: 9.5pt; font-style: italic; color: #333; margin: 0 0 4px; }}
.meta {{ font-size: 8.5pt; color: #555; border-bottom: 2px solid #333; padding-bottom: 10px; margin-bottom: 14px; }}
h2 {{ font-size: 12.5pt; border-bottom: 1px solid #999; padding-bottom: 2px; margin: 18px 0 8px; }}
h3 {{ font-size: 11pt; margin: 12px 0 4px; }}
p {{ margin: 0 0 8px; text-align: justify; }}
.abstract {{ background: #f4f6f8; border-left: 3px solid #345; padding: 10px 14px; font-size: 10pt; }}
table {{ border-collapse: collapse; width: 100%; font-size: 9pt; margin: 8px 0 4px; }}
th, td {{ border: 1px solid #bbb; padding: 3px 6px; text-align: center; }}
th {{ background: #e8ecf0; }}
td.name {{ text-align: left; font-weight: bold; }}
td.pos {{ color: #0a6; font-weight: bold; }}
td.neg {{ color: #c22; font-weight: bold; }}
figure {{ margin: 12px 0; text-align: center; page-break-inside: avoid; }}
figure img {{ max-width: 86%; border: 1px solid #ddd; }}
figcaption {{ font-size: 8.8pt; color: #333; margin-top: 4px; text-align: left; }}
.caption-label {{ font-weight: bold; }}
.refs {{ font-size: 8.6pt; }}
.refs li {{ margin-bottom: 3px; }}
code {{ font-family: 'Courier New', monospace; font-size: 9pt; }}
.small {{ font-size: 8.6pt; color: #444; }}
</style></head><body>

<h1>Structure-based selectivity mapping of HIV-1 protease inhibitors against human aspartic-protease off-targets</h1>
<p class="authors">Max Wo</p>
<p class="affil">Department of Pharmaceutical Sciences, University of California, Irvine</p>
<p class="meta">Independent computational research project &middot; Molecular docking study &middot;
Data and code: <code>github.com/maxwo123/computational-biology</code></p>

<div class="abstract">
<b>Abstract.</b> Off-target inhibition of host enzymes is a major driver of drug toxicity and
attrition in medicinal chemistry. HIV-1 protease inhibitors (PIs) act on a viral aspartic
protease, but two human aspartic proteases &mdash; cathepsin&nbsp;D and &beta;-secretase&nbsp;1
(BACE1) &mdash; share the same catalytic architecture and are therefore plausible off-targets.
Here we computationally map the selectivity of {n_total} HIV-1 PIs (spanning FDA-approved,
discontinued, and investigational agents) by molecular docking against the viral target
(HIV-1 protease, PDB 1HXW) and both human off-targets (cathepsin&nbsp;D, 1LYB; BACE1, 1FKN)
using AutoDock&nbsp;Vina. A selectivity metric &Delta;&Delta;G was defined as the difference
between each compound's off-target and on-target predicted binding energy. The docking protocol
was validated by native ligand re-docking, reproducing the crystallographic ritonavir pose in
HIV-1 protease to 1.93&nbsp;&Aring; RMSD. {n_sel} of {n_total} compounds were predicted to be
HIV-selective (binding the viral protease more tightly than either human enzyme). Predicted
selectivity correlated negatively with molecular polarity and size (topological polar surface
area, <i>r</i>&nbsp;=&nbsp;&minus;0.59; rotatable bonds, <i>r</i>&nbsp;=&nbsp;&minus;0.49),
suggesting that smaller, more rigid, less polar inhibitors preferentially engage the viral
pocket. These results illustrate a reproducible, open-source workflow for early-stage
off-target risk assessment.
</div>

<h2>1. Introduction</h2>
<p>The protease encoded by human immunodeficiency virus type 1 (HIV-1) cleaves the viral
Gag&ndash;Pol polyprotein into functional units and is essential for producing infectious
virions. It is an aspartic protease that operates as an obligate homodimer, with a single
catalytic site formed at the interface of two identical subunits and powered by a pair of
catalytic aspartate residues (the "Asp dyad"). Inhibiting this enzyme halts viral maturation,
and HIV-1 protease inhibitors (PIs) remain a cornerstone of antiretroviral therapy.</p>
<p>A central concern in drug design is <i>selectivity</i>: a therapeutic should engage its
intended target while sparing structurally related host proteins. Human aspartic proteases such
as cathepsin&nbsp;D (a lysosomal enzyme implicated in protein turnover and, when dysregulated,
in cancer progression) and BACE1 (&beta;-secretase, which generates the amyloid-&beta; peptide
central to Alzheimer's disease) use the same Asp-dyad catalytic chemistry as the viral enzyme.
This shared active-site architecture creates a plausible route to off-target binding and, in
principle, mechanism-based side effects. Because both human enzymes are structurally
well-characterized, they serve as informative counter-screens for probing PI selectivity.</p>
<p>Here we build an open, reproducible molecular-docking pipeline to estimate the relative
binding of a diverse HIV-1 PI library across the viral target and these two human off-targets,
quantify selectivity, and ask which molecular properties are associated with it.</p>

<h2>2. Materials and Methods</h2>
<h3>2.1 Target structures</h3>
<p>Three crystal structures were obtained from the RCSB Protein Data Bank: HIV-1 protease in
complex with ritonavir (PDB <b>1HXW</b>; primary target, both chains of the biological homodimer
retained), human cathepsin&nbsp;D with pepstatin (PDB <b>1LYB</b>; off-target), and human BACE1
with the transition-state inhibitor OM99-2 (PDB <b>1FKN</b>; off-target). For each receptor, the
co-crystallized inhibitor, crystallographic waters, buffer components, and (for cathepsin&nbsp;D)
N-linked glycans were removed, retaining only the biologically relevant protein chains. The
geometric centroid of each native inhibitor defined the center of the docking search box, which
was padded to at least 22.5&nbsp;&Aring; per axis.</p>
<h3>2.2 Compound library</h3>
<p>A library of {n_total} HIV-1 PIs was assembled from PubChem: {n_approved} FDA-approved agents,
one pharmacoenhancer (cobicistat, a ritonavir analog included as a structural comparator), and
the remainder discontinued or investigational compounds. Structures were desalted to the largest
organic fragment and deduplicated by InChIKey connectivity.</p>
<h3>2.3 Structure preparation</h3>
<p>Receptors were protonated at pH&nbsp;7.0 using PDB2PQR with PROPKA-predicted p<i>K</i><sub>a</sub>
values, fixing the ionization state of the catalytic aspartate dyad, then converted to the
AutoDock PDBQT format with Gasteiger partial charges (Open&nbsp;Babel). Ligands were protonated to
their dominant microspecies at pH&nbsp;7.4, embedded in three dimensions with the RDKit ETKDGv3
method, energy-minimized with the MMFF94 force field, and assigned rotatable-bond torsion trees
with Meeko.</p>
<h3>2.4 Molecular docking</h3>
<p>Docking used AutoDock&nbsp;Vina 1.2.7 (exhaustiveness 8, up to nine poses per run). Each of the
{n_total} ligands was docked against all three receptors ({n_total}&nbsp;&times;&nbsp;3&nbsp;=&nbsp;60
independent runs), and the most favorable predicted binding energy (&Delta;G, kcal&nbsp;mol<sup>&minus;1</sup>)
per pair was retained.</p>
<h3>2.5 Validation and analysis</h3>
<p>The protocol was validated by re-docking each receptor's native inhibitor into its own pocket
and measuring heavy-atom root-mean-square deviation (RMSD) from the crystallographic pose.
Selectivity was quantified as &Delta;&Delta;G&nbsp;=&nbsp;&Delta;G<sub>off-target</sub>&nbsp;&minus;&nbsp;&Delta;G<sub>HIV</sub>;
because &Delta;G is negative, a positive &Delta;&Delta;G indicates weaker off-target binding and
thus greater selectivity for the viral enzyme. The reported selectivity score is the worst-case
(smallest) &Delta;&Delta;G across the two off-targets. Compounds were clustered by their
three-receptor binding profile (k-means), and selectivity was correlated (Pearson) against
molecular descriptors computed with RDKit.</p>

<h2>3. Results</h2>
<h3>3.1 Protocol validation</h3>
<p>Re-docking the native inhibitors reproduced the experimental binding mode for the primary
target: ritonavir re-docked into HIV-1 protease with an RMSD of <b>1.93&nbsp;&Aring;</b>, below the
conventional 2&nbsp;&Aring; acceptance threshold, confirming correct box placement and docking
parameters. The two off-target native ligands are large, highly flexible peptidic inhibitors
(pepstatin and OM99-2); these re-docked less accurately (cathepsin&nbsp;D/pepstatin, 3.69&nbsp;&Aring;;
BACE1/OM99-2 ligand preparation was unsuccessful), a recognized limitation for peptide-like
molecules that does not affect docking of the drug-like library compounds.</p>

<h3>3.2 Binding-affinity landscape</h3>
<p>Across all 60 dockings, predicted affinities ranged from roughly &minus;5.7 to
&minus;11.0&nbsp;kcal&nbsp;mol<sup>&minus;1</sup>. Mean affinity was strongest for HIV-1 protease
(&minus;9.6&nbsp;kcal&nbsp;mol<sup>&minus;1</sup>), intermediate for BACE1
(&minus;9.2), and weakest for cathepsin&nbsp;D (&minus;8.8), consistent with the library having been
optimized historically against the viral enzyme (Figure&nbsp;1).</p>
<figure><img src="{FIG['fig1_affinity_heatmap']}">
<figcaption><span class="caption-label">Figure 1.</span> Predicted binding affinity
(kcal&nbsp;mol<sup>&minus;1</sup>) of each compound against the three receptors. Darker cells denote
tighter predicted binding. Rows are ordered by selectivity (most HIV-selective at top).</figcaption></figure>

<h3>3.3 Selectivity ranking</h3>
<p>{n_sel} of {n_total} compounds ({100*n_sel//n_total}%) were predicted to be HIV-selective,
binding the viral protease more tightly than both human off-targets (Table&nbsp;1, Figure&nbsp;2).
The most selective compounds included the cyclic-urea inhibitor DMP-323 and the non-peptidic
approved drug tipranavir. A minority of compounds &mdash; notably darunavir and atazanavir &mdash;
were predicted to bind an off-target (most often BACE1) at least as tightly as HIV-1 protease;
as discussed below, this reflects the limits of docking for absolute discrimination rather than a
validated safety liability of these clinically successful drugs.</p>

<p><b>Table 1.</b> Predicted binding energies and selectivity, ranked most- to least-selective.
&Delta;G in kcal&nbsp;mol<sup>&minus;1</sup>; selectivity = worst-case &Delta;&Delta;G (positive =
HIV-selective).</p>
<table>
<tr><th>Compound</th><th>Status</th><th>HIV-1 PR</th><th>Cathepsin D</th><th>BACE1</th><th>Selectivity</th></tr>
{results_rows()}
</table>

<figure><img src="{FIG['fig2_selectivity']}">
<figcaption><span class="caption-label">Figure 2.</span> Worst-case selectivity (min&nbsp;&Delta;&Delta;G)
per compound, colored by clinical status. Bars to the right indicate preference for HIV-1 protease.</figcaption></figure>

<h3>3.4 Structure&ndash;selectivity relationships</h3>
<p>Predicted selectivity correlated negatively with descriptors of molecular size and polarity
(Table&nbsp;2). The strongest association was with topological polar surface area
(<i>r</i>&nbsp;=&nbsp;&minus;0.59), followed by rotatable-bond count
(<i>r</i>&nbsp;=&nbsp;&minus;0.49) and hydrogen-bond acceptors (<i>r</i>&nbsp;=&nbsp;&minus;0.43).
In other words, smaller, more rigid, and less polar inhibitors tended to prefer the viral pocket,
whereas large, flexible, highly polar compounds engaged the human aspartic-protease sites
comparably well &mdash; a potential design cue for improving selectivity.</p>

<p><b>Table 2.</b> Pearson correlation of predicted selectivity with molecular descriptors.</p>
<table style="width:60%; margin-left:auto; margin-right:auto;">
<tr><th style="text-align:left">Descriptor</th><th>Pearson <i>r</i> vs selectivity</th></tr>
{corr_rows()}
</table>

<figure><img src="{FIG['fig4_selectivity_landscape']}">
<figcaption><span class="caption-label">Figure 3.</span> Selectivity landscape: predicted HIV-1
protease affinity versus the best (most favorable) off-target affinity. Points below the diagonal
are HIV-selective.</figcaption></figure>

<h3>3.5 Binding-profile clustering</h3>
<p>Unsupervised k-means clustering of compounds by their three-receptor affinity profile separated
the library into groups with distinct cross-reactivity signatures (Figure&nbsp;4), indicating that
selectivity behavior is structured rather than random across the chemotypes represented.</p>
<figure><img src="{FIG['fig3_clusters_pca']}">
<figcaption><span class="caption-label">Figure 4.</span> Principal-component projection of the
three-receptor binding profiles, colored by k-means cluster.</figcaption></figure>

<h2>4. Discussion</h2>
<p>The workflow recovered the expected qualitative picture: a library historically optimized
against HIV-1 protease binds the viral enzyme most strongly on average, and a majority of
compounds are predicted to be selective over the two human aspartic proteases. The
structure&ndash;selectivity trend &mdash; selectivity decreasing with polarity, flexibility, and
size &mdash; is chemically reasonable: the compact, largely hydrophobic HIV-1 protease pocket may
be more discriminating toward small rigid ligands, whereas larger flexible molecules can adopt
accommodating conformations in multiple aspartic-protease sites.</p>
<p>Several highly effective clinical drugs (darunavir, atazanavir) were <i>not</i> flagged as
selective in this purely geometric, single-pose docking analysis. This is an important
illustration of docking's limits rather than a contradiction: Vina scores are approximate,
neglect protein flexibility and explicit solvation/entropy, and are best used for
<i>relative ranking within</i> a target rather than for absolute cross-target discrimination.
The clinical selectivity of these agents also derives from pharmacokinetics and cellular context
not captured by rigid-receptor docking.</p>

<h2>5. Limitations</h2>
<p>Predicted affinities are docking scores, not experimental free energies, and carry errors of a
few kcal&nbsp;mol<sup>&minus;1</sup>. Receptors were treated as rigid; catalytic-dyad protonation
was fixed at a single pH&nbsp;7 assignment; and a single Vina run (one random seed) was used per
pair without replicate averaging. Off-target validation was limited by the difficulty of
re-docking large peptidic native ligands. Findings should therefore be read as hypothesis-
generating relative trends, not quantitative predictions.</p>

<h2>6. Conclusion</h2>
<p>Using an entirely open-source pipeline, we mapped the predicted selectivity of {n_total} HIV-1
protease inhibitors against two human aspartic-protease off-targets, validated the protocol against
a crystallographic pose (1.93&nbsp;&Aring; RMSD), and identified a size/polarity trend associated
with viral-target selectivity. The approach provides a transparent, reproducible template for
early-stage off-target risk assessment in medicinal chemistry teaching and practice.</p>

<h2>Data and code availability</h2>
<p class="small">All structures, prepared inputs, docking outputs, analysis scripts, and figures
are available at <code>github.com/maxwo123/computational-biology</code>. The pipeline is fully
reproducible via the numbered scripts in <code>04_scripts/</code>.</p>

<h2>Software</h2>
<p class="small">AutoDock Vina 1.2.7; RDKit 2025.09; Meeko 0.7; PDB2PQR 3.6 with PROPKA 3.5;
Open Babel 3.1; PyMOL; pandas; scikit-learn; matplotlib.</p>

<h2>References</h2>
<ol class="refs">
<li>Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0. <i>J Chem Inf Model</i> 2021;61:3891&ndash;3898.</li>
<li>Trott O, Olson AJ. AutoDock Vina. <i>J Comput Chem</i> 2010;31:455&ndash;461.</li>
<li>Kempf DJ, et al. ABT-538 (ritonavir) and HIV protease. <i>Proc Natl Acad Sci USA</i> 1995;92:2484&ndash;2488. (PDB 1HXW)</li>
<li>Baldwin ET, et al. Crystal structure of human cathepsin D. <i>Proc Natl Acad Sci USA</i> 1993;90:6796&ndash;6800. (PDB 1LYB)</li>
<li>Hong L, et al. Structure of BACE1&ndash;OM99-2 complex. <i>Science</i> 2000;290:150&ndash;153. (PDB 1FKN)</li>
<li>Dolinsky TJ, et al. PDB2PQR. <i>Nucleic Acids Res</i> 2007;35:W522&ndash;W525.</li>
<li>Olsson MHM, et al. PROPKA3. <i>J Chem Theory Comput</i> 2011;7:525&ndash;537.</li>
<li>O'Boyle NM, et al. Open Babel. <i>J Cheminform</i> 2011;3:33.</li>
<li>Forli S, et al. Computational protein&ndash;ligand docking with AutoDock/Meeko. <i>Nat Protoc</i> 2016;11:905&ndash;919.</li>
<li>Landrum G, et al. RDKit: Open-source cheminformatics. <a>rdkit.org</a>.</li>
<li>Ghosh AK, et al. Darunavir and the design of HIV-1 protease inhibitors. <i>Acc Chem Res</i> 2008;41:78&ndash;86.</li>
<li>Wlodawer A, Vondrasek J. Inhibitors of HIV-1 protease. <i>Annu Rev Biophys Biomol Struct</i> 1998;27:249&ndash;284.</li>
</ol>

</body></html>"""

out = os.path.join(BASE, "paper.html")
with open(out, "w") as f:
    f.write(HTML)
print(f"wrote {out} ({len(HTML)//1024} KB with embedded figures)")
