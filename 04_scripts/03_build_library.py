"""
Step 6: curate the screening library of HIV-1 protease inhibitors.

Fetches isomeric (stereochemistry-aware) SMILES + formula + MW from PubChem for a
curated set of protease inhibitors: FDA-approved agents plus discontinued/
investigational compounds ("clinical failures"). Writes 02_ligands/library.csv.

Run:  .venv/bin/python 04_scripts/03_build_library.py
"""
import csv, json, os, time, urllib.parse, urllib.request

BASE = "/Users/max/Desktop/HIV_Selectivity_Project"
OUT = os.path.join(BASE, "02_ligands", "library.csv")

# (query_name, display_name, status)  status in {approved, investigational, discontinued, booster}
COMPOUNDS = [
    ("saquinavir",   "Saquinavir",   "approved"),
    ("ritonavir",    "Ritonavir",    "approved"),
    ("indinavir",    "Indinavir",    "approved"),
    ("nelfinavir",   "Nelfinavir",   "approved"),
    ("amprenavir",   "Amprenavir",   "approved"),
    ("lopinavir",    "Lopinavir",    "approved"),
    ("atazanavir",   "Atazanavir",   "approved"),
    ("fosamprenavir","Fosamprenavir","approved"),
    ("tipranavir",   "Tipranavir",   "approved"),
    ("darunavir",    "Darunavir",    "approved"),
    ("cobicistat",   "Cobicistat",   "booster"),        # ritonavir-analog CYP3A booster, not a PI itself
    ("brecanavir",   "Brecanavir",   "discontinued"),
    ("mozenavir",    "Mozenavir",    "discontinued"),    # DMP-450, cyclic urea
    ("palinavir",    "Palinavir",    "discontinued"),
    ("telinavir",    "Telinavir",    "discontinued"),    # SC-52151
    ("KNI-272",      "KNI-272",      "investigational"),
    ("GS-8374",      "GS-8374",      "investigational"),  # phosphonate PI
    ("DMP-323",      "DMP-323",      "discontinued"),     # cyclic urea
    ("TMC-310911",   "TMC-310911",   "investigational"),  # aka ASC-09
    ("GRL-0519",     "GRL-0519",     "investigational"),
    ("Ro 31-8959",   "Ro-31-8959",   "approved"),         # saquinavir free base synonym (kept if distinct CID)
    ("darunavir ethanolate", "Darunavir-ethanolate", "approved"),
]

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def fetch(name):
    """Return dict(cid, smiles, formula, mw) or None."""
    # try isomeric SMILES; PubChem renamed some property tokens, so fall back to SMILES
    for prop in ("IsomericSMILES", "SMILES", "CanonicalSMILES"):
        url = f"{PUG}/compound/name/{urllib.parse.quote(name)}/property/{prop},MolecularFormula,MolecularWeight/JSON"
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.load(r)
            props = data["PropertyTable"]["Properties"][0]
            smi = props.get(prop) or props.get("SMILES") or props.get("IsomericSMILES")
            if smi:
                return {
                    "cid": props.get("CID"),
                    "smiles": smi,
                    "formula": props.get("MolecularFormula", ""),
                    "mw": props.get("MolecularWeight", ""),
                }
        except Exception:
            continue
    return None

rows, seen_cids = [], set()
for query, display, status in COMPOUNDS:
    info = fetch(query)
    time.sleep(0.25)  # be polite to PubChem
    if not info:
        print(f"  MISS  {display:18s} ({query})")
        continue
    if info["cid"] in seen_cids:
        print(f"  DUP   {display:18s} CID {info['cid']} already have — skipping")
        continue
    seen_cids.add(info["cid"])
    ligid = display.lower().replace(" ", "_").replace("-", "")
    rows.append({
        "ligand_id": ligid, "name": display, "status": status,
        "cid": info["cid"], "formula": info["formula"], "mw": info["mw"],
        "smiles": info["smiles"],
    })
    print(f"  OK    {display:18s} CID {info['cid']:>9}  {info['formula']}")

# structure-based dedup: collapse solvates/free-base-vs-salt to one entry per molecule
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
_chooser = rdMolStandardize.LargestFragmentChooser()

def desalt_and_key(smiles):
    """Strip counterions/solvents (keep largest fragment); return (clean_smiles, connectivity_key)."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    m = _chooser.choose(m)
    clean = Chem.MolToSmiles(m)
    key = Chem.InchiToInchiKey(Chem.MolToInchi(m)).split("-")[0]
    return clean, key

unique, seen_keys = [], {}
for r in rows:
    clean, key = desalt_and_key(r["smiles"])
    if clean:
        r["smiles"] = clean            # store the desalted structure used downstream
    if key and key in seen_keys:
        print(f"  DEDUP {r['name']:18s} same structure as {seen_keys[key]} — dropped")
        continue
    if key:
        seen_keys[key] = r["name"]
    unique.append(r)
rows = unique

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["ligand_id","name","status","cid","formula","mw","smiles"])
    w.writeheader()
    w.writerows(rows)

print(f"\nLibrary size: {len(rows)} unique compounds -> {OUT}")
from collections import Counter
print("By status:", dict(Counter(r["status"] for r in rows)))
