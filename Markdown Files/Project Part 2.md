# Phase 6: De Novo Drug Design & Ligand Optimization

## Project Tracking Updates
Please add these steps to your Notion Tracker or project management database to monitor progress for this phase:

* Step 27: Install Cheminformatics libraries (rdkit, pandas) into the Python environment.
* Step 28: Obtain the Canonical SMILES string for the most selective drug from Phase 4 (e.g., Tipranavir or Ritonavir).
* Step 29: Write an AI-assisted RDKit script to generate 10-20 structural analogs using standard medicinal chemistry substitutions.
* Step 30: Apply a Lipinski's Rule of 5 filter script to eliminate biologically inviable compounds.
* Step 31: Convert the passing 2D molecules into 3D .pdbqt files and run them through the AutoDock Vina pipeline.

---

## Project Context & Objective
In previous phases, we utilized AutoDock Vina to screen known HIV protease inhibitors against both the primary viral target and human off-targets (Cathepsin D, BACE1). By calculating the Selectivity Score, we identified which structural features correlate with high selectivity. 

In this phase, we act as molecular architects. We will take our most selective drug, translate its 3D chemical structure into a machine-readable format, and use Python (RDKit) and Generative AI (GitHub Copilot) to invent completely new, de novo chemical analogs. We will then mathematically filter these new drugs to ensure they can survive the human body's physiological environment.

---

## Step 27: Environment Setup & Cheminformatics Tools

To manipulate chemical structures programmatically, pharmaceutical companies rely on RDKit, an open-source cheminformatics toolkit for Python. 

### The Chemistry Concept
In Chem 51 (Organic Chemistry), you learned to draw structures using skeletal formulas (hexagons, zig-zags). Computers cannot parse images easily. Instead, we use SMILES (Simplified Molecular-Input Line-Entry System), which translates 3D molecules into a 1D string of text. (e.g., Benzene is c1ccccc1). RDKit allows Python to read these text strings and perform complex organic chemistry operations on them.

### Lab Procedure
1. Open your computational-biology project folder in VS Code.
2. Open your terminal (Terminal -> New Terminal).
3. Ensure your Python environment is active, then run the following command to install the required libraries:
   
   pip install rdkit pandas
   
   (Note: RDKit handles the chemical logic, while Pandas will organize our generated drugs into a clean data frame).

---

## Step 28: Defining the Seed Molecule (SMILES)

We need a starting point. Based on our previous docking data, we select a highly selective inhibitor as our "seed." 

1. Go to the PubChem database (pubchem.ncbi.nlm.nih.gov).
2. Search for your best-performing drug (e.g., Tipranavir or Ritonavir).
3. Scroll down to the "Computed Properties" section and copy the Canonical SMILES string.
   * Example (Tipranavir): CCC1=C(C(=O)O[C@@H]1CC2=CC=CC=C2)CC(CC)N(CC3=CC=CC=C3)S(=O)(=O)C4=CC(=CC=C4)C(F)(F)F
4. Save this string; it is the blueprint we will feed into our Python script.

---

## Step 29: AI-Assisted Analog Generation
### Lab Procedure
1. In VS Code, create a new file named 04_scripts/generate_analogs.py.
2. Use GitHub Copilot as your AI TA by opening the Copilot Chat or typing a comment block directly into the file. 

Prompt Copilot with the following instructions:
> "Write a Python script using RDKit and Pandas. Start with the following SMILES string: [PASTE YOUR SMILES HERE]. I want to generate 10 to 20 structural analogs by performing standard medicinal chemistry substitutions. For example, replace a benzene ring with a pyridine ring, swap a methyl group for a fluorine atom, or add a hydroxyl group. Output the results as a Pandas DataFrame with two columns: 'Analog_Name' and 'SMILES', and save it to a CSV."

### The Computational Concept
Behind the scenes, RDKit will parse the SMILES string into a molecular graph (nodes = atoms, edges = bonds). Copilot will generate code that finds specific sub-graphs (like a benzene ring) and digitally "cleaves" the bond to attach a new functional group, simulating an organic synthesis reaction in milliseconds.

---

## Step 30: The Physiology Filter (Lipinski's Rule of 5)

A drug that binds perfectly to a target is useless if it cannot be absorbed by the human body. From Bio Sci 99, we know that oral drugs must cross the hydrophobic phospholipid bilayer of the intestines while remaining soluble in hydrophilic blood. 

We will apply Lipinski's Rule of 5 to filter out biologically inviable compounds:
1. Molecular Weight < 500 Da (Too large = cannot cross membranes).
2. LogP < 5 (Too hydrophobic = trapped in fat tissue).
3. Hydrogen Bond Donors (OH, NH) < 5
4. Hydrogen Bond Acceptors (N, O) < 10

### Lab Procedure
1. Create a new script named 04_scripts/admet_filter.py.
2. Ask GitHub Copilot to write the filtering logic.

Prompt Copilot with the following instructions:
> "Write a Python script that reads the CSV of SMILES generated in the previous step. For each SMILES string, use RDKit's Descriptors module to calculate the Molecular Weight, MolLogP, NumHDonors, and NumHAcceptors. Filter the DataFrame to only keep molecules that strictly obey Lipinski's Rule of 5 (MW < 500, LogP < 5, HBD <= 5, HBA <= 10). Save the passing molecules to a new CSV called 'viable_analogs.csv'."

---

## Step 31: Moving Forward (To Docking)
Once your viable_analogs.csv is generated:
1. We will use a script (via OpenBabel or RDKit) to convert these 1D SMILES strings back into 3D .pdbqt files.
2. We will run these new files through your existing AutoDock Vina pipeline against both HIV Protease and Cathepsin D to see if you have successfully invented a drug with a better Selectivity Score than the FDA-approved ones!
