import os
import pandas as pd

# User settings
DATA_DIR = "data/all_files_new"
PSSM_DIR = "data/PSSM"
N_NONBINDERS = {
    "A0201": 237,
    "A0202": 130,
    "A1101": 139,
    "A3001": 16,
    "B0702": 42,
    "B1501": 36,
    "B5401": 17,
    "B5701": 3,
    # Add more as needed
}

for allele in os.listdir(DATA_DIR):
    allele_path = os.path.join(DATA_DIR, allele)
    if not os.path.isdir(allele_path):
        continue
    if allele not in N_NONBINDERS:
        print(f"Allele {allele} not in N_NONBINDERS, skipping.")
        continue
    n = N_NONBINDERS[allele]
    e000_path = os.path.join(allele_path, "e000.csv")
    if not os.path.isfile(e000_path):
        print(f"File not found: {e000_path}")
        continue

    # Use whitespace as separator
    df = pd.read_csv(e000_path, header=None, delim_whitespace=True)
    non_binders = df[df[1] < 0.462]
    selected = non_binders.head(n)

    # Now read the PSSM e000.csv for this allele
    pssm_allele_path = os.path.join(PSSM_DIR, allele)
    pssm_e000_path = os.path.join(pssm_allele_path, "e000.csv")
    if not os.path.isfile(pssm_e000_path):
        print(f"PSSM file not found: {pssm_e000_path}")
        continue

    # Read PSSM e000.csv (assume same format)
    pssm_df = pd.read_csv(pssm_e000_path, header=None, delim_whitespace=True)

    # Combine (append) the two DataFrames
    combined = pd.concat([pssm_df, selected], ignore_index=True)

    # Write to new file in PSSM/<allele>/
    out_path = os.path.join(pssm_allele_path, "e000_with_nonbinders.csv")
    combined.to_csv(out_path, index=False, header=False, sep=" ")
    print(f"Saved combined file to {out_path}")