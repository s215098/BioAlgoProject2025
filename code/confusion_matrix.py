import os
import glob
import pandas as pd
from sklearn.metrics import confusion_matrix

################################ PSSM ################################

pssm_results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/PSSM"
pssm_beta_values = ["10", "50"]  # Edit as needed
pssm_output_file = "results/confusion_matrix/pssm_confusion_matrices.txt"
os.makedirs(os.path.dirname(pssm_output_file), exist_ok=True)

with open(pssm_output_file, "w") as out_f:
    for allele_folder in os.listdir(pssm_results_dir):
        allele_path = os.path.join(pssm_results_dir, allele_folder)
        if not os.path.isdir(allele_path):
            continue
        for model_folder in os.listdir(allele_path):
            if not model_folder.startswith("b."):
                continue
            beta_val = model_folder[2:].strip()
            if beta_val not in pssm_beta_values:
                print(f"Skipping unknown beta value: {beta_val} in {model_folder}")
                continue
            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue
            for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]
                cm = confusion_matrix(bin_measured, bin_pred)
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                out_f.write(f"Confusion matrix for {allele_folder} beta={beta_val} {eval_base}:\n")
                out_f.write(str(cm) + "\n\n")



         
################################ SMM - Gradient Decent ################################
# Set the path to your results directory
results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/SMM/gradient_decent"

# List of lambda values to include (as strings matching your folder names)
# "0", "0.02", "0.1", "0.001"
lambda_values = ["0", "0.02", "0.1"]  # Edit this list as needed

output_file = "results/confusion_matrix/gd_confusion_matrices.txt"
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out_f:
    # Loop through all allele result folders
    for allele_folder in os.listdir(results_dir):
        allele_path = os.path.join(results_dir, allele_folder)
        if not os.path.isdir(allele_path):
            continue

        # Loop through all lambda/model folders
        for model_folder in os.listdir(allele_path):
            # Only process folders with the desired lambda values
            if not any(model_folder == f"l.{l}" for l in lambda_values):
                continue

            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue

            # Find all e00*.eval files
            for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
                # Read the eval file (assuming whitespace-separated)
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                # bin_pred = 3rd last column, bin_measured = 2nd last column
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]

                # Compute confusion matrix
                cm = confusion_matrix(bin_measured, bin_pred)
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                l_part = model_folder[2:] if model_folder.startswith("l.") else model_folder
                out_f.write(f"Confusion matrix for {allele_folder} l={l_part} {eval_base}:\n")
                out_f.write(str(cm) + "\n\n")
                



################################ SMM - Monte Carlo ################################
# Set the path to your results directory
results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/SMM/monte_carlo"

# List of lambda and t values to include (as strings)
# l: "0", "0.02", "0.1", "0.001"
# t: "10", "20", "50", "100"
lambda_values = ["0.001", "0.1"]      # Edit as needed
t_values = ["20", "100"]           # Edit as needed

output_file = "results/confusion_matrix/mc_confusion_matrices.txt"
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out_f:
    for allele_folder in os.listdir(results_dir):
        allele_path = os.path.join(results_dir, allele_folder)
        if not os.path.isdir(allele_path):
            continue

        for model_folder in os.listdir(allele_path):
            # Folder name pattern: l.<lambda>_t.<t>
            parts = model_folder.split("_")
            l_part = None
            t_part = None
            for part in parts:
                if part.startswith("l."):
                    l_part = part[2:]
                if part.startswith("t."):
                    t_part = part[2:]
            if l_part not in lambda_values or t_part not in t_values:
                continue

            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue

            for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]
                cm = confusion_matrix(bin_measured, bin_pred)
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                out_f.write(f"Confusion matrix for {allele_folder} l={l_part} t={t_part} {eval_base}:\n")
                out_f.write(str(cm) + "\n\n")




################################ ANN ################################
results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/ANN"

# Choose which encodings, hidden units, and epi values to include
encodings = ["blosum"]         # e.g. ["blosum", "sparse"]
hidden_units = ["5"]      # e.g. ["3", "5", "10"]
epi_values = ["0.05", "0.1"]   # e.g. ["0.1", "0.5"]

output_file = "results/confusion_matrix/ann_confusion_matrices.txt"
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out_f:
    for allele_folder in os.listdir(results_dir):
        allele_path = os.path.join(results_dir, allele_folder)
        if not os.path.isdir(allele_path):
            continue

        for encoding in os.listdir(allele_path):
            if encoding not in encodings:
                continue
            encoding_path = os.path.join(allele_path, encoding)
            if not os.path.isdir(encoding_path):
                continue

            for hidden_folder in os.listdir(encoding_path):
                if not hidden_folder.startswith("hidden."):
                    continue
                hidden_val = ".".join(hidden_folder.split(".")[1:])
                if not any(float(hidden_val) == float(h) for h in hidden_units):
                    continue
                hidden_path = os.path.join(encoding_path, hidden_folder)
                if not os.path.isdir(hidden_path):
                    continue


                for epi_folder in os.listdir(hidden_path):
                    if not epi_folder.startswith("epi."):
                        continue
                    epi_val = ".".join(epi_folder.split(".")[1:])
                    if epi_val not in epi_values:
                        continue
                    epi_path = os.path.join(hidden_path, epi_folder)
                    if not os.path.isdir(epi_path):
                        continue

                    for eval_file in glob.glob(os.path.join(epi_path, "e00*.eval")):
                        df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                        bin_pred = df.iloc[:, -2]
                        bin_measured = df.iloc[:, -1]
                        cm = confusion_matrix(bin_measured, bin_pred)
                        eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                        out_f.write(f"Confusion matrix for {allele_folder} encoding={encoding} hidden={hidden_val} epi={epi_val} {eval_base}:\n")
                        out_f.write(str(cm) + "\n\n")