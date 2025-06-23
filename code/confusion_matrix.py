import os
import glob
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score

################################ PSSM ################################

pssm_results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/PSSM"
pssm_beta_values = ["50"]  # Edit as needed
pssm_output_file = "results/metrics/pssm_metrics.txt"
os.makedirs(os.path.dirname(pssm_output_file), exist_ok=True)

with open(pssm_output_file, "w") as out_f:
    for allele_folder in os.listdir(pssm_results_dir):
        allele_path = os.path.join(pssm_results_dir, allele_folder)
        allele = allele_folder.split(".")[0]
        if not os.path.isdir(allele_path):
            continue
        for model_folder in os.listdir(allele_path):
            if not model_folder.startswith("b."):
                continue
            beta_val = model_folder[2:].strip()
            if beta_val not in pssm_beta_values:
                continue
            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue
            #for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
            for eval_file in glob.glob(os.path.join(model_path, "concat.eval")):
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]
                cm = confusion_matrix(bin_measured, bin_pred)
                f1 = f1_score(bin_measured, bin_pred, average="binary")
                acc = accuracy_score(bin_measured, bin_pred)
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                out_f.write(f"Confusion matrix for {allele} beta={beta_val} {eval_base}:\n")
                out_f.write(str(cm) + "\n")
                out_f.write(f"F1 score: {f1:.4f}\n")
                out_f.write(f"Accuracy: {acc:.4f}\n\n")

################################ SMM - Gradient Descent ################################

gd_results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/SMM/gradient_decent"
gd_lambda_values = ["0.1"]  # Edit as needed
gd_output_file = "results/metrics/gd_metrics.txt"
os.makedirs(os.path.dirname(gd_output_file), exist_ok=True)

with open(gd_output_file, "w") as out_f:
    for allele_folder in os.listdir(gd_results_dir):
        allele_path = os.path.join(gd_results_dir, allele_folder)
        allele = allele_folder.split(".")[0]
        if not os.path.isdir(allele_path):
            continue
        for model_folder in os.listdir(allele_path):
            if not any(model_folder == f"l.{l}" for l in gd_lambda_values):
                continue
            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue
            #for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
            for eval_file in glob.glob(os.path.join(model_path, "concat.eval")):
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]
                cm = confusion_matrix(bin_measured, bin_pred)
                f1 = f1_score(bin_measured, bin_pred, average="binary")
                acc = accuracy_score(bin_measured, bin_pred)
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                l_part = model_folder[2:] if model_folder.startswith("l.") else model_folder
                out_f.write(f"Confusion matrix for {allele} l={l_part} {eval_base}:\n")
                out_f.write(str(cm) + "\n")
                out_f.write(f"F1 score: {f1:.4f}\n")
                out_f.write(f"Accuracy: {acc:.4f}\n\n")

################################ SMM - Monte Carlo ################################

mc_results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/SMM/monte_carlo"
mc_lambda_values = ["0.1"]      # Edit as needed
mc_t_values = ["100"]           # Edit as needed
mc_output_file = "results/metrics/mc_metrics.txt"
os.makedirs(os.path.dirname(mc_output_file), exist_ok=True)

with open(mc_output_file, "w") as out_f:
    for allele_folder in os.listdir(mc_results_dir):
        allele_path = os.path.join(mc_results_dir, allele_folder)
        allele = allele_folder.split(".")[0]
        if not os.path.isdir(allele_path):
            continue
        for model_folder in os.listdir(allele_path):
            parts = model_folder.split("_")
            l_part = None
            t_part = None
            for part in parts:
                if part.startswith("l."):
                    l_part = part[2:]
                if part.startswith("t."):
                    t_part = part[2:]
            if l_part not in mc_lambda_values or t_part not in mc_t_values:
                continue
            model_path = os.path.join(allele_path, model_folder)
            if not os.path.isdir(model_path):
                continue
            #for eval_file in glob.glob(os.path.join(model_path, "e00*.eval")):
            for eval_file in glob.glob(os.path.join(model_path, "concat.eval")):
                df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                bin_pred = df.iloc[:, -2]
                bin_measured = df.iloc[:, -1]
                cm = confusion_matrix(bin_measured, bin_pred)
                acc = accuracy_score(bin_measured, bin_pred)
                f1 = f1_score(bin_measured, bin_pred, average="binary")
                eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                out_f.write(f"Confusion matrix for {allele} l={l_part} t={t_part} {eval_base}:\n")
                out_f.write(str(cm) + "\n")
                out_f.write(f"F1 score: {f1:.4f}\n")
                out_f.write(f"Accuracy: {acc:.4f}\n\n")

################################ ANN ################################

ann_results_dir = "/mnt/c/Users/nicol/OneDrive - Danmarks Tekniske Universitet/Algorithms in bioinformatics F25/BioAlgoProject2025/results/ANN"
encodings = ["blosum"]         # e.g. ["blosum", "sparse"]
hidden_units = ["5"]      # e.g. ["3", "5", "10"]
epi_values = ["0.05"]   # e.g. ["0.1", "0.5"]
ann_output_file = "results/metrics/ann_metrics.txt"
os.makedirs(os.path.dirname(ann_output_file), exist_ok=True)

with open(ann_output_file, "w") as out_f:
    for allele_folder in os.listdir(ann_results_dir):
        allele_path = os.path.join(ann_results_dir, allele_folder)
        allele = allele_folder.split(".")[0]
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
                    #for eval_file in glob.glob(os.path.join(epi_path, "e00*.eval")):
                    for eval_file in glob.glob(os.path.join(model_path, "concat.eval")):
                        df = pd.read_csv(eval_file, delim_whitespace=True, header=None)
                        bin_pred = df.iloc[:, -2]
                        bin_measured = df.iloc[:, -1]
                        cm = confusion_matrix(bin_measured, bin_pred)
                        f1 = f1_score(bin_measured, bin_pred, average="binary")
                        acc = accuracy_score(bin_measured, bin_pred)
                        eval_base = os.path.splitext(os.path.basename(eval_file))[0]
                        out_f.write(f"Confusion matrix for {allele} encoding={encoding} hidden={hidden_val} epi={epi_val} {eval_base}:\n")
                        out_f.write(str(cm) + "\n")
                        out_f.write(f"F1 score: {f1:.4f}\n")
                        out_f.write(f"Accuracy: {acc:.4f}\n\n")