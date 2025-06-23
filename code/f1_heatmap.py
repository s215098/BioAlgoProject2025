import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

metrics_dir = "results/metrics"
algorithms = ["ann", "gd", "mc", "pssm"]  # Edit as needed to match your files
f1_data = {}
acc_data = {}
algorithm_display = {
    "ann": "ANN",
    "gd": "SMM (GD)",
    "mc": "SMM (MC)",
    "pssm": "PSSM"
}

# Parse each algorithm's metrics file
for algo in algorithms:
    file_path = os.path.join(metrics_dir, f"{algo}_metrics.txt")
    if not os.path.isfile(file_path):
        continue
    with open(file_path) as f:
        lines = f.readlines()
    # Extract allele, F1 score, and accuracy from each block
    for i, line in enumerate(lines):
        if line.startswith("Confusion matrix for"):
            match = re.search(r"for\s+([A-Za-z0-9_.-]+).*e(\d+):", line)
            if match:
                allele = match.group(1)
                eval_id = match.group(2)
                # Next lines: look for F1 and Accuracy
                f1 = None
                acc = None
                for j in range(i+1, min(i+8, len(lines))):
                    if lines[j].startswith("F1 score:"):
                        f1 = float(lines[j].split(":")[1].strip())
                    if lines[j].startswith("Accuracy:"):
                        acc = float(lines[j].split(":")[1].strip())
                    if f1 is not None and acc is not None:
                        break
                if f1 is not None:
                    f1_data.setdefault(algo, {}).setdefault(allele, []).append(f1)
                if acc is not None:
                    acc_data.setdefault(algo, {}).setdefault(allele, []).append(acc)

allele_order = ["A0201", "A1101", "A0202", "B0702", "B1501", "A3001", "B5401", "B5701"]

# Compute average F1 and accuracy for each allele/algorithm
alleles = sorted({allele for algo_dict in f1_data.values() for allele in algo_dict} |
                 {allele for algo_dict in acc_data.values() for allele in algo_dict})
heatmap_data_f1 = pd.DataFrame(index=allele_order, columns=algorithms)
heatmap_data_acc = pd.DataFrame(index=allele_order, columns=algorithms)
for algo in algorithms:
    for allele in alleles:
        scores_f1 = f1_data.get(algo, {}).get(allele, [])
        scores_acc = acc_data.get(algo, {}).get(allele, [])
        heatmap_data_f1.loc[allele, algo] = sum(scores_f1) / len(scores_f1) if scores_f1 else None
        heatmap_data_acc.loc[allele, algo] = sum(scores_acc) / len(scores_acc) if scores_acc else None
        
# Rename columns for display
heatmap_data_f1_disp = heatmap_data_f1.rename(columns=algorithm_display)
heatmap_data_acc_disp = heatmap_data_acc.rename(columns=algorithm_display)

# Plot F1 heatmap
plt.figure(figsize=(10, max(4, len(alleles)*0.5)))
sns.heatmap(heatmap_data_f1_disp.astype(float), annot=True, cmap="Blues", fmt=".3f")
plt.title("Averaged F1 Scores per Allele and Algorithm")
plt.ylabel("Allele")
plt.xlabel("Algorithm")
plt.tight_layout()
plt.savefig("results/metrics/f1_heatmap.png")  # Save F1 heatmap
#plt.show()

# Plot Accuracy heatmap
plt.figure(figsize=(10, max(4, len(alleles)*0.5)))
sns.heatmap(heatmap_data_acc_disp.astype(float), annot=True, cmap="Blues", fmt=".3f")
plt.title("Averaged Accuracy per Allele and Algorithm")
plt.ylabel("Allele")
plt.xlabel("Algorithm")
plt.tight_layout()
plt.savefig("results/metrics/accuracy_heatmap.png")  # Save accuracy heatmap
#plt.show()