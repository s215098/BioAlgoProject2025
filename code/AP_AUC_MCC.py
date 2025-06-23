# %% [markdown]
# ## ROC-AUC Plotting

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef

from argparse import ArgumentParser

parser = ArgumentParser(description="Plot ROC and PR curves from prediction data.")
parser.add_argument('-pf', action="store", dest="pred_file", type=str, help='Path to the prediction file.')
args = parser.parse_args()
pred_file = args.pred_file


# # pred_file = "/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/ANN/B5701.res/blosum/hidden.5/epi.0.1/concat.eval"
# pred_data = np.loadtxt(pred_file, dtype=str)

# y_eval_new = np.array(pred_data[:, 2], dtype=float)
# y_preds_eval_new = np.array(pred_data[:, 1], dtype=float)

# # Set threshold for binary classification
# BINDER_THRESHOLD = 0.426

# # Convert predictions and true values into binary classes

# y_eval_class = np.where(y_eval_new.flatten() >= BINDER_THRESHOLD, 1, 0)
# y_preds_class = np.where(np.array(y_preds_eval_new).flatten() >= BINDER_THRESHOLD, 1, 0)

# # print(y_eval_class)
# # print(y_preds_class)

# # Combine into dataframe for peptide length-specific ROC curves
# # Ensure peptides are correctly defined
# peptide_lengths = np.array([len(p) for p in pred_data[:, 0]])  # Use the first column of pred_data for peptide sequences
# pred_per_len = pd.DataFrame({
#     'peptide_length': peptide_lengths,
#     'target': y_eval_class,
#     'prediction': np.array(y_preds_eval_new)
# })
# print(pred_per_len)

# # Function to plot ROC curve for a specific peptide length
# def plot_roc_curve(fpr, tpr, roc_auc, length):
#     plt.plot(fpr, tpr, label=f'AUC = {roc_auc:.2f} ({length}-mer)')

# # Plot ROC curve for each peptide length
# plt.figure(figsize=(7, 7))
# for length, group in pred_per_len.groupby('peptide_length'):
#     fpr, tpr, _ = roc_curve(group['target'], group['prediction'])
#     roc_auc = auc(fpr, tpr)
#     plot_roc_curve(fpr, tpr, roc_auc, length)

# print("AUC:", roc_auc)


# # Function to plot PR curve for a specific peptide length
# def plot_pr_curve(precision, recall, ap, length):
#     plt.plot(recall, precision, label=f'AP = {ap:.2f} ({length}-mer)')

# # Plot PR curve for each peptide length
# plt.figure(figsize=(7, 7))
# for length, group in pred_per_len.groupby('peptide_length'):
#     precision, recall, _ = precision_recall_curve(group['target'], group['prediction'])
#     ap = average_precision_score(group['target'], group['prediction'])
#     plot_pr_curve(precision, recall, ap, length)

# print("AP:", ap)

# # # Plot formatting
# # plt.title('Precision-Recall Curve')
# # plt.xlabel('Recall')
# # plt.ylabel('Precision')
# # plt.plot([0, 1], [0.5, 0.5], linestyle='--', color='gray', label="Baseline")
# # plt.legend(loc='lower left')
# # plt.tight_layout()
# # plt.savefig("pr_curve.png")
# # plt.show()


# # Add to your peptide-length loop
# for length, group in pred_per_len.groupby('peptide_length'):
#     precision, recall, _ = precision_recall_curve(group['target'], group['prediction'])
#     ap = average_precision_score(group['target'], group['prediction'])

#     # MCC needs binary predictions — apply threshold
#     preds_bin = (group['prediction'] >= BINDER_THRESHOLD).astype(int)
#     mcc = matthews_corrcoef(group['target'], preds_bin)

#     plot_pr_curve(precision, recall, ap, length)


# Load prediction data
pred_data = np.loadtxt(pred_file, dtype=str)
y_true = np.array(pred_data[:, 2], dtype=float)
y_scores = np.array(pred_data[:, 1], dtype=float)

# Threshold for classification
BINDER_THRESHOLD = 0.426
y_bin_true = (y_true >= BINDER_THRESHOLD).astype(int)
y_bin_pred = (y_scores >= BINDER_THRESHOLD).astype(int)

# Compute metrics
auc_score = roc_auc_score(y_bin_true, y_scores)
ap_score = average_precision_score(y_bin_true, y_scores)
mcc_score = matthews_corrcoef(y_bin_true, y_bin_pred)

# Output in CSV format
print(f"{auc_score:.4f},{ap_score:.4f},{mcc_score:.4f}")

