{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "ename": "ValueError",
     "evalue": "continuous format is not supported",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mValueError\u001b[0m                                Traceback (most recent call last)",
      "Cell \u001b[0;32mIn[1], line 13\u001b[0m\n\u001b[1;32m     10\u001b[0m y_true \u001b[38;5;241m=\u001b[39m data[\u001b[38;5;241m2\u001b[39m]\n\u001b[1;32m     12\u001b[0m \u001b[38;5;66;03m# Compute ROC curve and AUC\u001b[39;00m\n\u001b[0;32m---> 13\u001b[0m fpr, tpr, thresholds \u001b[38;5;241m=\u001b[39m \u001b[43mroc_curve\u001b[49m\u001b[43m(\u001b[49m\u001b[43my_true\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43my_scores\u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m     14\u001b[0m roc_auc \u001b[38;5;241m=\u001b[39m auc(fpr, tpr)\n\u001b[1;32m     16\u001b[0m \u001b[38;5;66;03m# Plot the ROC curve\u001b[39;00m\n",
      "File \u001b[0;32m/opt/anaconda3/envs/algo/lib/python3.9/site-packages/sklearn/utils/_param_validation.py:216\u001b[0m, in \u001b[0;36mvalidate_params.<locals>.decorator.<locals>.wrapper\u001b[0;34m(*args, **kwargs)\u001b[0m\n\u001b[1;32m    210\u001b[0m \u001b[38;5;28;01mtry\u001b[39;00m:\n\u001b[1;32m    211\u001b[0m     \u001b[38;5;28;01mwith\u001b[39;00m config_context(\n\u001b[1;32m    212\u001b[0m         skip_parameter_validation\u001b[38;5;241m=\u001b[39m(\n\u001b[1;32m    213\u001b[0m             prefer_skip_nested_validation \u001b[38;5;129;01mor\u001b[39;00m global_skip_validation\n\u001b[1;32m    214\u001b[0m         )\n\u001b[1;32m    215\u001b[0m     ):\n\u001b[0;32m--> 216\u001b[0m         \u001b[38;5;28;01mreturn\u001b[39;00m \u001b[43mfunc\u001b[49m\u001b[43m(\u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[43margs\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[38;5;241;43m*\u001b[39;49m\u001b[43mkwargs\u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m    217\u001b[0m \u001b[38;5;28;01mexcept\u001b[39;00m InvalidParameterError \u001b[38;5;28;01mas\u001b[39;00m e:\n\u001b[1;32m    218\u001b[0m     \u001b[38;5;66;03m# When the function is just a wrapper around an estimator, we allow\u001b[39;00m\n\u001b[1;32m    219\u001b[0m     \u001b[38;5;66;03m# the function to delegate validation to the estimator, but we replace\u001b[39;00m\n\u001b[1;32m    220\u001b[0m     \u001b[38;5;66;03m# the name of the estimator by the name of the function in the error\u001b[39;00m\n\u001b[1;32m    221\u001b[0m     \u001b[38;5;66;03m# message to avoid confusion.\u001b[39;00m\n\u001b[1;32m    222\u001b[0m     msg \u001b[38;5;241m=\u001b[39m re\u001b[38;5;241m.\u001b[39msub(\n\u001b[1;32m    223\u001b[0m         \u001b[38;5;124mr\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mparameter of \u001b[39m\u001b[38;5;124m\\\u001b[39m\u001b[38;5;124mw+ must be\u001b[39m\u001b[38;5;124m\"\u001b[39m,\n\u001b[1;32m    224\u001b[0m         \u001b[38;5;124mf\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mparameter of \u001b[39m\u001b[38;5;132;01m{\u001b[39;00mfunc\u001b[38;5;241m.\u001b[39m\u001b[38;5;18m__qualname__\u001b[39m\u001b[38;5;132;01m}\u001b[39;00m\u001b[38;5;124m must be\u001b[39m\u001b[38;5;124m\"\u001b[39m,\n\u001b[1;32m    225\u001b[0m         \u001b[38;5;28mstr\u001b[39m(e),\n\u001b[1;32m    226\u001b[0m     )\n",
      "File \u001b[0;32m/opt/anaconda3/envs/algo/lib/python3.9/site-packages/sklearn/metrics/_ranking.py:1150\u001b[0m, in \u001b[0;36mroc_curve\u001b[0;34m(y_true, y_score, pos_label, sample_weight, drop_intermediate)\u001b[0m\n\u001b[1;32m   1046\u001b[0m \u001b[38;5;129m@validate_params\u001b[39m(\n\u001b[1;32m   1047\u001b[0m     {\n\u001b[1;32m   1048\u001b[0m         \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124my_true\u001b[39m\u001b[38;5;124m\"\u001b[39m: [\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124marray-like\u001b[39m\u001b[38;5;124m\"\u001b[39m],\n\u001b[0;32m   (...)\u001b[0m\n\u001b[1;32m   1057\u001b[0m     y_true, y_score, \u001b[38;5;241m*\u001b[39m, pos_label\u001b[38;5;241m=\u001b[39m\u001b[38;5;28;01mNone\u001b[39;00m, sample_weight\u001b[38;5;241m=\u001b[39m\u001b[38;5;28;01mNone\u001b[39;00m, drop_intermediate\u001b[38;5;241m=\u001b[39m\u001b[38;5;28;01mTrue\u001b[39;00m\n\u001b[1;32m   1058\u001b[0m ):\n\u001b[1;32m   1059\u001b[0m \u001b[38;5;250m    \u001b[39m\u001b[38;5;124;03m\"\"\"Compute Receiver operating characteristic (ROC).\u001b[39;00m\n\u001b[1;32m   1060\u001b[0m \n\u001b[1;32m   1061\u001b[0m \u001b[38;5;124;03m    Note: this implementation is restricted to the binary classification task.\u001b[39;00m\n\u001b[0;32m   (...)\u001b[0m\n\u001b[1;32m   1148\u001b[0m \u001b[38;5;124;03m    array([ inf, 0.8 , 0.4 , 0.35, 0.1 ])\u001b[39;00m\n\u001b[1;32m   1149\u001b[0m \u001b[38;5;124;03m    \"\"\"\u001b[39;00m\n\u001b[0;32m-> 1150\u001b[0m     fps, tps, thresholds \u001b[38;5;241m=\u001b[39m \u001b[43m_binary_clf_curve\u001b[49m\u001b[43m(\u001b[49m\n\u001b[1;32m   1151\u001b[0m \u001b[43m        \u001b[49m\u001b[43my_true\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43my_score\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43mpos_label\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43mpos_label\u001b[49m\u001b[43m,\u001b[49m\u001b[43m \u001b[49m\u001b[43msample_weight\u001b[49m\u001b[38;5;241;43m=\u001b[39;49m\u001b[43msample_weight\u001b[49m\n\u001b[1;32m   1152\u001b[0m \u001b[43m    \u001b[49m\u001b[43m)\u001b[49m\n\u001b[1;32m   1154\u001b[0m     \u001b[38;5;66;03m# Attempt to drop thresholds corresponding to points in between and\u001b[39;00m\n\u001b[1;32m   1155\u001b[0m     \u001b[38;5;66;03m# collinear with other points. These are always suboptimal and do not\u001b[39;00m\n\u001b[1;32m   1156\u001b[0m     \u001b[38;5;66;03m# appear on a plotted ROC curve (and thus do not affect the AUC).\u001b[39;00m\n\u001b[0;32m   (...)\u001b[0m\n\u001b[1;32m   1161\u001b[0m     \u001b[38;5;66;03m# but does not drop more complicated cases like fps = [1, 3, 7],\u001b[39;00m\n\u001b[1;32m   1162\u001b[0m     \u001b[38;5;66;03m# tps = [1, 2, 4]; there is no harm in keeping too many thresholds.\u001b[39;00m\n\u001b[1;32m   1163\u001b[0m     \u001b[38;5;28;01mif\u001b[39;00m drop_intermediate \u001b[38;5;129;01mand\u001b[39;00m \u001b[38;5;28mlen\u001b[39m(fps) \u001b[38;5;241m>\u001b[39m \u001b[38;5;241m2\u001b[39m:\n",
      "File \u001b[0;32m/opt/anaconda3/envs/algo/lib/python3.9/site-packages/sklearn/metrics/_ranking.py:818\u001b[0m, in \u001b[0;36m_binary_clf_curve\u001b[0;34m(y_true, y_score, pos_label, sample_weight)\u001b[0m\n\u001b[1;32m    816\u001b[0m y_type \u001b[38;5;241m=\u001b[39m type_of_target(y_true, input_name\u001b[38;5;241m=\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;124my_true\u001b[39m\u001b[38;5;124m\"\u001b[39m)\n\u001b[1;32m    817\u001b[0m \u001b[38;5;28;01mif\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m (y_type \u001b[38;5;241m==\u001b[39m \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mbinary\u001b[39m\u001b[38;5;124m\"\u001b[39m \u001b[38;5;129;01mor\u001b[39;00m (y_type \u001b[38;5;241m==\u001b[39m \u001b[38;5;124m\"\u001b[39m\u001b[38;5;124mmulticlass\u001b[39m\u001b[38;5;124m\"\u001b[39m \u001b[38;5;129;01mand\u001b[39;00m pos_label \u001b[38;5;129;01mis\u001b[39;00m \u001b[38;5;129;01mnot\u001b[39;00m \u001b[38;5;28;01mNone\u001b[39;00m)):\n\u001b[0;32m--> 818\u001b[0m     \u001b[38;5;28;01mraise\u001b[39;00m \u001b[38;5;167;01mValueError\u001b[39;00m(\u001b[38;5;124m\"\u001b[39m\u001b[38;5;132;01m{0}\u001b[39;00m\u001b[38;5;124m format is not supported\u001b[39m\u001b[38;5;124m\"\u001b[39m\u001b[38;5;241m.\u001b[39mformat(y_type))\n\u001b[1;32m    820\u001b[0m check_consistent_length(y_true, y_score, sample_weight)\n\u001b[1;32m    821\u001b[0m y_true \u001b[38;5;241m=\u001b[39m column_or_1d(y_true)\n",
      "\u001b[0;31mValueError\u001b[0m: continuous format is not supported"
     ]
    }
   ],
   "source": [
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "from sklearn.metrics import roc_curve, auc\n",
    "\n",
    "# Load the data (no header, whitespace-delimited)\n",
    "data = pd.read_csv(\"/Users/alberteenglund/Documents/DTU/8_Semester/22125_algorithms_in_bioinformatics/algorithms/BioAlgoProject2025/results/ANN/A0201.res/blosum/hidden.5/epi.0.05/e001.eval\", sep=\"\\s+\", header=None)\n",
    "\n",
    "# Column 2 = predicted scores, Column 3 = true labels\n",
    "y_scores = data[1]\n",
    "y_true = data[2]\n",
    "\n",
    "# Compute ROC curve and AUC\n",
    "fpr, tpr, thresholds = roc_curve(y_true, y_scores)\n",
    "roc_auc = auc(fpr, tpr)\n",
    "\n",
    "# Plot the ROC curve\n",
    "plt.figure(figsize=(7, 6))\n",
    "plt.plot(fpr, tpr, color=\"darkorange\", lw=2, label=f\"ROC curve (AUC = {roc_auc:.3f})\")\n",
    "plt.plot([0, 1], [0, 1], color=\"gray\", lw=1, linestyle=\"--\")\n",
    "plt.xlim([0.0, 1.0])\n",
    "plt.ylim([0.0, 1.05])\n",
    "plt.xlabel(\"False Positive Rate\")\n",
    "plt.ylabel(\"True Positive Rate\")\n",
    "plt.title(\"ROC Curve\")\n",
    "plt.legend(loc=\"lower right\")\n",
    "plt.grid(True)\n",
    "plt.tight_layout()\n",
    "plt.show()\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "algo",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.13"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
