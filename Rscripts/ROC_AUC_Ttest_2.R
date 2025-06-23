################################################################################
############################### ROC, AUC, T test ###############################
################################################################################

################################ LOADING DATA ##################################

library(ggplot2)
library(dplyr)
library(pROC)

##### PARAMETERS #####
dir <- "/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/results"
alleles <- c("A0201", "A0202", "A1101", "A3001", "B0702", "B1501", "B5401", "B5701")
files <- c("e001.eval", "e002.eval", "e003.eval", "e004.eval", "concat.eval")

##### FUNCTION TO LOAD DATA #####
load_results <- function(folder, model, allele_list, hyperparameters) {
  data_list <- list()
  for (allele in allele_list) {
    for (file in files) {
      file_path <- file.path(folder, model, paste0(allele, ".res"), hyperparameters, file)
      if (file.exists(file_path)) {
        df <- read.csv(file_path, sep = " ", header = FALSE)
        if (ncol(df) >= 5) {
          temp <- data.frame(
            Measured_con = as.numeric(df[[3]]),
            Predicted_con = as.numeric(df[[2]]),
            Measured_bin = as.integer(df[[5]]),
            Predicted_bin = as.integer(df[[4]]),
            Allele = allele,
            Model = file
          )
          data_list[[length(data_list) + 1]] <- temp
        }
      } else {
        warning(paste("File does not exist:", file_path))
      }
    }
  }
  return(do.call(rbind, data_list))
}

##### FUNCTION TO PLOT ROC AND CALCULATE AUC #####
plot_roc_auc <- function(data, model_name) {
  roc_obj <- roc(data$Measured_bin, data$Predicted_con)
  auc_value <- auc(roc_obj)
  plot(roc_obj, col = "blue", main = paste("ROC Curve -", model_name))
  print(paste("AUC for", model_name, ":", round(auc_value, 4)))
  return(auc_value)
}

################################################################################
################################### RUN MODELS #################################
################################################################################

# ANN
data_ann <- load_results(dir, "ANN", alleles, "blosum/hidden.5/epi.0.05")
auc_ann <- plot_roc_auc(data_ann, "ANN")

# SMM - Gradient Descent
data_smm_gd <- load_results(dir, "SMM/gradient_decent", alleles, "l.0.1")
auc_smm_gd <- plot_roc_auc(data_smm_gd, "SMM GD")

# SMM - Monte Carlo
data_smm_mc <- load_results(dir, "SMM/monte_carlo", alleles, "l.0.1_t.100")
auc_smm_mc <- plot_roc_auc(data_smm_mc, "SMM MC")

# PSSM
data_pssm <- load_results(dir, "PSSM", alleles, "b.50")
auc_pssm <- plot_roc_auc(data_pssm, "PSSM")
