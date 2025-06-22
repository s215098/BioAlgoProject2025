################################################################################
############################### ROC, AUC, T test ###############################
################################################################################


################################ LOADING DATA ##################################

library(ggplot2)
library(dplyr)
library(pROC)


##### Importing data per parameter #####

dir <- "/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/results"
#List with our chosen alleles
alleles <- c("A0201", "A0202", "A1101", "A3001", "B0702", "B1501", "B5401", "B5701")
files <- c("e001.eval", "e002.eval", "e003.eval", "e004.eval", "concat.eval")


# ##### FUNCTION FOR IMPORTING DATA (per model) #####
load_results <- function(folder, model, allele_list, hyperparameters){
  
  for (allele in alleles){
    
    for (file in files) {
      # Construct full path to summary.txt
      file_path <- file.path(dir, model, paste0(allele, ".res"), hyperparameters, file)
      
      # Loading files
      print(paste("Reading:", file_path))
      
      if (file.exists(file_path)) {
        
        data <- read.csv(file_path, sep = " ", header = FALSE)
        
        var_name <- paste0("data_", allele, "_", file)
        assign(var_name, data, envir = .GlobalEnv)
        print(paste("Saved data as", var_name))
      } else {
        warning(paste("File does not exist:", file_path))
      }
    }
  }
}





################################################################################
################################### ANN ########################################
################################################################################


load_results(folder = dir,
             model = "ANN",
             allele_list = alleles,
             hyperparameters = "blosum/hidden.5/epi.0.05"
)



##### Preparing data ######
plot_data_list <- list()

for (file in files) {
  for (allele in alleles) {
    var_name <- paste0("data_", allele, "_", file)
    if (exists(var_name)) {
      df <- get(var_name)
      if (ncol(df) >= 3) {
        temp <- data.frame(
          Measured_con = as.numeric(df[[3]]),
          Predicted_con = as.numeric(df[[2]]),
          Measured_bin = as.integer(df[[5]]),
          Predicted_bin = as.integer(df[[4]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

roc_obj <- roc(plot_df$Measured_bin, plot_df$Predicted_con)
auc_value <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve - ANN")
auc_value



################################################################################
################################## SMM GD ######################################
################################################################################



load_results(folder = dir,
             model = "SMM/gradient_decent",
             allele_list = alleles,
             hyperparameters = "l.0.1"
)



##### Preparing data ######
plot_data_list <- list()

for (file in files) {
  for (allele in alleles) {
    var_name <- paste0("data_", allele, "_", file)
    if (exists(var_name)) {
      df <- get(var_name)
      if (ncol(df) >= 3) {
        temp <- data.frame(
          Measured_con = as.numeric(df[[3]]),
          Predicted_con = as.numeric(df[[2]]),
          Measured_bin = as.integer(df[[5]]),
          Predicted_bin = as.integer(df[[4]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

roc_obj <- roc(plot_df$Measured_bin, plot_df$Predicted_con)
auc_value <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve - SMM GD")
auc_value











################################################################################
################################## SMM MC ######################################
################################################################################



load_results(folder = dir,
             model = "SMM/monte_carlo",
             allele_list = alleles,
             hyperparameters = "l.0.1_t.100"
)



##### Preparing data ######
plot_data_list <- list()

for (file in files) {
  for (allele in alleles) {
    var_name <- paste0("data_", allele, "_", file)
    if (exists(var_name)) {
      df <- get(var_name)
      if (ncol(df) >= 3) {
        temp <- data.frame(
          Measured_con = as.numeric(df[[3]]),
          Predicted_con = as.numeric(df[[2]]),
          Measured_bin = as.integer(df[[5]]),
          Predicted_bin = as.integer(df[[4]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

roc_obj <- roc(plot_df$Measured_bin, plot_df$Predicted_con)
auc_value <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve - SMM MC")
auc_value









################################################################################
#################################### PSSM ######################################
################################################################################


load_results(folder = dir,
             model = "PSSM",
             allele_list = alleles,
             hyperparameters = "b.50"
)


##### Preparing data ######
plot_data_list <- list()

for (file in files) {
  for (allele in alleles) {
    var_name <- paste0("data_", allele, "_", file)
    if (exists(var_name)) {
      df <- get(var_name)
      if (ncol(df) >= 3) {
        temp <- data.frame(
          Measured_con = as.numeric(df[[3]]),
          Predicted_con = as.numeric(df[[2]]),
          Measured_bin = as.integer(df[[5]]),
          Predicted_bin = as.integer(df[[4]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

roc_obj <- roc(plot_df$Measured_bin, plot_df$Predicted_con)
auc_value <- auc(roc_obj)
plot(roc_obj, col = "blue", main = "ROC Curve - PSSM")
auc_value






