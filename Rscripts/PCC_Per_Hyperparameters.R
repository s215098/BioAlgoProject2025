##### IMPORTING LIBRARIES ####
library(ggplot2)
library(dplyr)

# ##### IMPORTING DATA (per parameter) #####
results_dir <- "/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/results"

# ##### FUNCTION FOR IMPORTING DATA (per allele) #####
load_results <- function(folder, allele_list){
  for (allele in allele_list){
    # Construct full path to summary.txt
    file_path <- file.path(results_dir, folder, paste0(allele, ".res"), "summary.txt")
    
    print(paste("Reading:", file_path))
    
    if (file.exists(file_path)) {
      data <- read.csv(file_path, sep = " ", header = FALSE)
      var_name <- paste0("data_", allele)
      assign(var_name, data, envir = .GlobalEnv)
      print(paste("Saved data as", var_name))
    } else {
      warning(paste("File does not exist:", file_path))
    }
  }
}

#List with our chosen alleles
alleles <- c("A0201", "A0202", "A1101", "A3001", "B0702", "B1501", "B5401", "B5701")


######### PLotting PCCs for the different hyperparameters. #########
# # Adding column names.
# colnames(data) = c("Hyperparameter", "PCC", "MSE")
# 
# # Making a plot.
# ggplot(data = data, aes(x = Hyperparameter, y=PCC)) + 
#   geom_point() +
#   labs(title = "PCC vs Hyperparameters",
#        x = "Hyperparameters",
#        y = "PCC")


######### BAR PLOT PEARSON CORRELATION #########
# To visualize which hyperparameter performed best throughout the alleles. 


### For PSSM ###
load_results("PSSM", allele_list = alleles)

# Combine data frames and include Allele as a column
data_list <- lapply(alleles, function(a) {
  df <- get(paste0("data_", a))
  return(df)
})

# Combine all into one tidy data frame
all_data <- do.call(rbind, data_list)

#subsetting
all_data_PSSM <- do.call(rbind, data_list) %>% select(V1, V2, V4, V6, V8)
colnames(all_data_PSSM) <- c("Allele", "Hyperparameter", "N", "PCC", "MSE")

all_data_PSSM$Hyperparameter <- as.factor(all_data_PSSM$Hyperparameter)

#Making the plot:
ggplot(all_data_PSSM, aes(x = Hyperparameter, y = PCC)) +
  geom_bar(stat = "identity", fill = "#8BCDF9") +
  facet_wrap(~ Allele) +
  theme_minimal() +
  labs(
    title = "PCC values across hyperparameter thresholds - PSSM",
    x = "Beta",
    y = "PCC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Saving the plot with width = 1600 and height = 900
ggsave("~/Desktop/Algo/BioAlgoProject2025/Plots/PSSM/PCC_vs_Hyperparameter_PSSM_beta.png",
       device = "png",
       width = 10,
       height = 8
       )




### For SMM Gradient Decent ###

load_results("SMM/gradient_decent", allele_list = alleles)

# Combine data frames and include Allele as a column
data_list <- lapply(alleles, function(a) {
  df <- get(paste0("data_", a))
  return(df)
})

# Combine all into one tidy data frame
all_data <- do.call(rbind, data_list)


all_data_SMMGD <- do.call(rbind, data_list) %>% select(V1, V2, V4, V6, V8)
colnames(all_data_SMMGD) <- c("Allele", "Hyperparameter", "N", "PCC", "MSE")

all_data_SMMGD$Hyperparameter <- as.factor(all_data_SMMGD$Hyperparameter)

#Making the plot:
ggplot(all_data_SMMGD, aes(x = Hyperparameter, y = PCC)) +
  geom_bar(stat = "identity", fill = "#8BCDF9") +
  facet_wrap(~ Allele) +
  theme_minimal() +
  labs(
    title = "PCC values across hyperparameter thresholds - SMM GD",
    x = "Lambda",
    y = "PCC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Saving the plot with width = 1600 and height = 900
ggsave("~/Desktop/Algo/BioAlgoProject2025/Plots/SMM_GD/PCC_vs_Hyperparameter_SMMGD_Lambda.png",
       device = "png",
       width = 10,
       height = 8
)





### For SMM Monte Carlo ###

load_results("SMM/monte_carlo", allele_list = alleles)

# Combine data frames and include Allele as a column
data_list <- lapply(alleles, function(a) {
  df <- get(paste0("data_", a))
  return(df)
})

# Combine all into one tidy data frame
all_data <- do.call(rbind, data_list)


all_data_SMMMC <- do.call(rbind, data_list) %>%
  select(V1, V2, V3, V5, V7, V9) %>%
  mutate(Combined = paste(V2, V3, sep = " & ")) %>% #combine to new column
  relocate(Combined, .before = V2) %>% #move before the original columns
  select(-V2, -V3) # remove original columns

colnames(all_data_SMMMC) <- c("Allele", "Hyperparameter", "N", "PCC", "MSE")

all_data_SMMMC$Hyperparameter <- as.factor(all_data_SMMMC$Hyperparameter)

#Making the plot:
ggplot(all_data_SMMMC, aes(x = Hyperparameter, y = PCC)) +
  geom_bar(stat = "identity", fill = "#8BCDF9") +
  facet_wrap(~ Allele) +
  theme_minimal() +
  labs(
    title = "PCC values across hyperparameter thresholds - SMM MC",
    x = "Lambda & Tsteps",
    y = "PCC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Saving the plot with width = 1600 and height = 900
ggsave("~/Desktop/Algo/BioAlgoProject2025/Plots/SMM_MC/PCC_vs_Hyperparameter_SMMMC_Lambda&Tsteps.png",
       device = "png",
       width = 10,
       height = 8
)









### For ANN ###

load_results("ANN", allele_list = alleles)

# Combine data frames and include Allele as a column
data_list <- lapply(alleles, function(a) {
  df <- get(paste0("data_", a))
  return(df)
})

# Combine all into one tidy data frame
all_data <- do.call(rbind, data_list)


all_data_ANN <- do.call(rbind, data_list) %>%
  select(V1, V2, V5, V7, V9, V11, V13) %>%
  mutate(Combined = paste(V2, V5, V7, sep = " & ")) %>% #combine to new column
  relocate(Combined, .before = V2) %>% #move before the original columns
  select(-V2, -V5, -V7) # remove original columns

colnames(all_data_ANN) <- c("Allele", "Hyperparameter", "N", "PCC", "MSE")

all_data_ANN$Hyperparameter <- as.factor(all_data_ANN$Hyperparameter)

#Making the plot:
ggplot(all_data_ANN, aes(x = Hyperparameter, y = PCC)) +
  geom_bar(stat = "identity", fill = "#8BCDF9") +
  # geom_label()
  facet_wrap(~ Allele) +
  theme_minimal() +
  labs(
    title = "PCC values across hyperparameter thresholds - ANN",
    x = "Encoding (sparse / blosum) & N Hidden Neurons & Epsilon",
    y = "PCC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

# Saving the plot with width = 1600 and height = 900
ggsave("~/Desktop/Algo/BioAlgoProject2025/Plots/ANN/PCC_vs_Hyperparameter_SMMMC_Lambda&Tsteps.png",
       device = "png",
       width = 10,
       height = 8
)

