################################################################################
################################## PCC PLOT ####################################
################################################################################

# Plots predicted values against measured values
# and lines with the models
# in "eval.000" (for all four + for concatenated)

################################ LOADING DATA ##################################

library(ggplot2)
library(dplyr)




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
          Measured = as.numeric(df[[3]]),
          Predicted = as.numeric(df[[2]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

# Clean model labels (optional but tidy)
plot_df$Model <- factor(plot_df$Model, levels = files,
                        labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated"))
# Separate datasets
points_df <- plot_df[plot_df$Model != "Concatenated", ] #don't plot points for the concatenated as it will overwrite the others. 
lines_df  <- plot_df  # use full data for lines


# Calculate correlation. 
pcc_labels <- lines_df %>%
  group_by(Allele, Model) %>%
  summarise(pcc = cor(Predicted, Measured), .groups = "drop") %>%
  mutate(label = sprintf("PCC = %.2f", pcc))

# Add y-offsets to avoid overlap
# model_levels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated")
model_levels <- c("Concatenated", "Model 4", "Model 3", "Model 2", "Model 1")
pcc_labels$Model <- factor(pcc_labels$Model, levels = model_levels) # Set Model as factor with levels in pcc_labels

pcc_labels <- pcc_labels %>%
  group_by(Allele) %>%
  arrange(Model) %>%  # orders by factor levels, not alphabetically
  mutate(y = 0.01 + 0.05 * row_number(),  # space out lines top-down
         x = 0.9)                         # fixed x-position (adjust if needed)

#### PLOTTING ####

ggplot() +
  geom_point(data = points_df, aes(x = Predicted, y = Measured, color = Model), alpha = 0.6) +
  geom_smooth(data = lines_df, aes(x = Predicted, y = Measured, color = Model), method = "lm", se = FALSE) +
  # geom_text(data = pcc_labels, aes(x = x, y = y, label = label, color = Model), hjust = 0, size = 3) +
  geom_label(
    data = pcc_labels,
    aes(x = x, y = y, label = label, color = Model),
    fill = "white",      # background color of label
    label.size = 0,    # border thickness around label box
    # alpha = 0.8,          # transparency of label background
    show.legend = FALSE
  ) +
  # facet_wrap(~ Allele, scales = "free") +
  facet_wrap(~ Allele) +
  labs(
    title = "Predicted vs Measured Values per Allele",
    subtitle = "ANN",
    x = "Predicted",
    y = "Measured",
    color = "Model"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/ANN/PCC_plot_ANN.png"), width = 16, height = 14)








################################################################################
################################### SMM GD #####################################
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
          Measured = as.numeric(df[[3]]),
          Predicted = as.numeric(df[[2]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

# Clean model labels (optional but tidy)
plot_df$Model <- factor(plot_df$Model, levels = files,
                        labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated"))
# Separate datasets
points_df <- plot_df[plot_df$Model != "Concatenated", ] #don't plot points for the concatenated as it will overwrite the others. 
lines_df  <- plot_df  # use full data for lines


# Calculate correlation. 
pcc_labels <- lines_df %>%
  group_by(Allele, Model) %>%
  summarise(pcc = cor(Predicted, Measured), .groups = "drop") %>%
  mutate(label = sprintf("PCC = %.2f", pcc))

# Add y-offsets to avoid overlap
# model_levels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated")
model_levels <- c("Concatenated", "Model 4", "Model 3", "Model 2", "Model 1")
pcc_labels$Model <- factor(pcc_labels$Model, levels = model_levels) # Set Model as factor with levels in pcc_labels

pcc_labels <- pcc_labels %>%
  group_by(Allele) %>%
  arrange(Model) %>%  # orders by factor levels, not alphabetically
  mutate(y = -0.4 + 0.08 * row_number(),  # space out lines top-down
         x = 0.99)                         # fixed x-position (adjust if needed)

#### PLOTTING ####

ggplot() +
  geom_point(data = points_df, aes(x = Predicted, y = Measured, color = Model), alpha = 0.6) +
  geom_smooth(data = lines_df, aes(x = Predicted, y = Measured, color = Model), method = "lm", se = FALSE) +
  # geom_text(data = pcc_labels, aes(x = x, y = y, label = label, color = Model), hjust = 0, size = 3) +
  geom_label(
    data = pcc_labels,
    aes(x = x, y = y, label = label, color = Model),
    fill = "white",      # background color of label
    label.size = 0,    # border thickness around label box
    # alpha = 0.8,          # transparency of label background
    show.legend = FALSE
  ) +
  # facet_wrap(~ Allele, scales = "free") +
  facet_wrap(~ Allele) +
  labs(
    title = "Predicted vs Measured Values per Allele",
    subtitle = "SMM Gradient Decent",
    x = "Predicted",
    y = "Measured",
    color = "Model"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/SMM_GD/PCC_plot_SMM_GD.png"), width = 16, height = 14)









################################################################################
################################### SMM MC #####################################
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
          Measured = as.numeric(df[[3]]),
          Predicted = as.numeric(df[[2]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

# Clean model labels (optional but tidy)
plot_df$Model <- factor(plot_df$Model, levels = files,
                        labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated"))
# Separate datasets
points_df <- plot_df[plot_df$Model != "Concatenated", ] #don't plot points for the concatenated as it will overwrite the others. 
lines_df  <- plot_df  # use full data for lines


# Calculate correlation. 
pcc_labels <- lines_df %>%
  group_by(Allele, Model) %>%
  summarise(pcc = cor(Predicted, Measured), .groups = "drop") %>%
  mutate(label = sprintf("PCC = %.2f", pcc))

# Add y-offsets to avoid overlap
# model_levels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated")
model_levels <- c("Concatenated", "Model 4", "Model 3", "Model 2", "Model 1")
pcc_labels$Model <- factor(pcc_labels$Model, levels = model_levels) # Set Model as factor with levels in pcc_labels

pcc_labels <- pcc_labels %>%
  group_by(Allele) %>%
  arrange(Model) %>%  # orders by factor levels, not alphabetically
  mutate(y = -0.6 + 0.09 * row_number(),  # space out lines top-down
         x = 0.5)                         # fixed x-position (adjust if needed)

#### PLOTTING ####

ggplot() +
  geom_point(data = points_df, aes(x = Predicted, y = Measured, color = Model), alpha = 0.6) +
  geom_smooth(data = lines_df, aes(x = Predicted, y = Measured, color = Model), method = "lm", se = FALSE) +
  # geom_text(data = pcc_labels, aes(x = x, y = y, label = label, color = Model), hjust = 0, size = 3) +
  geom_label(
    data = pcc_labels,
    aes(x = x, y = y, label = label, color = Model),
    fill = "white",      # background color of label
    label.size = 0,    # border thickness around label box
    # alpha = 0.8,          # transparency of label background
    show.legend = FALSE
  ) +
  # facet_wrap(~ Allele) +
  facet_wrap(~ Allele, scale = "free_x") +
  labs(
    title = "Predicted vs Measured Values per Allele",
    subtitle = "SMM Monte Carlo",
    x = "Predicted",
    y = "Measured",
    color = "Model"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))


ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/SMM_MC/PCC_plot_SMM_MC.png"), width = 16, height = 14)








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
          Measured = as.numeric(df[[3]]),
          Predicted = as.numeric(df[[2]]),
          Allele = allele,
          Model = file
        )
        plot_data_list[[length(plot_data_list) + 1]] <- temp
      }
    }
  }
}

plot_df <- do.call(rbind, plot_data_list)

# Clean model labels (optional but tidy)
plot_df$Model <- factor(plot_df$Model, levels = files,
                        labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated"))
# Separate datasets
points_df <- plot_df[plot_df$Model != "Concatenated", ] #don't plot points for the concatenated as it will overwrite the others. 
lines_df  <- plot_df  # use full data for lines


# Calculate correlation. 
pcc_labels <- lines_df %>%
  group_by(Allele, Model) %>%
  summarise(pcc = cor(Predicted, Measured), .groups = "drop") %>%
  mutate(label = sprintf("PCC = %.2f", pcc))

# Add y-offsets to avoid overlap
# model_levels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Concatenated")
model_levels <- c("Concatenated", "Model 4", "Model 3", "Model 2", "Model 1")
pcc_labels$Model <- factor(pcc_labels$Model, levels = model_levels) # Set Model as factor with levels in pcc_labels

pcc_labels <- pcc_labels %>%
  group_by(Allele) %>%
  arrange(Model) %>%  # orders by factor levels, not alphabetically
  mutate(y = 0.01 + 0.05 * row_number(),  # space out lines top-down
         x = 25)                         # fixed x-position (adjust if needed)

#### PLOTTING ####

ggplot() +
  geom_point(data = points_df, aes(x = Predicted, y = Measured, color = Model), alpha = 0.6) +
  geom_smooth(data = lines_df, aes(x = Predicted, y = Measured, color = Model), method = "lm", se = FALSE) +
  # geom_text(data = pcc_labels, aes(x = x, y = y, label = label, color = Model), hjust = 0, size = 3) +
  geom_label(
    data = pcc_labels,
    aes(x = x, y = y, label = label, color = Model),
    fill = "white",      # background color of label
    label.size = 0,    # border thickness around label box
    # alpha = 0.8,          # transparency of label background
    show.legend = FALSE
  ) +
  # facet_wrap(~ Allele, scales = "free") +
  facet_wrap(~ Allele) +
  labs(
    title = "Predicted vs Measured Values per Allele",
    subtitle = "PSSM",
    x = "Predicted",
    y = "Measured",
    color = "Model"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))



ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/PSSM/PCC_plot_PSSM.png"), width = 12, height = 9)

