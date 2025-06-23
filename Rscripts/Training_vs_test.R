
#### Training performance vs test performance ####

################################################################################
##################### FUNCTION FOR IMPORTING DATA ##############################
################################################################################


### Importing libraries:
library(ggplot2)
library(dplyr)
library(tidyr)


##### Importing data per parameter #####

dir <- "/Users/kristinetoftjohansen/Desktop/Algo/BioAlgoProject2025/results"
#List with our chosen alleles
alleles <- c("A0201", "A0202", "A1101", "A3001", "B0702", "B1501", "B5401", "B5701")
files <- c("train.log.1", "train.log.2", "train.log.3", "train.log.4")
# files <- c("train.log.4")

# ##### FUNCTION FOR IMPORTING DATA (per model) #####
load_results <- function(folder, model, allele_list, hyperparameters){
  
  # files <- c("train.log.1", "train.log.2", "train.log.3", "train.log.4")
  
  for (allele in alleles){
    
    for (file in files) {
      # Construct full path to summary.txt
      file_path <- file.path(dir, model, paste0(allele, ".res"), hyperparameters, file)
      
      # Loading files
      print(paste("Reading:", file_path))
      
      if (file.exists(file_path)) {
        
        # read all lines
        all_lines <- readLines(file_path)
        
        # filter away the hashtagged lines:
        clean_lines <- all_lines[!grepl("^#", all_lines)]
        clean_lines <- clean_lines[-(1:7)]
        
        # Split each line on whitespace
        split_lines <- strsplit(clean_lines, "\\s+")
        
        # Convert to data frame
        df <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE) %>%
          select(V5, V9, V13, V17) %>%
          mutate(Allele = allele)
        # and assign column names:
        colnames(df) <- c("Train error", "Train perf", "Test error", "Test perf", "Allele")
        
        # add x meaning iterations:
        df$x <- seq_len(nrow(df))
        
        var_name <- paste0("data_", allele, "_", file)
        assign(var_name, df, envir = .GlobalEnv)
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



# Collect data for all alleles and train logs into a list
data_list <- list()

for (allele in alleles) {
  for (file in files) {
    var_name <- paste0("data_", allele, "_", file)
    if (exists(var_name)) {
      df <- get(var_name)
      df$File <- file  # Add which training run this is
      data_list[[length(data_list) + 1]] <- df
    }
  }
}

# Combine all data frames
all_data <- bind_rows(data_list)




### Convert to long format for plotting train/test errors together

# For errors 
long_data_error <- all_data %>%
  pivot_longer(cols = c("Train error", "Test error"),
               names_to = "Error_type",
               values_to = "Error") %>%
  mutate(Error = as.numeric(Error))  # Make sure Error column is numeric


# For performance
long_data_perf <- all_data %>%
  pivot_longer(cols = c("Train perf", "Test perf"),
               names_to = "Perf_type",
               values_to = "Performance") %>%
  mutate(Performance = as.numeric(Performance))  # Make sure Error column is numeric




### Plotting
library(patchwork)



# for plotting each individual train.log - remember to redefine in the start of the script.
p_error <- ggplot(long_data_error, aes(x = x, y = Error, color = Error_type)) +
  geom_line() +
  facet_wrap(~ Allele) +
  labs(title = "Error", x = "Iteration", y = "Error") +
  theme_minimal()

p_perf <- ggplot(long_data_perf, aes(x = x, y = Performance, color = Perf_type)) +
  geom_line() +
  facet_wrap(~ Allele) +
  labs(title = "Performance", x = "Iteration", y = "Performance") +
  theme_minimal()

# p_error
# p_perf 



# for plotting all together:
p_error_all <- ggplot(long_data_error, aes(x = x, y = Error, color = File, linetype = Error_type)) +
  geom_line() +
  facet_wrap(~ Allele) +
  labs(title = "Error over Iterations",
       x = "Iteration",
       y = "Error",
       linetype = "Dataset",
       color = "Model") +
  
  scale_linetype_manual(Dataset = c("Train error" = "dotted", "Test error" = "solid"),
                        Model = c("Model 1", "Model 2", "Model 3", "Model 4")) +
  theme_minimal()

p_error_all


# laver nye legends



ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/ANN/Training_vs_test_error_ANN_all.png"), width = 12, height = 9)


p_perf_all <- ggplot(long_data_perf, aes(x = x, y = Performance, color = File, linetype = Perf_type)) +
  geom_line() +
  facet_wrap(~ Allele) +
  labs(title = "Performance over Iterations", x = "Iteration", y = "Performance") +
  scale_linetype_manual(values = c("Train perf" = "dotted", "Test perf" = "solid")) +
  theme_minimal()
p_perf_all
ggsave(paste0("~/Desktop/Algo/BioAlgoProject2025/Plots/ANN/Training_vs_test_perf_ANN_all.png"), width = 12, height = 9)







##### NYE INDIVIDUELLE PLOTS ######



# Choose two alleles to plot
selected_alleles <- c("A0202")

# Filter for just these two
filtered_error <- long_data_error %>%
  filter(Allele %in% selected_alleles)

filtered_perf <- long_data_perf %>%
  filter(Allele %in% selected_alleles)

p_perf_selected <- ggplot(filtered_perf,
                          aes(x = x,
                              y = Performance,
                              color = File,
                              linetype = Perf_type)) +
  geom_line() +
  facet_wrap(~ Allele) +
  labs(title = "Performance - training vs test",
       x = "Iteration",
       y = "Performance (PCC)",
       color = "Model",
       linetype = "Dataset") +
  scale_linetype_manual(values = c("Train perf" = "dotted", "Test perf" = "solid")) +
  scale_color_manual(values = c("train.log.1" = "#8BCDF9",
                                "train.log.2" = "#FC7322",
                                "train.log.3" = "#478978",
                                "train.log.4" = "#171749"),
                     labels = c("train.log.1" = "Model 1",
                                "train.log.2" = "Model 2",
                                "train.log.3" = "Model 3",
                                "train.log.4" = "Model 4")) +
  theme_minimal()
p_perf_selected


