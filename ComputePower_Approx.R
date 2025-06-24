# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(parallel,pbapply,tidyverse,progress,matrixStats) 

# Set the number of cores to use
numCores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(numCores)     # Create a cluster
scenario <- 1
subscenario <- "Alpha1T1Alpha2T2"

# Load your required scripts/functions on each cluster node
clusterEvalQ(cl, {
  source("LogOR_Function.R")
  source("PowerFunction_SumBeta.R")
  # Treatment Preference Proportions
  pP1=0.5 #  desired proportion of individuals expressing Preference in stage 1
  pP2=0.5 # desired proportion of patients expressing Preference in stage 2 (among non-responders)
  
  # prior shape parameters
  response_rate_prior_shapes <- rep(1/3, 12)
  
  # linkage parameters
  #linkage_parameter_values <- rep(1, 6) # S1
  linkage_parameter_values <- c(1, 1, 1.1, 1.02, 1.05, 1.2) # S2
  
  
  # indifference response rates
  pi_A <- 0.6                                          
  pi_B <- 0.45  
  pi_AC <- 0.5                                          
  pi_AD <- 0.4                                         
  pi_BC <- 0.3                                         
  pi_BD <- 0.2   
  
  # Specify theta targets
  pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
  pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)
  
  # List of sample sizes to iterate over
  sample_sizes <- seq(300, 1000, 100)
  
})


# Set seed for reproducibility in each cluster
clusterSetRNGStream(cl, 12)

# Define the function for power calculation
compute_power_wrapper <- function(N_samp) {
  powerFunction(sample_size=N_samp, 
                Treatment_Preference_stage1_prob = pP1,
                Treatment_Preference_stage2_prob = pP2, 
                PreferA_prob = pTheta_A, 
                PreferC_prob = pTheta_C,
                stage2_NR_outcome_prob = c(pi_AC, pi_AD, pi_BC, pi_BD),
                stage1_trtA_outcome_prob = pi_A,
                stage1_trtB_outcome_prob = pi_B,
                linkage_params = linkage_parameter_values,
                prior_shape_params = response_rate_prior_shapes,
                threshold = 0.2,
                alpha_type1 = 0.05)  # Pass sample size as an argument
}

# Run in parallel with a progress bar
sample_sizes <- seq(300, 1000, 100)

start_time <- Sys.time()
power_results <- pblapply(sample_sizes, compute_power_wrapper, cl = cl)
end_time <- Sys.time()

# Stop the cluster
stopCluster(cl)
closeAllConnections()

# Assuming power_results is structured as follows:
# power_results[[i]] is a list of the form: list(power = ..., other1 = ..., other2 = ...)

# Convert the list to a data frame
# Extract the first component (power) and sample size
power_df <- data.frame(
  N = sample_sizes,  # Sample sizes
  Power = sapply(power_results, function(res) res[[1]])  # Extract the power from each list
)

# For example, if you want to include other components, assuming they are at index 2 and 3:
power_df$SetBestCorrect <- sapply(power_results, function(res) res[[2]])
power_df$BestIncluded <- sapply(power_results, function(res) res[[3]])
power_df$ComputationTime <- end_time - start_time


# View the resulting data frame
print(power_df)

power_df$ComputationTime <- round(power_df$ComputationTime,2)
# gt(power_df[,c(1,2,4,5)]) %>% cols_label(ComputationTime = "Computation Time (Sec)",
#                                          BestIncluded = "Pr(Best Included in Set)")

### Save Settings
# set date and time for file saving 
st<-format(Sys.time(), "%Y_%m_%d_%H_%M")

# Define the folder path where you want your results saved
# folder_path <- "path to your folder where you want to store results"

# Example: 
folder_path <- paste0("/Users/mariwank/Documents/UMICH/dissertation/SampleSize/FinalSimResults/SumBeta/Scenario", scenario, "/", subscenario)

# Create the folder if it doesn't exist
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at", folder_path, "\n")
} else {
  cat("Folder already exists at", folder_path, "\n")
}

# Define the file name
file_name <- paste0("Power_results_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(power_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

