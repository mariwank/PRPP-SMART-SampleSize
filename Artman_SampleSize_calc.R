# Artman traditional SMART sample size 

library(SMARTbayesR)
source("LogOR_Function.R")

# Corresponding PRPP-SMART input parameters for traditional SMART
beta1_link <- c(1,1)    # Beta1A, Beta1B
pi_A <- 0.6                                          
pi_B <- 0.45  
pi_AC <- 0.5                                          
pi_AD <- 0.4                                         
pi_BC <- 0.3                                         
pi_BD <- 0.2   
pA0A <- pi_A * beta1_link[1]                        
pB0B <- pi_B * beta1_link[2] 

delta <- 0.2 # threshold of exclusion 
logORThreshold <- LogOR(response_prob = c(pA0A,pB0B,pi_AC,pi_AD,pi_BC,pi_BD),
                        stage_one_trt_one_response_prob = pi_A,
                        stage_one_trt_two_response_prob = pi_B)

# n_min for PRPP-SMART sample size calculation
n_min <- PowerBayesian(design = "design-1",
                   sample_size = 300,
                   response_prob = c(pA0A, pi_AC, pi_AD, pB0B, pi_BC, pi_BD),
                   stage_one_trt_one_response_prob = pi_A,
                   stage_one_trt_two_response_prob = pi_B,
                   type = "log-OR",
                   threshold = 0.2,
                   alpha = 0.05) 

