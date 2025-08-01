# Code adapted from Artman et al 2022. 


library(R2jags)
library(coda)

source("PRPPSMART_DataGen.R") # data generation function
source("LogOR_Function.R")

I <- 100 # number of monte carlo integrations to perform for the data
scenario <- 1
subscenario <- "Link1" # Link1 or LinkNot1
seed_spec <- 1 # seed for reproducibility 

# Treatment Preference Proportions
if (scenario == 1){
  pP1=0.50 #  desired proportion of individuals expressing Preference in stage 1
  pP2=0.50 # desired proportion of patients expressing Preference in stage 2 (among non-responders)
  
} else if (scenario == 2){
  pP1=2/3
  pP2=2/3
  
}  else if (scenario == 3){
  pP1=0.5
  pP2=2/3
  
}

# Specify treatment preferences
PreferA_prob=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
PreferC_prob=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)

pNP1 <- 1 - pP1 #  desired proportion of individuals expressing no Preference in stage 1
pNP2 <- 1 - pP2 # desired proportion of patients expressing no Preference in stage 2 (among non-responders)



# linkage parameters
if (subscenario == "LinkNot1"){
  
  linkage_parameter_values <- c(1, 1, 1.1, 1.02, 1.05, 1.2) 
  
} else if (subscenario == "Link1"){
  
  linkage_parameter_values <- c(1, 1, 1, 1, 1, 1) 
  
}

# indifference response rates
pi_A <- 0.6           # stage 1 outcome rate to randomize A: Pr(Y1=1|T1=A,P1=0)                                 
pi_B <- 0.45          # stage 1 outcome rate to randomize B: Pr(Y1=1|T1=B,P1=0)     
pi_AC <- 0.5          # stage 2 outcome rate of non-responders to randomized A and randomize C in the second stage: Pr(Y2=1|T1=A,P1=0,Y1=0,T2=C,P2=0)                                
pi_AD <- 0.4          # stage 2 outcome rate of non-responders to randomized A and randomize D in the second stage: Pr(Y2=1|T1=A,P1=0,Y1=0,T2=D,P2=0)                                
pi_BC <- 0.3          # stage 2 outcome rate of non-responders to randomized B and randomize C in the second stage: Pr(Y2=1|T1=B,P1=0,Y1=0,T2=C,P2=0)                               
pi_BD <- 0.2          # stage 2 outcome rate of non-responders to randomized B and randomize D in the second stage: Pr(Y2=1|T1=B,P1=0,Y1=0,T2=D,P2=0)


# List of sample sizes to iterate over
sample_sizes <- seq(300, 1000, 100)
threshold = 0.2 # threshold of exclusion from set of best (minimum detectable difference)
alpha_type1 = 0.05 # probability of excluding best DTR 



# linkage parameters
beta1A <- linkage_parameter_values[1] # responders to A second stage response 
beta1B <- linkage_parameter_values[2] # responders to B second stage response
alpha1A <- linkage_parameter_values[3] 
alpha1B <- linkage_parameter_values[4]
alpha2C <- linkage_parameter_values[5]
alpha2D <- linkage_parameter_values[6]


# First stage preference treatment response rates
pi_A1 <- pi_A*alpha1A                                # stage 1 outcome rate to prefer A: Pr(Y1=1|T1=A,P1=1)
pi_B1 <- pi_B*alpha1B                                # stage 1 outcome rate to prefer B: Pr(Y1=1|T1=B,P1=1)
# Second stage preference treatment response rates
pA0A <- pi_A * beta1A                                # Second stage outcome rate of responders to randomized A: Pr(Y2=1|T1=A,P1=0,Y1=1)
pB0B <- pi_B * beta1B                                # Second stage outcome rate of responders to randomized B: Pr(Y2=1|T1=B,P1=0,Y1=1)
pA1A <- pi_A * alpha1A * beta1A                      # Second stage outcome rate of responders to preferred A: Pr(Y2=1|T1=A,P1=1,Y1=1)
pB1B <- pi_B * alpha1B * beta1B                      # Second stage response rate of responders to preferred B: Pr(Y2=1|T1=B,P1=1,Y1=1)
pA0C1 <- alpha2C * pi_AC                             # stage 2 outcome rate of non-responders to randomized A and prefer C in the second stage: Pr(Y2=1|T1=A,P1=0,Y1=0,T2=C,P2=1)
pA0D1 <- alpha2D * pi_AD                             # stage 2 outcome rate of non-responders to randomized A and prefer D in the second stage: Pr(Y2=1|T1=A,P1=0,Y1=0,T2=D,P2=1)
pA1C0 <- alpha1A * pi_AC                             # stage 2 outcome rate of non-responders to prefer A and randomize C in the second stage: Pr(Y2=1|T1=A,P1=1,Y1=0,T2=C,P2=0)
pA1D0 <- alpha1A * pi_AD                             # stage 2 outcome rate of non-responders to prefer A and randomize D in the second stage: Pr(Y2=1|T1=A,P1=1,Y1=0,T2=D,P2=0)
pA1C1 <- alpha1A * alpha2C * pi_AC                   # stage 2 outcome rate of non-responders to prefer A and prefer C in the second stage: Pr(Y2=1|T1=A,P1=1,Y1=0,T2=C,P2=1)
pA1D1 <- alpha1A * alpha2D * pi_AD                   # stage 2 outcome rate of non-responders to prefer A and prefer C in the second stage: Pr(Y2=1|T1=A,P1=1,Y1=0,T2=D,P2=1)
pB0C1 <- alpha2C * pi_BC                             # stage 2 outcome rate of non-responders to randomize B and prefer C in the second stage: Pr(Y2=1|T1=B,P1=0,Y1=0,T2=C,P2=1)
pB0D1 <- alpha2D * pi_BD                             # stage 2 outcome rate of non-responders to randomize B and prefer D in the second stage: Pr(Y2=1|T1=B,P1=0,Y1=0,T2=D,P2=1)
pB1C0 <- alpha1B * pi_BC                             # stage 2 outcome rate of non-responders to prefer B and randomize C in the second stage: Pr(Y2=1|T1=B,P1=1,Y1=0,T2=C,P2=0)
pB1D0 <- alpha1B * pi_BD                             # stage 2 outcome rate of non-responders to prefer B and randomize D in the second stage: Pr(Y2=1|T1=B,P1=1,Y1=0,T2=D,P2=0)
pB1C1 <- alpha1B * alpha2C * pi_BC                   # stage 2 outcome rate of non-responders to prefer B and prefer C in the second stage: Pr(Y2=1|T1=B,P1=1,Y1=0,T2=C,P2=1)
pB1D1 <- alpha1B * alpha2D * pi_BD                   #  stage 2 outcome rate of non-responders to prefer B and prefer D in the second stage: Pr(Y2=1|T1=B,P1=1,Y1=0,T2=D,P2=1)




NUM_PATHS <- 6 # number of treatment paths independent of preference in SMART (AA, BB, AC, AD, BC, BD) - max(treatment_stageII)
MCMC_SAMPLE <- 10500
BURN.IN <- 500
n_MCMC_chain <- 1
n_MCMC_thin <- 1
saved_parms = c('pi','beta', 'alphaP','DTR', 'pi_pref') # parameters to monitor 
pi_prior.a <- c(rep(1/3,6))  # mean of pi_j, pi_jk = a / a + b = (1/3) / (1/3+1/3) = 0.5, pi ~ beta(a,b)
pi_prior.b <- c(rep(1/3,6)) 
beta1_prior.r <- 1/2           # mean of beta1 = r/mu = 1/1 = 1, beta1 ~ gamma(r,mu) gives a gamma prior with more variance (smaller values more variance)
beta1_prior.mu <- 1/2     
alphaP_prior.r <- 1/2         # mean of alphaP = r/mu = 1/1 = 1, alphaP ~ gamma(r,mu) gives a gamma prior with more variance 
alphaP_prior.mu <- 1/2

# Jags model
bjsm_model = function()
{ 
  for (i in 1:N){   # n is total sample size
    # likelihood
    Y1[i]~dbern(pi_1[i])
    Y2[i]~dbern(pi_2[i])
    # explaining
    pi_1[i] <- ifelse(preference_stageI[i] == 0, pi[treatment_stageI[i]], # A0/B0
                      ifelse(preference_stageI[i] == 1 && treatment_stageI[i] == 1, alphaP[1] * pi[treatment_stageI[i]], alphaP[2] * pi[treatment_stageI[i]])) # A1 and B1
    
    
    pi_2[i] <- ifelse(Y1[i] == 1 && preference_stageI[i] == 0, 
                      pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # A0A/B0B
                      ifelse(Y1[i] == 1 && preference_stageI[i] == 1 && treatment_stageI[i] == 1,  alphaP[1] * pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # A1A
                             ifelse(Y1[i] == 1 && preference_stageI[i] == 1 && treatment_stageI[i] == 2, alphaP[2] * pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # B1B
                                    ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 0, pi[treatment_stageII[i]], # A0C0/A0D0/B0C0/B0D0
                                           ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 1 && treatmentCD_stageII[i] == 3, alphaP[3] * pi[treatment_stageII[i]], # A0C1/B0C1
                                                  ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 1 && treatmentCD_stageII[i] == 4, alphaP[4] * pi[treatment_stageII[i]], # A0D1/B0D1
                                                         ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 0 && treatment_stageI[i] == 1, alphaP[1] * pi[treatment_stageII[i]], # A1C0/A1D0
                                                                ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 0 && treatment_stageI[i] == 2, alphaP[2] * pi[treatment_stageII[i]], # B1C0/B1D0 
                                                                       ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 1 && treatmentCD_stageII[i] == 3, alphaP[1] * alphaP[3] * pi[treatment_stageII[i]], # A1C1
                                                                              ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 1 && treatmentCD_stageII[i] == 4, alphaP[1] * alphaP[4] * pi[treatment_stageII[i]], # A1D1
                                                                                     ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 2 && treatmentCD_stageII[i] == 3, alphaP[2] * alphaP[3] * pi[treatment_stageII[i]], alphaP[2] * alphaP[4] * pi[treatment_stageII[i]]))))))))))) # B1C1 and B1D1
    
    
    
    
  }
  
  alphaP[1] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha1A
  alphaP[2] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha1B
  alphaP[3] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha2C
  alphaP[4] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha2D
  beta[1] ~ dgamma(beta1_prior_a,beta1_prior_b)     # beta1A
  beta[2] ~ dgamma(beta1_prior_a,beta1_prior_b)     # beta1B
  
  for (j in 1:num_paths){
    pi[j] ~ dbeta(pi_prior_a[j],pi_prior_b[j]); T(0.001,0.999)
  }
  
  # DTR
  DTR[1] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * pi[3] ## AAC00
  DTR[2] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * pi[4] ## AAD00
  DTR[3] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * pi[5] ## BBC00
  DTR[4] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * pi[6] ## BBD00
  DTR[5] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * (alphaP[3]*pi[3]) ## AAC01
  DTR[6] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * (alphaP[4]*pi[4]) ## AAD01
  DTR[7] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * (alphaP[3]*pi[5]) ## BBC01
  DTR[8] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * (alphaP[4]*pi[6])  ## BBD01
  DTR[9] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1]) + (1 - alphaP[1]*pi[1]) * (alphaP[1]*pi[3]) ## AAC10
  DTR[10] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1]) + (1 - alphaP[1]*pi[1]) * (alphaP[1]*pi[4]) ## AAD10
  DTR[11] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2]) + (1 - alphaP[2]*pi[2]) * (alphaP[2]*pi[5]) ## BBC10
  DTR[12] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2]) + (1 - alphaP[2]*pi[2]) * (alphaP[2]*pi[6]) ## BBD10
  DTR[13] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1])  + (1 - alphaP[1]*pi[1]) * (alphaP[1]*alphaP[3]*pi[3]) ## AAC11
  DTR[14] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1])  + (1 - alphaP[1]*pi[1]) * (alphaP[1]*alphaP[4]*pi[4]) ## AAD11
  DTR[15] =  (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2])  + (1 - alphaP[2]*pi[2]) * (alphaP[2]*alphaP[3]*pi[5]) ## BBC11
  DTR[16] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2])  + (1 - alphaP[2]*pi[2]) * (alphaP[2]*alphaP[4]*pi[6]) ## BBD11
  
  # Nuisance main effects
  pi_pref[1] <- alphaP[1]*pi[1] # pi_A1
  pi_pref[2] <- alphaP[2]*pi[2] # pi_B1
  pi_pref[3] <- alphaP[3]*pi[3] # pi_AC01
  pi_pref[4] <- alphaP[4]*pi[4] # pi_AD01
  pi_pref[5] <- alphaP[3]*pi[5] # pi_BC01
  pi_pref[6] <- alphaP[4]*pi[6] # pi_BD01
  pi_pref[7] <- alphaP[1]*pi[3] # pi_AC10
  pi_pref[8] <- alphaP[1]*pi[4] # pi_AD10
  pi_pref[9] <- alphaP[2]*pi[5] # pi_BC10
  pi_pref[10] <- alphaP[2]*pi[6] # pi_BD10
  pi_pref[11] <- alphaP[1]*alphaP[3]*pi[3] # pi_AC11
  pi_pref[12] <- alphaP[1]*alphaP[4]*pi[4] # pi_AD11
  pi_pref[13] <- alphaP[2]*alphaP[3]*pi[5] # pi_BC11
  pi_pref[14] <- alphaP[2]*alphaP[4]*pi[6] # pi_BD11
  pi_pref[15] <- beta[1]*pi[1] # pi_A0A
  pi_pref[16] <- beta[2]*pi[2] # pi_B0B
  
}


upper_limit <- array(NA, dim = c(I, 10, 4)) # where to store CI upper limit for indifference DTRs
power <- rep(NA, length(sample_sizes)) # store power per N
bestinset <- rep(NA, length(sample_sizes))


# Set the desired number of simulations per N
total_simulations <- I

iteration_number_store <- rep(NA, length(sample_sizes))
successful_sim_store <- rep(NA, length(sample_sizes))

start_time <- Sys.time()
for (n in 1:length(sample_sizes)){
  
  # Initialize the number of successful simulations
  successful_simulations <- 0
  iteration_number <- 0
  
  N <- sample_sizes[n]
  
  while (successful_simulations < total_simulations){
    
    
    iteration_number <- iteration_number + 1
    
    set.seed(iteration_number + seed_spec)
    data <- generate_data(N=N, pNP1=pNP1, pTheta_A=PreferA_prob, pNP2=pNP2, pTheta_C=PreferC_prob, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
    
    bjsm_df <- data[[2]] # BJSM formatted data
    
    # check to make sure simulated data has at least 1 subjects per treatment path
    trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 1)
    if (nrow(trialpath_df) < 20) {next} # check to make sure data in each pathway 
    
    data_list <- list(N = nrow(bjsm_df),
                      num_paths = NUM_PATHS, # max(treatment_stageII)
                      Y1 = bjsm_df$response_stageI,
                      Y2 = bjsm_df$response_stageII,
                      treatment_stageI = bjsm_df$treatment_stageI,
                      treatment_stageII = bjsm_df$treatment_stageII,
                      preference_stageI = bjsm_df$preference_stageI,
                      preference_stageII = bjsm_df$preference_stageII,
                      treatmentCD_stageII = bjsm_df$treatmentCD_stageII,
                      pi_prior_a = pi_prior.a, # beta
                      pi_prior_b = pi_prior.b, # beta 
                      alphaP_prior_a = alphaP_prior.r, # gamma
                      alphaP_prior_b = alphaP_prior.mu, # gamma
                      beta1_prior_a = beta1_prior.r,   # gamma
                      beta1_prior_b = beta1_prior.mu  # gamma
    )
    
    # Increment the successful simulation counter
    successful_simulations <- successful_simulations + 1
    
    jags.out <- jags(data=data_list, model.file=bjsm_model,parameters.to.save=saved_parms,n.chains=n_MCMC_chain,
                     n.thin=n_MCMC_thin, n.iter=MCMC_SAMPLE,n.burnin=BURN.IN, progress.bar = "none", quiet=TRUE)
    
    
    store <- print(jags.out)
    posterior_piA <- store$sims.list$pi[,1]
    posterior_piB <- store$sims.list$pi[,2]
    posterior_piAC <- store$sims.list$pi[,3]
    posterior_piAD <- store$sims.list$pi[,4]
    posterior_piBC <- store$sims.list$pi[,5]
    posterior_piBD <- store$sims.list$pi[,6]
    posterior_piA0A <- store$sims.list$pi_pref[,15]
    posterior_piB0B <- store$sims.list$pi_pref[,16]
    
    
    group_ids <- rep(1:10, each = 1000)
    
    piA_post_bjsm_sample <- split(posterior_piA, group_ids) # Split draws into J=10 samples of 1000 
    piB_post_bjsm_sample <- split(posterior_piB, group_ids) # Split draws into J=10 samples of 1000 
    piAC_post_bjsm_sample <- split(posterior_piAC, group_ids) # Split draws into J=10 samples of 1000 
    piAD_post_bjsm_sample <- split(posterior_piAD, group_ids) # Split draws into J=10 samples of 1000 
    piBC_post_bjsm_sample <- split(posterior_piBC, group_ids) # Split draws into J=10 samples of 1000 
    piBD_post_bjsm_sample <- split(posterior_piBD, group_ids) # Split draws into J=10 samples of 1000 
    piA0A_post_bjsm_sample <- split(posterior_piA0A, group_ids) # Split draws into J=10 samples of 1000 
    piB0B_post_bjsm_sample <- split(posterior_piB0B, group_ids) # Split draws into J=10 samples of 1000 
    
    
    
    for (j in 1:10){ # number of parameter draws for Monte Carlo integration for the parameter integral 
      
      # Draw 1000 draws from the posterior of the probability of outcome at the end of the study for each of the 6 indifference embedded treatment sequences.
      # A0A (thetaA), B0B (thetaB), AC, AD, BC, BD
      thetaA_draws <-  piA0A_post_bjsm_sample[[j]] # A0A
      thetaB_draws <-  piB0B_post_bjsm_sample[[j]] # B0B
      piAC_draws <- piAC_post_bjsm_sample[[j]] # A0C0
      piAD_draws <- piAD_post_bjsm_sample[[j]] # A0D0
      piBC_draws <- piBC_post_bjsm_sample[[j]] # B0C0
      piBD_draws <- piBD_post_bjsm_sample[[j]] # B0D0
      
      # Draw 1000 draws from the posterior of the probability of outcome at end of stage 1 for each of the indifference stage 1 treatment sequences 
      A0_draws <- piA_post_bjsm_sample[[j]]
      B0_draws <- piB_post_bjsm_sample[[j]]
      
      
      # Compute embedded DTR response probabilities using Robins' G-computation method
      dtrdraws <- cbind(thetaA_draws * A0_draws + piAC_draws * (1 - A0_draws), # AAC00
                          thetaA_draws * A0_draws + piAD_draws * (1 - A0_draws), # AAD00
                          thetaB_draws * B0_draws + piBC_draws * (1 - B0_draws), # BBC00
                          thetaB_draws * B0_draws + piBD_draws * (1 - B0_draws)) # BBD00
      
      
          
      logORThreshold <- LogOR(response_prob = c(pA0A,pB0B,pi_AC,pi_AD,pi_BC,pi_BD),
                                stage_one_trt_one_response_prob = pi_A,
                                stage_one_trt_two_response_prob = pi_B)

        
      dtrdraws_log_odds <- log(dtrdraws / (1 - dtrdraws))
      max_odds_ind <- which.max(colMeans(dtrdraws_log_odds))
      
      log_OR_matrix <- dtrdraws_log_odds - matrix(dtrdraws_log_odds[, max_odds_ind], nrow = 1000, ncol = 4)
      rank_matrix <- apply(log_OR_matrix[, -max_odds_ind], 2, rank, ties.method = 'random')
      rank_max <- apply(rank_matrix, 1, max)
      
      new_dat <- apply(log_OR_matrix[,], 2, sort)
      ranks_quantile <- floor(quantile(rank_max, probs = 1 - alpha_type1))
      
      upper_limit[successful_simulations,j,] <- new_dat[ranks_quantile,]
      
      
      
      
    }
    
  }
  
  rejection_indices <- which(abs(logORThreshold)>=threshold)
  
  
  if (length(rejection_indices)==1) {
    power[n] <- mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0)
    
  } else {
    power[n] <- mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod))
  }
  
  
  # best index
  best_index <- which(abs(logORThreshold)==0)
  bestinset[n] <- mean(apply(upper_limit, 3, function(x) x)[, best_index] >= 0)
  
  iteration_number_store[n] <- iteration_number
  successful_sim_store[n] <- successful_simulations

  
}

end_time <- Sys.time()


power_df <- data.frame(
  N = sample_sizes,  # Sample sizes
  Predicted_Power = power  # Extract the power from each list
)


power_df$BestIncluded <- bestinset
power_df$ComputationTime <- end_time - start_time


### Save Settings
# set date and time for file saving 
st<-format(Sys.time(), "%Y_%m_%d_%H_%M")

# Define the folder path where you want your results saved
# folder_path <- "path to your folder where you want to store results"

# Example: 
folder_path <- paste0("/home/mariwank/Final_SampleSize_Results/SimResults/Scenario", scenario, "/", subscenario, "/I", I, "/seed", seed_spec)

# Create the folder if it doesn't exist
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at", folder_path, "\n")
} else {
  cat("Folder already exists at", folder_path, "\n")
}

### Power ###
# Define the file name
file_name <- paste0("PredictedPower_results_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(power_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")


### Positivity ###

pos_df <- data.frame(N = sample_sizes, success_sim = successful_sim_store, total_sim = iteration_number_store)

# Define the file name
file_name <- paste0("Positivity_results_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(pos_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")



