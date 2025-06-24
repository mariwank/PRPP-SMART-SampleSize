# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(R2jags, coda, tidyverse, utils) 


# Load in functions
source("DataGeneration_BJSM.R") # data generation function
source("LogOR_Function.R") # compute Log-OR differences for indifference DTRs

###################### Input parameter #######################

scenario <- 1
subscenario <- "LinkNot1" # Link1 or LinkNot1

if (scenario == 1){
  pNP1=0.50 #  desired proportion of individuals expressing No Preference in stage 1
  pNP2=0.50 # desired proportion of patients expressing No Preference in stage 2 (among non-responders)
  
} else if (scenario == 2){
  pNP1=1/3
  pNP2=1/3
  
}  else if (scenario == 3){
  pNP1=0.5
  pNP2=1/3
  
}


if (subscenario == "LinkNot1"){
  
  source("S_LinkNot1.R") # treatment outcome rates
  
  
} else if (subscenario == "Link1") {
  
  source("S_Link1.R") # treatment outcome rates

}

# List of sample sizes to iterate over
sample_sizes <- seq(300, 1000, 100)

# Specify theta targets
pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)

threshold = 0.2 # in this scenario it produces set of best AAC00 and AAD00
alpha_type1 = 0.05 

#################################################################

### TRUE DTRS ###
expected_pref <- c()  # expected DTR response rates from our simulated data
expected_pref[1] <- pi_A * pA0A + (1 - pi_A) * pi_AC  #AAC00
expected_pref[2] <- pi_A * pA0A + (1 - pi_A) * pi_AD  #AAD00
expected_pref[3] <- pi_B * pB0B + (1 - pi_B) * pi_BC  #BBC00
expected_pref[4] <- pi_B * pB0B + (1 - pi_B) * pi_BD  #BBD00
expected_pref[5] <- pi_A * pA0A + (1 - pi_A) * pA0C1 #AAC01
expected_pref[6] <- pi_A * pA0A + (1 - pi_A) * pA0D1 #AAD01
expected_pref[7] <- pi_B * pB0B + (1 - pi_B) * pB0C1 #BBC01
expected_pref[8] <- pi_B * pB0B + (1 - pi_B) * pB0D1 #BBD01
expected_pref[9] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C0 #AAC10
expected_pref[10] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D0 #AAD10
expected_pref[11] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C0 #BBC10
expected_pref[12] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D0 #BBD10
expected_pref[13] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C1 #AAC11
expected_pref[14] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D1 #AAD11
expected_pref[15] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C1 #BBC11
expected_pref[16] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D1 #BBD11

# Create to check nominal coverage
true_piDTR_mat <- matrix(c(pi_A, pi_B, pi_AC, pi_AD, pi_BC, pi_BD, expected_pref, alphaP_link, beta1_link, pi_A1, pi_B1, pA0C1, pA0D1, pB0C1, pB0D1, pA1C0, pA1D0, pB1C0, pB1D0, pA1C1, pA1D1, pB1C1, pB1D1, pA0A, pB0B, pA1A, pB1B))
rownames(true_piDTR_mat) <- c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD", "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11", "alpha1A", "alpha1B", "alpha2C", "alpha2D", "Beta1A", "Beta1B", "PiA1", "PiB1", "PiAC01", "PiAD01", "PiBC01", "PiBD01", "PiAC10", "PiAD10", "PiBC10", "PiBD10", "PiAC11", "PiAD11", "PiBC11", "PiBD11", "ThetaA", "ThetaB", "ThetaA1", "ThetaB1")


Emp_Power <- c() # store empirical power per N
Emp_BestIncuded <- c() # store empirical best DTR included in set

start_time <- Sys.time()
for (j in 1:length(sample_sizes)){
  N_j <- sample_sizes[j]
  # find number of sims needed to get 500 given the above settings 
  source("NSimToGet500.R")
  iterations_needed <- count_iterations()
  print(iterations_needed)
  
  ######################### BJSM Set-Up #################################
  
  ##### Estimation #####
  NUM_PATHS <- 6 # number of treatment paths independent of preference in SMART (AA, BB, AC, AD, BC, BD) = max(treatment_stageII)
  MCMC_SAMPLE <- 5500
  BURN.IN <- 500
  n_MCMC_chain <- 1
  n_MCMC_thin <- 1
  saved_parms = c('pi','beta', 'alphaP','DTR', 'pi_pref') # parameters to monitor 
  
  # prior setting 4 for BJSM (more variance)
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
    
    alphaP[1] ~ dgamma(alphaP_prior_a,alphaP_prior_b);  # alpha1A 
    alphaP[2] ~ dgamma(alphaP_prior_a,alphaP_prior_b);  # alpha1B 
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
    pi_pref[17] <- beta[1]*alphaP[1]*pi[1] # pi_A1A
    pi_pref[18] <- beta[2]*alphaP[2]*pi[2] # pi_B1B
  }
  
  ######################### SIMULATION #############################
  n.sim <- iterations_needed
  num_skip <- rep(0, n.sim) # number simulations skipped
  
  # store results
  DTR_hat_bjsm <- c()            # store posterior means 
  DTR_hat_bjsm_var <- c()        # store posterior variances
  pi_hat_bjsm <- c()             # store posterior means
  pi_hat_bjsm_var <- c()         # store posterior variances
  pi_pref_hat_bjsm <- c()        # store posterior means
  pi_pref_hat_bjsm_var <- c()    # store posterior variances
  beta_hat_bjsm <- c()           # store posterior means
  beta_hat_bjsm_var <- c()       # store posterior variances
  alpha_hat_bjsm <- c()          # store posterior means
  alpha_hat_bjsm_var <- c()      # store posterior variances
  theta_hat_bjsm <- c()
  theta_hat_bjsm_var <- c()
  ci_hat <- c()
  CheckExcludedDTRExclude <- c() # store whether the inferior DTRs are excluded from set of best (empirical power)
  BestIncluded <- c()            # store whether the best DTR was included in set of best (i.e., upper limit CI >= 0)
  SetBestCorrect <- c()          # store whether the computed set of best was correct based on threshold 
  
  
  for (i in 1:n.sim) {
    set.seed(i+100000)
    data <- generate_data(N=N_j, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
    
    bjsm_df <- data[[2]] # BJSM formatted data
    
    # check to make sure simulated data has at least 3 subjects per treatment path
    trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 1)
    if (nrow(trialpath_df) < 20) {num_skip[i]=1 ; next} # check to make sure data in each pathway 
    
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
    
    jags.out <- jags(data=data_list, model.file=bjsm_model,parameters.to.save=saved_parms,n.chains=n_MCMC_chain,
                     n.thin=n_MCMC_thin, n.iter=MCMC_SAMPLE,n.burnin=BURN.IN, progress.bar = "none", quiet=TRUE)
    
    
    store <- print(jags.out)
    mcmc.out <- summary(coda::as.mcmc(jags.out))
    
    # pi estimates posterior mean
    pi_hat_bjsm <- rbind(pi_hat_bjsm, store$mean$pi) # 6 columns n.sim rows "PiA", "PiB", "PiC", "PiD"
    
    # pi estimates posterior variance
    pi_hat_bjsm_var <- rbind(pi_hat_bjsm_var, store$sd$pi^2) # 6 columns n.sim rows "PiA", "PiB", "PiC", "PiD"
    
    # theta estimates posterior mean 
    theta_hat_bjsm <- rbind(theta_hat_bjsm, store$mean$pi_pref[c(15,16,17,18)]) # 4 columns n.sim rows "thetaA", "thetaB", "thetaA1", "thetaB1"
     
    # theta estimates posterior variance
    theta_hat_bjsm_var <- rbind(theta_hat_bjsm_var, store$sd$pi_pref[c(15,16,17,18)]^2) # 4 columns n.sim rows "thetaA", "thetaB", "thetaA1", "thetaB1"
  
    # nuisance pi estimates posterior mean
    pi_pref_hat_bjsm <- rbind(pi_pref_hat_bjsm, store$mean$pi_pref[1:14]) # 14 columns n.sim rows 
    
    # nuisance pi estimates posterior variance
    pi_pref_hat_bjsm_var <- rbind(pi_pref_hat_bjsm_var, store$sd$pi_pref[1:14]^2) # 14 columns n.sim rows 
    
    # DTR estimates posterior mean
    DTR_hat_bjsm <- rbind(DTR_hat_bjsm, store$mean$DTR) # 16 columsn n.sim rows "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"
    
    # DTR estimates posterior variance
    DTR_hat_bjsm_var <- rbind(DTR_hat_bjsm_var, store$sd$DTR^2) # 16 columsn n.sim rows "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"
    
    # beta estimates posterior mean
    beta_hat_bjsm <- rbind(beta_hat_bjsm, store$mean$beta) # 2 columsn n.sim rows "Beta1A", "Beta1B"
    
    # beta estimates posterior variances
    beta_hat_bjsm_var <- rbind(beta_hat_bjsm_var, store$sd$beta^2) # 2 columsn n.sim rows "Beta1A", "Beta1B"
    
    # alpha estimates posterior mean
    alpha_hat_bjsm <- rbind(alpha_hat_bjsm, store$mean$alphaP) # 4 columsn n.sim rows "AlphaP1" (1st stage preference), "AlphaP2 (2nd stage preference)
    
    # alpha estimates posterior variances
    alpha_hat_bjsm_var <- rbind(alpha_hat_bjsm_var, store$sd$alphaP^2) # 4 columsn n.sim rows "AlphaP1" (1st stage preference), "AlphaP2 (2nd stage preference)
    
    # Coverage 
    quantile_mat <- mcmc.out$quantiles[c(1:6,8:47),c(1,5)] 
    # reorder quantile mat to be in the same order as true_piDTR_mat
    rowname_order <- c("pi[1]", "pi[2]", "pi[3]", "pi[4]", "pi[5]", "pi[6]", "DTR[1]", "DTR[2]", "DTR[3]", "DTR[4]", "DTR[5]", "DTR[6]", "DTR[7]", "DTR[8]", "DTR[9]", "DTR[10]", "DTR[11]", "DTR[12]", "DTR[13]", "DTR[14]", "DTR[15]", "DTR[16]", "alphaP[1]", "alphaP[2]", "alphaP[3]", "alphaP[4]", "beta[1]", "beta[2]", "pi_pref[1]", "pi_pref[2]", "pi_pref[3]", "pi_pref[4]", "pi_pref[5]", "pi_pref[6]", "pi_pref[7]", "pi_pref[8]", "pi_pref[9]", "pi_pref[10]", "pi_pref[11]", "pi_pref[12]", "pi_pref[13]", "pi_pref[14]", "pi_pref[15]", "pi_pref[16]", "pi_pref[17]", "pi_pref[18]")
    quantile_mat <- quantile_mat[rowname_order,,drop=FALSE]
    param_in_ci <- data.table::between(true_piDTR_mat[,1], quantile_mat[, 1], quantile_mat[, 2])
    ci_hat <- rbind(ci_hat, param_in_ci) 
    
    ### MCB with indifference DTRs #####

      
    logORThreshold <- LogOR(response_prob = c(pA0A,pB0B,pi_AC,pi_AD,pi_BC,pi_BD),
                              stage_one_trt_one_response_prob = pi_A,
                              stage_one_trt_two_response_prob = pi_B)
      

    
    thetadraws <- store$sims.matrix[,1:4] # estimated indifference DTRs
    
    #Compute log-OR
    thetadraws_log_odds <- log(thetadraws/(1-thetadraws)) # estimated log-OR indifference DTRs
    
    #Compute index of best indifference DTR
    max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
    
    #Compute log-odds ratios between each indifference DTR and best
    Log_OR_matrix <- thetadraws_log_odds-matrix(thetadraws_log_odds[,max_odds_ind],nrow=nrow(thetadraws),ncol=4)
    
    #Rank log-OR
    rank_matrix <- apply(Log_OR_matrix[,-max_odds_ind],2,rank,ties.method = 'random')
    
    #Find max rank
    rank_max <- apply(rank_matrix,1,max)
    
    #Create sorted log-OR
    new_dat <- apply(Log_OR_matrix[,],2,sort)
    
    #Compute 100(1-alpha)% upper quantile
    ranks_quantile <- ceiling(quantile(rank_max,1-alpha_type1))
    
    # Compute upper limit of credible interval. One for each log-OR which determines the set of best.
    upper_limit <- new_dat[ranks_quantile,] 
    
    # indicies 
    rejection_indices <- which(abs(logORThreshold)>=threshold) # DTRs that should be excluded from set of best based on threshold 
    best_index <- which(abs(logORThreshold)==0)

    CheckExcludedDTRExclude <- rbind(CheckExcludedDTRExclude, all(upper_limit[rejection_indices] < 0)) # all DTRs that should be excluded should have upper limit CI < 0
    
    BestIncluded <- rbind(BestIncluded, upper_limit[best_index] >= 0)
    
  }
  
  ## Empirical Power per N
  Emp_Power[j] <- mean(CheckExcludedDTRExclude, na.rm=TRUE)
  Emp_BestIncuded[j] <- mean(BestIncluded, na.rm=TRUE)

  
  ######################### EVALUATION ############################
  
  ### Save Settings
  # set date and time for file saving 
  st<-format(Sys.time(), "%Y_%m_%d_%H_%M")
  
  # Define the folder path where you want your results saved
  # folder_path <- "path to your folder where you want to store results"
  
  # Example: 
  folder_path <- paste0("/home/mariwank/Final_SampleSize_Results/Empirical/SimResults/Scenario", scenario, "/", subscenario, "/N", sample_sizes[j])
  
  # Create the folder if it doesn't exist
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created at", folder_path, "\n")
  } else {
    cat("Folder already exists at", folder_path, "\n")
  }
  
  # Power per N
  Power_perN <- data.frame(N = sample_sizes[j], 
                      EmpiricalPower = mean(CheckExcludedDTRExclude, na.rm=TRUE),
                      BestIncuded = mean(BestIncluded, na.rm=TRUE),
                      SetBestCorrect = mean(SetBestCorrect, na.rm=TRUE))
  
  
  #Define the file name
  file_name <- paste0("EmpiricalPower_Results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(Power_perN, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## Linkage Parameters (alpha/beta)
  beta_alpha_effect_output <- data.frame(link_param = c("Beta1A", "Beta1B", "Alpha1A", "Alpha1B", "Alpha2C", "Alpha2D"),
                                         true_param = c(beta1_link[1], beta1_link[2], alphaP_link[1], alphaP_link[2], alphaP_link[3], alphaP_link[4]),
                                         param_hat = c(apply(beta_hat_bjsm,2,mean, na.rm = TRUE), apply(alpha_hat_bjsm,2,mean, na.rm = TRUE)),
                                         sd_param_hat = c(apply(beta_hat_bjsm,2,sd, na.rm = TRUE), apply(alpha_hat_bjsm,2,sd, na.rm = TRUE)),
                                         avg_sd_param_hat = c(sqrt(apply(beta_hat_bjsm_var,2,mean, na.rm = TRUE)), sqrt(apply(alpha_hat_bjsm_var,2,mean, na.rm = TRUE))))
  beta_alpha_effect_output$bias <- beta_alpha_effect_output$param_hat - beta_alpha_effect_output$true_param
  beta_alpha_effect_output$rMSE <- sqrt(beta_alpha_effect_output$sd_param_hat^2 + beta_alpha_effect_output$bias^2)
  
  # Define the file name
  file_name <- paste0("LinkageParam_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(beta_alpha_effect_output, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## No Preference Main Effects
  trt_effect_output <- data.frame(pi = c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD"),
                                  true_pi = c(pi_A,pi_B,pi_AC,pi_AD,pi_BC,pi_BD),
                                  pi_hat = apply(pi_hat_bjsm,2,mean, na.rm = TRUE),
                                  sd_pi_hat = apply(pi_hat_bjsm,2,sd, na.rm = TRUE),
                                  avg_sd_pi_hat = sqrt(apply(pi_hat_bjsm_var,2,mean, na.rm = TRUE)))
  trt_effect_output$bias <- trt_effect_output$pi_hat - trt_effect_output$true_pi
  trt_effect_output$rMSE <- sqrt(trt_effect_output$sd_pi_hat^2 + trt_effect_output$bias^2)
  
  # Define the file name
  file_name <- paste0("MainEffectNoPref_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(trt_effect_output, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## ThetaJ
  theta_effect_output <- data.frame(pi = c("ThetaA", "ThetaB", "ThetaA1", "ThetaB1"),
                                true_theta = c(pA0A,pB0B,pA1A,pB1B),
                                theta_hat = apply(theta_hat_bjsm,2,mean, na.rm = TRUE),
                                sd_theta_hat = apply(theta_hat_bjsm,2,sd, na.rm = TRUE),
                                avg_sd_theta_hat = sqrt(apply(theta_hat_bjsm_var,2,mean, na.rm = TRUE)))

  theta_effect_output$bias <- theta_effect_output$theta_hat - theta_effect_output$true_theta
  theta_effect_output$rMSE <- sqrt(theta_effect_output$sd_theta_hat^2 + theta_effect_output$bias^2)
  
  # Define the file name
  file_name <- paste0("theta_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(theta_effect_output, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")

  
  ## Nuisance Main Effects
  trt_effect_output_pref <- data.frame(pi = c("PiA1", "PiB1", "PiAC01", "PiAD01", "PiBC01", "PiBD01", "PiAC10", "PiAD10", "PiBC10", "PiBD10", "PiAC11", "PiAD11", "PiBC11", "PiBD11"),
                                       true_pi = c(pi_A1, pi_B1, pA0C1, pA0D1, pB0C1, pB0D1, pA1C0, pA1D0, pB1C0, pB1D0, pA1C1, pA1D1, pB1C1, pB1D1),
                                       pi_hat = apply(pi_pref_hat_bjsm,2,mean, na.rm = TRUE),
                                       sd_pi_hat = apply(pi_pref_hat_bjsm,2,sd, na.rm = TRUE),
                                       avg_sd_pi_hat = sqrt(apply(pi_pref_hat_bjsm_var,2,mean, na.rm = TRUE)))
  trt_effect_output_pref$bias <- trt_effect_output_pref$pi_hat - trt_effect_output_pref$true_pi
  trt_effect_output_pref$rMSE <- sqrt(trt_effect_output_pref$sd_pi_hat^2 + trt_effect_output_pref$bias^2)
  
  # Define the file name
  file_name <- paste0("Nuisance_MainEffect_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(trt_effect_output_pref, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## DTR
  dtr_effect_output <- data.frame(dtr = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"),
                                  true_dtr = expected_pref,
                                  DTR_hat = apply(DTR_hat_bjsm,2, mean, na.rm = TRUE),
                                  sd_DTR_hat = apply(DTR_hat_bjsm,2,sd, na.rm = TRUE),
                                  avg_sd_DTR_hat = sqrt(apply(DTR_hat_bjsm_var,2, mean, na.rm = TRUE)))
  dtr_effect_output$bias <- dtr_effect_output$DTR_hat - dtr_effect_output$true_dtr
  dtr_effect_output$rMSE <- sqrt(dtr_effect_output$sd_DTR_hat^2 + dtr_effect_output$bias^2)
  
  # Define the file name
  file_name <- paste0("DTR_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(dtr_effect_output, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
  ## coverage 
  params <- c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD", "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11", "alphaP1A",  "alphaP1B", "alphaP2C", "alphaP2D", "Beta1A", "Beta1B", "PiA1", "PiB1", "PiAC01", "PiAD01", "PiBC01", "PiBD01", "PiAC10", "PiAD10", "PiBC10", "PiBD10", "PiAC11", "PiAD11", "PiBC11", "PiBD11", "PiA0A", "PiB0B", "PiA1A", "PiB1B")
  coverage_output <- data.frame(Parameter = params, 
                                Coverage = apply(ci_hat, 2, mean, na.rm = TRUE))
  
  # Define the file name
  file_name <- paste0("CI_results_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(coverage_output, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## Sims skipped
  num_skip_total <- sum(num_skip)
  mean_skip <- mean(num_skip)
  skip_df <- data.frame(num_sim = n.sim, num_skip_total = num_skip_total, avg_num_skip = mean_skip)
  
  # Define the file name
  file_name <- paste0("Num_Sim_Skip_Summary_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(skip_df, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
}
end_time <- Sys.time()

############## Overall Power Results ###################
### Save Settings
# set date and time for file saving 
st<-format(Sys.time(), "%Y_%m_%d_%H_%M")

# Define the folder path where you want your results saved
# folder_path <- "path to your folder where you want to store results"

# Example: 
folder_path <- paste0("/home/mariwank/Final_SampleSize_Results/Empirical/SimResults/Scenario", scenario, "/", subscenario)

# Create the folder if it doesn't exist
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at", folder_path, "\n")
} else {
  cat("Folder already exists at", folder_path, "\n")
}


Power_AllN <- data.frame(N = sample_sizes, 
                    EmpiricalPower = Emp_Power,
                    BestIncuded = Emp_BestIncuded,
                    ComputationTime = as.numeric(end_time - start_time))


#Define the file name
file_name <- paste0("EmpiricalPower_AllN_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(Power_AllN, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")
