###########################################################################################################################

# Data Simulation code for PRPP-SMART

# Author: Mari Wank
###########################################################################################################################
# Description: 
#    Code simulates a dataset from a two-stage PRPP-SMART design with a binary end of stage outcomes. 

###########################################################################################################################

#Function: gendata


#Purpose: This function generates one simulation dataset for a two stage PRPP-SMART trial with a binary outcome. 



# Required Parameters: 
#         method: Method you will use to analyze the data (BJSM or WRRM). This only impacts the output format of the simulated data. 
#         N: Number of individuals in the generated dataset
#         pNP1: Desired proportion of individuals expressing No Preference in stage 1
#         pTheta_A: Desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
#         pNP2: Desired proportion of patients expressing No Preference in stage 2 (among non-responders)
#         pTheta_C: Desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)
#         pi_A:  Probability of responding in stage 1 given randomized treatment A
#         pi_B:  Probability of responding in stage 1 given randomized treatment B
#         pi_A1: Probability of responding in stage 1 given preferred treatment A
#         pi_B1: Probability of responding in stage 1 given preferred treatment B
#         pi_AC: Probability of Y=1 given treatment path A0C0
#         pi_AD: Probability of Y=1 given treatment path A0D0
#         pi_BC: Probability of Y=1 given treatment path B0C0
#         pi_BD: Probability of Y=1 given treatment path B0D0
#         pA0A: Probability of Y=1 given treatment path A0A
#         pB0B: Probability of Y=1 given treatment path B0B
#         pA1A: Probability of Y=1 given treatment path A1A
#         pB1B: Probability of Y=1 given treatment path B1B
#         pA1C1: Probability of Y=1 given treatment path A1C1
#         pA1C0: Probability of Y=1 given treatment path A1C0
#         pA1D1: Probability of Y=1 given treatment path A1D1
#         pA1D0: Probability of Y=1 given treatment path A1D0
#         pA0C1: Probability of Y=1 given treatment path A0C1
#         pA0D1: Probability of Y=1 given treatment path A0D1
#         pB1C1: Probability of Y=1 given treatment path B1C1
#         pB1C0: Probability of Y=1 given treatment path B1C0
#         pB1D1: Probability of Y=1 given treatment path B1D1
#         pB1D0: Probability of Y=1 given treatment path B1D0
#         pB0C1: Probability of Y=1 given treatment path B0C1
#         pB0D1: Probability of Y=1 given treatment path B0D1

#Output:   Method: WRRM
#       (1) A two-stage PRPP-SMART dataset of N subjects 
#             The dataset has the following variables:
#               ID: Numeric subject ID variable 
#               X1: A continuous baseline variable generated from N(0,1)
#               X2: A binary baseline variable generated from Bern(0.5)
#               S1_Preference: Categorical variable; Preference for first stage treatment; takes values A, B, NP; A-Prefer A, B-Prefer B, NP-No preference
#               P1: Binary Variable: Indicator for whether an individual had a preference in stage 1. 1-had a preference, 0-had no preference
#               T1: Binary variable; Individuals treatment assigned at Stage 1; takes Values A(1) or B(-1). 
#               R: Binary variable; Individual's stage 1 response. 1 denotes response, 0 denotes no response
#               S2_Preference: Categorical variable; Preference for second stage treatment; takes values C, D, NP, 999. C-Prefer C, D-Prefer D, NP-No preference, 999-responders (stage 1 responders are not asked stage 2 preference)
#               P2: Binary Variable: Indicator for whether an individual had a preference in stage 2. 1-had a preference, 0-had no preference, 999-responder in stage 1
#               T2: Categorical variable; Individuals treatment assigned at Stage 2; takes values C(1), D(-1), 999-responders (stage 1 resopnder continue on initial treatment) 
#               Y: Binary variable; Individual's binary endpoint. 1 denotes response, 0 denotes no response
#               Treatment_Path: Treatment path the participant follows

#       (2) A weighted and replicated dataset to be used in a WRRM model
#             The dataset has the following variables:
#               ID: Numeric subject ID variable 
#               xalpha00: intercept alpha_00 in WRRM
#               xalpha01: intercept alpha_01 in WRRM
#               xalpha10: intercept alpha_10 in WRRM
#               xalpha11: intercept alpha_11 in WRRM
#               xbeta0: stage 1 randomized treatment beta0 in WRRM
#               xbeta1: stage 1 preference treatment beta1 in WRRM
#               xtheta0: stage 2 randomized treatment theta0 in WRRM
#               xtheta1: stage 2 preference treatment theta1 in WRRM
#               xgamma: stage 1 and 2 treatment interaction gamma in WRRM
#               Y: Binary variable; Individual's binary endpoint. 1 denotes response, 0 denotes no response
#               w: PRPP-SMART weight
#               Treatment_Path: Treatment path the participant follows
#          
#       (3) A dataset of sample sizes for each embedded DTR in the PRPP-SMART. That is, the number of participants (responders and non-responders) that constitute each embedded DTR in the simulated data
#Output:   Method: BJSM
#       (1) A two-stage PRPP-SMART dataset of N subjects 
#             The dataset has the following variables:
#               ID: Numeric subject ID variable 
#               X1: A continuous baseline variable generated from N(0,1)
#               X2: A binary baseline variable generated from Bern(0.5)
#               S1_Preference: Categorical variable; Preference for first stage treatment; takes values A, B, NP; A-Prefer A, B-Prefer B, NP-No preference
#               P1: Binary Variable: Indicator for whether an individual had a preference in stage 1. 1-had a preference, 0-had no preference
#               T1: Binary variable; Individuals treatment assigned at Stage 1; takes Values A(1) or B(-1). 
#               Y1: Binary variable; Individual's stage 1 response. 1 denotes response, 0 denotes no response
#               S2_Preference: Categorical variable; Preference for second stage treatment; takes values C, D, NP, 999. C-Prefer C, D-Prefer D, NP-No preference, 999-responders (stage 1 responders are not asked stage 2 preference)
#               P2: Binary Variable: Indicator for whether an individual had a preference in stage 2. 1-had a preference, 0-had no preference, 999-responder in stage 1
#               T2: Categorical variable; Individuals treatment assigned at Stage 2; takes values C(1), D(-1), 999-responders (stage 1 resopnder continue on initial treatment) 
#               Y2: Binary variable; Individual's stage 2 response. 1 denotes response, 0 denotes no response
#               Treatment_Path: Treatment path the participant follows

#       (2) A dataset formatted to be used in a BJSM 
#             The dataset has the following variables:
#               preference_stageI: Indicator for whether an individual had a preference in stage 1. 1-had a preference, 0-had no preference 
#               treatment_stageI: Stage 1 treatment coded for BJSM 1-A, 2-B
#               response_stageI: Individual's stage 1 response. 1 denotes response, 0 denotes no response
#               preference_stageII: Indicator for whether an individual had a preference in stage 2. 1-had a preference, 0-had no preference, 999-responder in stage 1
#               treatment_stageII: Stage 2 treatment coded for BJSM
# 1: Treatment Paths AA (A1A, A0A)
# 2: Treatment Paths BB (B1B, B0B)
# 3: Treatment Paths AC (A0C0, A0C1, A1C0, A1C1)
# 4: Treatment Paths AD (A0D0, A0D1, A1D0, A1D1) 
# 5: Treatment Paths BC (B0C0, B0C1, B1C0, B1C1)
# 6: Treatment Paths BD (B0D0, B0D1, B1D0, B1D1)
#               response_stageII: Individual's stage 2 response. 1 denotes response, 0 denotes no response
#               treatment_path: Treatment path the participant follows




###########################################################################################################################

generate_data <- function(N, pNP1, pTheta_A, pNP2, pTheta_C, pi_A, pi_B, pi_A1, pi_B1, pi_AC, pi_AD, pi_BC, pi_BD, pA0A, pB0B, pA1A, pB1B, pA1C1, pA1C0, pA1D1, pA1D0, pA0C1, pA0D1, pB1C1, pB1C0, pB1D1, pB1D0, pB0C1, pB0D1){
  
  # load libraries 
  library(Rlab)
  library(tidyverse)
  
  expit <- function(y) exp(y)/(1+exp(y))
  
  
  #Function: preference 
  #Purpose: These functions assign the preference of each individual based on the true propensity for exhibiting 
  #         a preference for treatment A, B, or having no preference in stage 1 and treatment C, D, or having no
  #         no preference in stage 2. Note, pref_probx is a Nx3 matrix with the
  #         true probabilities of preference for each individual
  
  
  preference_stage1 <- function(i){
    sample(c("A", "B", "NP"), size=1, replace = TRUE, prob=pref_prob1[i,])
  }
  
  preference_stage2 <- function(i){
    sample(c("C", "D", "NP"), size=1, replace = TRUE, prob=pref_prob2[i,])
  }
  
  
  # generate baseline covariates 
  X1<-rnorm(N,mean=0,sd=1)  # generate X1
  X2<-rbinom(N,1,0.5)  # generate X2: Home vs clinic/doctor’s office for treatment
  
  
  ### Preference Model: Stage 1 Propensity Scores for Preference ###
  #    Here we generate preference for first stage treatment for each patient. The true propensity for exhibiting a 
  #    preference for treatment A, B, or having no preference for each subject in the simulation population 
  #    is modeled using a logit model approach where first-stage propensity scores for preference are conditional on the 
  #    first-stage baseline covariates X1 and X2. In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of subjects exhibiting 
  #    no preference/preference for A/B
  
  
  ## No Preference ##
  
  # Coefficients on X1 and X2 
  a1 <- 0.2
  a2 <- 0.133 
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP1 <- function(a0)
  {
    alpha <- rbind(a0, a1, a2)
    lnp <- cbind(1, X1, X2) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP1) # want this to be close to 0
  }
  a0_star <- uniroot(search0_NP1, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  alpha <- c(a0_star, a1, a2)
  prob_NP <- expit(cbind(1, X1, X2) %*% alpha) # calculate no preference probability for each patient
  
  
  ## Model Theta: Prefer A among those with a preference ##
  pA_marginal <- pTheta_A*(1-pNP1) # marginal A probability in simulated data
  
  # Coefficients on X1 and X2 
  b1 <- 0.05
  b2 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_ThetaA <- function(b0)
  {
    beta <- rbind(b0, b1, b2)
    lptheta<- cbind(1, X1, X2) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta_A) # want this to be close to 0
  }
  
  b0_star <- uniroot(search0_ThetaA, c(-4,4))$root 
  beta <- c(b0_star, b1, b2)
  prob_A <- expit(cbind(1, X1, X2) %*% beta) * (1-prob_NP)
  prob_B <- (1-expit(cbind(1, X1, X2) %*% beta)) * (1-prob_NP) 
  
  ## generate preference of first stage treatment of each individual ##
  
  # put probabilities into a matrix
  pref_prob1 <- cbind(prob_A,prob_B,prob_NP) 
  colnames(pref_prob1) <- c("Prob Prefer A", "Prob Prefer B", "Prob No Preference") # each row is a subject's probability to prefer A, prefer B, or have no preference
  # apply(pref_prob1,1,sum) # check
  
  # sample based on probabilities 
  S1_Preference <- sapply(1:N, preference_stage1) # stage 1 preference 
  
  ### CHECKS ###
  # sum(prop.table(table(S1_Preference))) # check to make sure prob_A, prob_B, prob_NP sum to 1
  # prop.table(table(S1_Preference)) # check if observed proportions close to targets
  # length(which(S1_Preference == "A")) / length(which(S1_Preference == "A" | S1_Preference == "B")) # should be close to pTheta_A
  
  
  ### Stage 1 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 1 assigned treatment
  #    based on the true propensities for exhibiting a preference for treatment A, B, or having no preference
  #    Prefer A: get A
  #    Prefer B: get B
  #    No preference: equal probability of getting treatment A/B
  
  
  # generate actual stage 1 treatment that each individual received
  # 1-treatment A, 0-treatment B
  T1 <- 
    (S1_Preference == "A") * 1 + # prefer A get A
    (S1_Preference == "B") * 0 + # prefer B get B
    (S1_Preference == "NP") * sample(c(1, 0), N, replace = T, prob=c(0.5,0.5)) # no preference randomly assign A/B
  
  
  ### Generate Preference Indicator ###
  #    Here we create an indicator, P1, for whether an individual had a preference in stage 1 or not. Specifically, 1
  #    if an individual had a preference in stage 1, 0 if an individual had no preference in stage 1
  
  P1 <- ifelse(S1_Preference == "A" | S1_Preference == "B", 1, 0)
  
  ### Generate Stage 1 Response variable ###
  #    For each subject we generate a binary response where 1 indicates the subjects responds
  #    to the assigned Stage 1 treatment and 0 indicates the subject does not response to the assigned 
  #    stage 1 treatment. 
  
  # Generate stage 1 response variable 
  response_stageI <- rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  response_stageI[P1 == 1 & T1==1] <- rbinom(sum(P1 == 1 & T1==1),1, pi_A1) # got preferred treatment A in Stage 1
  response_stageI[P1 == 0 & T1==1] <- rbinom(sum(P1 == 0 & T1==1),1, pi_A) # got randomized treatment A in Stage 1
  response_stageI[P1 == 1 & T1==0] <- rbinom(sum(P1 == 1 & T1==0),1, pi_B1) # got preferred treatment B in Stage 1
  response_stageI[P1 == 0 & T1==0] <- rbinom(sum(P1 == 0 & T1==0),1, pi_B) # got randomized treatment B in Stage 1
  
  ### Prefernce Model: Stage 2 Propensity Scores for Preference ###
  #    Here we generate preference of second stage treatment only for first-stage non-responders. Patients who respond
  #    continue on the stage 1 treatment. The true propensity for exhibiting a 
  #    preference for treatment C, D, or having no preference for non-responders in the simulation population 
  #    is modeled using a logit model approach where second-stage propensity scores for preference are conditional on whether
  #    the patient had a preference in stage 1. In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of non-responder subjects 
  #    exhibiting no preference/preference for C/D
  
  # find non-responders 
  nr_index <- which(response_stageI == 0)
  
  ## No Preference ##
  
  # Coefficients on P1 
  c1 <- -0.1
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP2 <- function(c0){
    alpha <- rbind(c0, c1)
    lnp <- cbind(1, P1[nr_index]) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP2) # want this to be close to 0
  }
  c0_star <- uniroot(search0_NP2, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  gamma <- c(c0_star, c1)
  prob_NP2 <- expit(cbind(1, P1[nr_index]) %*% gamma) # calculate no preference probability for each patient
  
  
  ## Model Theta: Prefer C among those with a preference ##
  pC_marginal <- pTheta_C*(1-pNP2) # marginal C probability in simulated data
  
  # Coefficients on S1_P
  d1 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_ThetaC <- function(d0)
  {
    beta <- rbind(d0, d1)
    lptheta<- cbind(1, P1[nr_index]) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta_C) # want this to be close to 0
  }
  
  d0_star <- uniroot(search0_ThetaC, c(-4,4))$root 
  phi <- c(d0_star, d1)
  prob_C <- expit(cbind(1, P1[nr_index]) %*% phi) * (1-prob_NP2)
  prob_D <- (1-expit(cbind(1, P1[nr_index]) %*% phi)) * (1-prob_NP2) 
  
  ## generate preference of second stage treatment for non-responders ##
  
  # put probabilities into a matrix
  pref_prob2 <- cbind(prob_C,prob_D,prob_NP2) 
  colnames(pref_prob2) <- c("Prob Prefer C", "Prob Prefer D", "Prob No Preference") # each row is a subject's probability to prefer C, prefer D, or have no preference
  # apply(pref_prob2,1,sum) # check
  
  # sample based on probabilities 
  n2 <- length(nr_index)
  S2_Preference <- sapply(1:n2, preference_stage2) # stage 2 preference only generated for stage 1 non-responders
  
  ### CHECKS ###
  #sum(prop.table(table(S2_Preference))) # check to make sure prob_C, prob_D, prob_NP2 sum to 1
  # prop.table(table(S2_Preference)) # check if observed proportions close to targets
  # length(which(S2_Preference == "C")) / length(which(S2_Preference == "C" | S2_Preference == "D")) # should be close to pTheta_C
  
  
  # create final stage 2 preference variable 
  S2_Preference_final <- rep(999, N) # 999 for stage 1 responders since they continue on stage 1 treatment 
  
  S2_Preference_final[nr_index] <- S2_Preference
  
  
  ### Stage 2 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 2 assigned treatment.
  #    For non-responders this is based off the true propensities for exhibiting a preference for treatment C, D, or having no preference
  #    For responders they continue on their stage 1 assigned treatment.
  #    Non-responders:
  #      Prefer C: get C
  #      Prefer D: get D
  #      No preference: equal probability of getting treatment C/D
  
  
  # generate actual stage 2 treatment that each individual received
  # 1-treatment C, 0-treatment D
  T2 <- c() 
  for(i in 1:N){
    
    if (response_stageI[i]==1){ # responds to Stage 1 treatment
      T2[i] <- T1[i] # continues on initial treatment received 
      
    }
    if (response_stageI[i] == 0 & S2_Preference_final[i] == "NP"){# no response to Stage 1 treatment & has no preference in stage 2
      T2[i] <- rbern(1,0.5) # equal chance for treatment C/D
    }
    
    if (response_stageI[i] == 0 &  S2_Preference_final[i] == "C"){# no response to Stage 1 treatment & prefer C
      T2[i] <- 1 # get preferred treatment 
    }
    
    if (response_stageI[i] == 0 &  S2_Preference_final[i] == "D"){# no response to Stage 1 treatment & prefer D
      T2[i] <- 0 # get preferred treatment 
    }
  }
  
  
  ### Generate Stage 2 Preference Indicator ###
  #    Here we create an indicator, S2_P, for whether an individual had a preference in stage 2 or not. Specifically, 1
  #    if an individual had a preference in stage 2, 0 if an individual had no preference in stage 2
  
  P2 <- ifelse(S2_Preference_final == "C" | S2_Preference_final == "D", 1, 0)
  
  P2[which(response_stageI==1)] <- 999
  
  # Relabel T1 and T2 with actual treatment 
  T1[T1==1] <- "A"
  T1[T1==0] <- "B"
  
  T2[which(response_stageI == 1)] <- T1[which(response_stageI == 1)] # assign stage 1 responders stage 1 treatment
  T2[T2==1] <- "C"
  T2[T2==0] <- "D"
  
  
  # Create final treatment_stage1 and treatment_stage2 variables for BJSM code
  # 1: Treatment Paths AA (A1A or A0A)
  # 2: Treatment Paths BB (B1B or B0B)
  # 3: Treatment Paths AC (A0C0, A0C1, A1C0, A1C1)
  # 4: Treatment Paths AD (A0D0, A0D1, A1D0, A1D1) 
  # 5: Treatment Paths BC (B0C0, B0C1, B1C0, B1C1)
  # 6: Treatment Paths BD (B0D0, B0D1, B1D0, B1D1)
  
  treatment_stageI <- rep(NA, N)
  treatment_stageI[which(T1 == "A")] <- 1
  treatment_stageI[which(T1 == "B")] <- 2
  
  treatment_stageII <- rep(NA, N)
  treatment_stageII[which(T1 == "A" & response_stageI == 1 & T2 == "A")] <- 1
  treatment_stageII[which(T1 == "B" & response_stageI == 1 & T2 == "B")] <- 2
  treatment_stageII[which(T1 == "A" & response_stageI == 0 & T2 == "C")] <- 3
  treatment_stageII[which(T1 == "A" & response_stageI == 0 & T2 == "D")] <- 4
  treatment_stageII[which(T1 == "B" & response_stageI == 0 & T2 == "C")] <- 5
  treatment_stageII[which(T1 == "B" & response_stageI == 0 & T2 == "D")] <- 6
  
  
  treatmentCD_stageII <- rep(0,N)
  treatmentCD_stageII[which(response_stageI == 0 & T2 == "C")] <- 3
  treatmentCD_stageII[which(response_stageI == 0 & T2 == "D")] <- 4
  
  ### Generate Trial Outcome Variable ###
  #    For each subject we generate the final binary trial outcome given their study path where 1 indicates the subjects responds
  #    0 indicates the subject does not respond.
  
  
  response_stageII <- rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  response_stageII[P1 == 1 & T1=="A" & response_stageI==1]<-rbinom(sum(P1 == 1 & T1=="A" & response_stageI==1),1,pA1A)
  response_stageII[P1 == 0 & T1=="A" & response_stageI==1]<-rbinom(sum(P1 == 0 & T1=="A" & response_stageI==1),1,pA0A)
  response_stageII[P1 == 1 & T1=="B" & response_stageI==1]<-rbinom(sum(P1 == 1 & T1=="B" & response_stageI==1),1,pB1B)
  response_stageII[P1 == 0 & T1=="B" & response_stageI==1]<-rbinom(sum(P1 == 0 & T1=="B" & response_stageI==1),1,pB0B)
  
  # paths that begin with prefer A in stage 1
  response_stageII[P1 == 1 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "C"),1,pA1C1)
  response_stageII[P1 == 1 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "C"),1,pA1C0)
  response_stageII[P1 == 1 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "D"),1,pA1D1)
  response_stageII[P1 == 1 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "D"),1,pA1D0)
  
  # paths that begin with randomize A in stage 1
  response_stageII[P1 == 0 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "C"),1,pA0C1)
  response_stageII[P1 == 0 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "C"),1,pi_AC)
  response_stageII[P1 == 0 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="A" & response_stageI==0 & P2 == 1 & T2 == "D"),1,pA0D1)
  response_stageII[P1 == 0 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="A" & response_stageI==0 & P2 == 0 & T2 == "D"),1,pi_AD)
  
  # paths that begin with prefer B in stage 1
  response_stageII[P1 == 1 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "C"),1,pB1C1)
  response_stageII[P1 == 1 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "C"),1,pB1C0)
  response_stageII[P1 == 1 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "D"),1,pB1D1)
  response_stageII[P1 == 1 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "D"),1,pB1D0)
  
  
  # paths that begin with randomize B in stage 1
  response_stageII[P1 == 0 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "C"),1,pB0C1)
  response_stageII[P1 == 0 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "C"),1,pi_BC)
  response_stageII[P1 == 0 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="B" & response_stageI==0 & P2 == 1 & T2 == "D"),1,pB0D1)
  response_stageII[P1 == 0 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="B" & response_stageI==0 & P2 == 0 & T2 == "D"),1,pi_BD)
  
  ### Create treatment path variable ###
  treatment_path <- rep(NA, N)
  
  # responder paths
  treatment_path[which(P1==1 & T1=="A" & response_stageI==1)] <- "A1A"
  treatment_path[which(P1==0 & T1=="A" & response_stageI==1)] <- "A0A"
  treatment_path[which(P1==1 & T1=="B" & response_stageI==1)] <- "B1B"
  treatment_path[which(P1==0 & T1=="B" & response_stageI==1)] <- "B0B"
  
  # non-responder paths
  treatment_path[which(P1==1 & T1=="A" & response_stageI==0 & P2==1 & T2=="C")] <- "A1C1"
  treatment_path[which(P1==1 & T1=="A" & response_stageI==0 & P2==0 & T2=="C")] <- "A1C0"
  treatment_path[which(P1==1 & T1=="A" & response_stageI==0 & P2==1 & T2=="D")] <- "A1D1"
  treatment_path[which(P1==1 & T1=="A" & response_stageI==0 & P2==0 & T2=="D")] <- "A1D0"
  
  treatment_path[which(P1==1 & T1=="B" & response_stageI==0 & P2==1 & T2=="C")] <- "B1C1"
  treatment_path[which(P1==1 & T1=="B" & response_stageI==0 & P2==0 & T2=="C")] <- "B1C0"
  treatment_path[which(P1==1 & T1=="B" & response_stageI==0 & P2==1 & T2=="D")] <- "B1D1"
  treatment_path[which(P1==1 & T1=="B" & response_stageI==0 & P2==0 & T2=="D")] <- "B1D0"
  
  treatment_path[which(P1==0 & T1=="A" & response_stageI==0 & P2==1 & T2=="C")] <- "A0C1"
  treatment_path[which(P1==0 & T1=="A" & response_stageI==0 & P2==0 & T2=="C")] <- "A0C0"
  treatment_path[which(P1==0 & T1=="A" & response_stageI==0 & P2==1 & T2=="D")] <- "A0D1"
  treatment_path[which(P1==0 & T1=="A" & response_stageI==0 & P2==0 & T2=="D")] <- "A0D0"
  
  treatment_path[which(P1==0 & T1=="B" & response_stageI==0 & P2==1 & T2=="C")] <- "B0C1"
  treatment_path[which(P1==0 & T1=="B" & response_stageI==0 & P2==0 & T2=="C")] <- "B0C0"
  treatment_path[which(P1==0 & T1=="B" & response_stageI==0 & P2==1 & T2=="D")] <- "B0D1"
  treatment_path[which(P1==0 & T1=="B" & response_stageI==0 & P2==0 & T2=="D")] <- "B0D0"
  
  
  # Create PRPP-SMART data 
  data_output <- 
    tibble(
      id = 1:N,
      X1 = X1,
      X2 = X2,
      S1_Preference = S1_Preference,
      P1 = P1,
      T1 = T1,
      Y1 = response_stageI,
      S2_Preference = S2_Preference_final,
      P2 = P2,
      T2 = T2,
      Y2 = response_stageII,
      Treatment_Path = treatment_path
    )
  
  # Create dataset formatted to be used in BJSM
  bjsm_data <- data.frame(preference_stageI = P1, # stage 1 preference indicator # S1_P --> P1
                          treatment_stageI, # stage 1 treatment coded for BJSM 1-A, 2-B
                          response_stageI,  # Response at end of stage 1 
                          preference_stageII = P2, # stage 2 preference indicator among non-responders S2_P --> P2
                          treatment_stageII, # stage 2 treatment coded for BJSM 1, 2, 3, 4, 5, 6
                          treatmentCD_stageII, # stage 2 treatment ended with C (3) or D (4)
                          response_stageII, # Response at end of stage 2 
                          treatment_path) # trial path participant follows
  
  return(list(data_output, bjsm_data))
  
}
