# Arguments:
# sample_size: total sample size in PRPP-SMART study
# stage2_NR_outcome_prob: probability of stage 2 outocme for each of 4 non-responder embedded treatment sequences in the following order (AC, AD, BC, BD)
# stage1_trtA_outcome_prob: probability of outcome at the end of stage-1 to treatment T1=A
# stage1_trtB_outcome_prob: probability of outcome at the end of stage-1 to treatment T1=B
# linkage_params: values of 11 linkage parameters in the following order (beta1A, beta1B, alpha1A, alpha1B, alpha2C, alpha2D, alphaAC, alphaAD, alphaBC, alphaBD)
# prior_shape_params: prior shape parameters for pi_j and pi_jk in the order pi_A(a,b), pi_B(a,b), pi_AC(a,b), pi_AD(a,b), pi_BC(a,b), pi_BD(a,b)
# threshold: threshold of log-OR for exclusion from the set of best
# alpha_type1: probability of excluding true best EDTR from the set of best  



powerFunction <- function(sample_size,
                          Treatment_Preference_stage1_prob,
                          Treatment_Preference_stage2_prob,
                          PreferA_prob, 
                          PreferC_prob,
                          stage2_NR_outcome_prob,
                          stage1_trtA_outcome_prob,
                          stage1_trtB_outcome_prob,
                          linkage_params,
                          prior_shape_params,
                          threshold,
                          alpha_type1 = 0.05){
  
  
  # progress bar 
  pb <- txtProgressBar(min=0,
                       max=1000,
                       initial = 0,
                       style=3) 
  
  upper_limit <- array(NA, dim = c(1000, 10, 4)) # where to store CI upper limit for indifference DTRs
  
  pNP1 <- 1 - Treatment_Preference_stage1_prob
  pNP2 <- 1 - Treatment_Preference_stage2_prob
  N <- sample_size
  pA <- PreferA_prob
  pB <- 1-pA 
  pC <- PreferC_prob 
  pD <- 1-pC 
  
  # linkage parameters
  beta1A <- linkage_params[1] # responders to A second stage response 
  beta1B <- linkage_params[2] # responders to B second stage response
  alpha1A <- linkage_params[3] 
  alpha1B <- linkage_params[4]
  alpha2C <- linkage_params[5]
  alpha2D <- linkage_params[6]

  
  # First stage treatment response rates
  pi_A <- stage1_trtA_outcome_prob                                          # stage 1 response rate to randomize A: Pr(R=1|T1=A,P1=0)
  pi_B <- stage1_trtB_outcome_prob                                         # stage 1 response rate to randomize B: Pr(R=1|T1=B,P1=0)
  pi_A1 <- pi_A*alpha1A                                # stage 1 response rate to prefer A: Pr(R=1|T1=A,P1=1)
  pi_B1 <- pi_B*alpha1B                                # stage 1 response rate to prefer B: Pr(R=1|T1=B,P1=1)
  # Second stage treatment response rates
  pi_AC <- stage2_NR_outcome_prob[1]                                         # Second stage response rate of non-responders to randomized A who receive randomized C in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=C)        
  pi_AD <- stage2_NR_outcome_prob[2]                                         # Second stage response rate of non-responders to randomized A who receive randomized D in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=D)
  pi_BC <- stage2_NR_outcome_prob[3]                                          # Second stage response rate of non-responders to randomized B who receive randomized C in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=C)
  pi_BD <- stage2_NR_outcome_prob[4]                                          # Second stage response rate of non-responders to randomized B who receive randomized D in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=D)
  pA0A <- pi_A * beta1A                                # Second stage response rate of responders to randomized A
  pB0B <- pi_B * beta1B                                # Second stage response rate of responders to randomized B
  pA1A <- pi_A * alpha1A * beta1A                      # Second stage response rate of responders to preferred A
  pB1B <- pi_B * alpha1B * beta1B                      # Second stage response rate of responders to preferred B
  pA0C1 <- alpha2C * pi_AC                             # Second stage response rate of non-responders to randomized A who receive preferred C in the second stage
  pA0D1 <- alpha2D * pi_AD                             # Second stage response rate of non-responders to randomized A who receive preferred D in the second stage
  pA1C0 <- alpha1A * pi_AC                             # Second stage response rate of non-responders to preferred A who receive randomized C in the second stage
  pA1D0 <- alpha1A * pi_AD                             # Second stage response rate of non-responders to preferred A who receive randomized D in the second stage
  pA1C1 <-  alpha1A * alpha2C * pi_AC                            # Second stage response rate of non-responders to preferred A who receive preferred C in the second stage
  pA1D1 <- alpha1A * alpha2D * pi_AD                             # Second stage response rate of non-responders to preferred A who receive preferred D in the second stage
  pB0C1 <- alpha2C * pi_BC                             # Second stage response rate of non-responders to randomized B who receive preferred C in the second stage
  pB0D1 <- alpha2D * pi_BD                             # Second stage response rate of non-responders to randomized B who receive preferred D in the second stage
  pB1C0 <- alpha1B * pi_BC                             # Second stage response rate of non-responders to preferred B who receive randomized C in the second stage
  pB1D0 <- alpha1B * pi_BD                             # Second stage response rate of non-responders to preferred B who receive randomized D in the second stage
  pB1C1 <- alpha1B * alpha2C * pi_BC                             # Second stage response rate of non-responders to preferred B who receive preferred C in the second stage
  pB1D1 <-  alpha1B * alpha2D * pi_BD                            # Second stage response rate of non-responders to preferred B who receive preferred D in the second stage
  
  
  # prior shape parameters on randomized response rates
  a_A <- prior_shape_params[1] # pi_A ~ Beta(a_A, b_A)
  b_A <- prior_shape_params[2]
  a_B <- prior_shape_params[3] # pi_B ~ Beta(a_B, b_B)
  b_B <- prior_shape_params[4]
  
  c_AC <- prior_shape_params[5] # pi_AC ~ Beta(c_AC, d_AC)
  d_AC <- prior_shape_params[6]
  c_AD <- prior_shape_params[7] # pi_AD ~ Beta(c_AD, d_AD)
  d_AD <- prior_shape_params[8]
  c_BC <- prior_shape_params[9] # pi_BC ~ Beta(c_BC, d_BC)
  d_BC <- prior_shape_params[10]
  c_BD <- prior_shape_params[11] # pi_BD ~ Beta(c_BD, d_BD)
  d_BD <- prior_shape_params[12]
  
  # Initialize the number of successful simulations
  successful_simulations <- 0
  
  # Set the desired number of simulations
  total_simulations <- 1000
  iteration_number <- 0
  
  while (successful_simulations < total_simulations){
    
    
    iteration_number <- iteration_number + 1
    
    
    ##### Compute summary statistics (data generation) ####
    # number of subjects randomized to A
    T_A0 <- ceiling(N*pNP1*0.5)
    # number of subjects prefer A
    T_A1 <- ceiling(N*(1-pNP1)*pA)
    # number of subjects randomized to B
    T_B0 <- ceiling(N*pNP1*0.5)
    # number of subjects prefer B
    T_B1 <- ceiling(N*(1-pNP1)*(pB))
    
    
    # number of stage 1 responders to A P1 = 0
    r_A0 <- rbinom(n=1, size = T_A0, pi_A)
    # number of stage 1 responders to B P1 = 0
    r_B0 <- rbinom(n=1, size = T_B0, pi_B)
    # number of stage 1 responders to A P1 = 1
    s_A1 <- rbinom(n=1, size = T_A1, pi_A1)
    # number of stage 1 responders to B P1 = 1
    s_B1 <- rbinom(n=1, size = T_B1, pi_B1)
    
    # number of individuals that respond to A in stages 1 and 2 with P1 = 0
    X_AA0 <- rbinom(n=1, size = r_A0, pA0A)
    # number of individuals that respond to B in stages 1 and 2 with P1 = 0
    X_BB0 <- rbinom(n=1, size = r_B0, pB0B)
    # number of individuals that respond to A in stages 1 and 2 with P1 = 1
    Y_AA1 <- rbinom(n=1, size = s_A1, pA1A)
    # number of individuals that respond to B in stages 1 and 2 with P1 = 1
    Y_BB1 <- rbinom(n=1, size = s_B1, pB1B)
    
    # number of non-responders to stage 1 treatment A that receive C in stage 2 with P1 = 0 and P2 = 0
    n_AC00 <- ceiling(N*pNP1*0.5*(1-pi_A)*pNP2*0.5)
    # number of non-responders to stage 1 treatment A that receive D in stage 2 with P1 = 0 and P2 = 0
    n_AD00 <- ceiling(N*pNP1*0.5*(1-pi_A)*pNP2*0.5)
    # number of non-responders to stage 1 treatment B that receive C in stage 2 with P1 = 0 and P2 = 0
    n_BC00 <- ceiling(N*pNP1*0.5*(1-pi_B)*pNP2*0.5)
    # number of non-responders to stage 1 treatment B that receive D in stage 2 with P1 = 0 and P2 = 0
    n_BD00 <- ceiling(N*pNP1*0.5*(1-pi_B)*pNP2*0.5)
    
    # number of non-responders to stage 1 treatment A that receive C in stage 2 with P1 = 0 and P2 = 1
    m_AC01 <- ceiling(N*pNP1*0.5*(1-pi_A)*(1-pNP2)*pC)
    # number of non-responders to stage 1 treatment A that receive D in stage 2 with P1 = 0 and P2 = 1
    m_AD01 <- ceiling(N*pNP1*0.5*(1-pi_A)*(1-pNP2)*(pD))
    # number of non-responders to stage 1 treatment B that receive C in stage 2 with P1 = 0 and P2 = 1
    m_BC01 <- ceiling(N*pNP1*0.5*(1-pi_B)*(1-pNP2)*pC)
    # number of non-responders to stage 1 treatment B that receive D in stage 2 with P1 = 0 and P2 = 1
    m_BD01 <- ceiling(N*pNP1*0.5*(1-pi_B)*(1-pNP2)*(pD))
    
    # number of non-responders to stage 1 treatment A that receive C in stage 2 with P1 = 1 and P2 = 0
    l_AC10 <-  ceiling(N*(1-pNP1)*pA*(1-pi_A1)*pNP2*0.5)
    # number of non-responders to stage 1 treatment A that receive D in stage 2 with P1 = 1 and P2 = 0
    l_AD10 <- ceiling(N*(1-pNP1)*pA*(1-pi_A1)*pNP2*0.5)
    # number of non-responders to stage 1 treatment B that receive C in stage 2 with P1 = 1 and P2 = 0
    l_BC10 <- ceiling(N*(1-pNP1)*(pB)*(1-pi_B1)*pNP2*0.5)
    # number of non-responders to stage 1 treatment B that receive D in stage 2 with P1 = 1 and P2 = 0
    l_BD10 <- ceiling(N*(1-pNP1)*(pB)*(1-pi_B1)*pNP2*0.5)
    
    # number of non-responders to stage 1 treatment A that receive C in stage 2 with P1 = 1 and P2 = 1
    h_AC11 <- ceiling(N*(1-pNP1)*pA*(1-pi_A1)*(1-pNP2)*pC)
    # number of non-responders to stage 1 treatment A that receive D in stage 2 with P1 = 1 and P2 = 1
    h_AD11 <- ceiling(N*(1-pNP1)*pA*(1-pi_A1)*(1-pNP2)*(pD))
    # number of non-responders to stage 1 treatment B that receive C in stage 2 with P1 = 1 and P2 = 1
    h_BC11 <- ceiling(N*(1-pNP1)*(pB)*(1-pi_B1)*(1-pNP2)*pC)
    # number of non-responders to stage 1 treatment B that receive D in stage 2 with P1 = 1 and P2 = 1
    h_BD11 <-ceiling(N*(1-pNP1)*(pB)*(1-pi_B1)*(1-pNP2)*(pD))
    
    # number of individuals that do not respond to A in stage 1 but respond to C in stage 2 with P1 = 0 and P2 = 0
    X_AC00 <- rbinom(n=1, size = n_AC00, pi_AC)
    # number of individuals that do not respond to A in stage 1 but respond to D in stage 2 with P1 = 0 and P2 = 0
    X_AD00 <- rbinom(n=1, size = n_AD00, pi_AD)
    # number of individuals that do not respond to B in stage 1 but respond to C in stage 2 with P1 = 0 and P2 = 0
    X_BC00 <- rbinom(n=1, size = n_BC00, pi_BC)
    # number of individuals that do not respond to B in stage 1 but respond to D in stage 2 with P1 = 0 and P2 = 0
    X_BD00 <- rbinom(n=1, size = n_BD00, pi_BD)
    
    # number of individuals that do not respond to A in stage 1 but respond to C in stage 2 with P1 = 0 and P2 = 1
    Y_AC01 <- rbinom(n=1, size = m_AC01, pA0C1)
    # number of individuals that do not respond to A in stage 1 but respond to D in stage 2 with P1 = 0 and P2 = 1
    Y_AD01 <- rbinom(n=1, size = m_AD01, pA0D1)
    # number of individuals that do not respond to B in stage 1 but respond to C in stage 2 with P1 = 0 and P2 = 1
    Y_BC01 <- rbinom(n=1, size = m_BC01, pB0C1)
    # number of individuals that do not respond to B in stage 1 but respond to D in stage 2 with P1 = 0 and P2 = 1
    Y_BD01 <- rbinom(n=1, size = m_BD01, pB0D1)
    
    # number of individuals that do not respond to A in stage 1 but respond to C in stage 2 with P1 = 1 and P2 = 0
    Z_AC10 <- rbinom(n=1, size = l_AC10, pA1C0)
    # number of individuals that do not respond to A in stage 1 but respond to D in stage 2 with P1 = 1 and P2 = 0
    Z_AD10 <- rbinom(n=1, size = l_AD10, pA1D0)
    # number of individuals that do not respond to B in stage 1 but respond to C in stage 2 with P1 = 1 and P2 = 0
    Z_BC10 <- rbinom(n=1, size = l_BC10, pB1C0)
    # number of individuals that do not respond to B in stage 1 but respond to D in stage 2 with P1 = 1 and P2 = 0
    Z_BD10 <- rbinom(n=1, size = l_BD10, pB1D0)
    
    # number of individuals that do not respond to A in stage 1 but respond to C in stage 2 with P1 = 1 and P2 = 1
    W_AC11 <- rbinom(n=1, size = h_AC11, pA1C1)
    # number of individuals that do not respond to A in stage 1 but respond to D in stage 2 with P1 = 1 and P2 = 1
    W_AD11 <- rbinom(n=1, size = h_AD11, pA1D1)
    # number of individuals that do not respond to B in stage 1 but respond to C in stage 2 with P1 = 1 and P2 = 1
    W_BC11 <- rbinom(n=1, size = h_BC11, pB1C1)
    # number of individuals that do not respond to B in stage 1 but respond to D in stage 2 with P1 = 1 and P2 = 1
    W_BD11 <- rbinom(n=1, size = h_BD11, pB1D1)
    
    # Check if any variable is zero
    if (any(c(n_AC00, n_AD00, n_BC00, n_BD00,
              m_AC01, m_AD01, m_BC01, m_BD01,
              l_AC10,l_AD10, l_BC10, l_BD10,
              h_AC11, h_AD11, h_BC11, h_BD11,
              r_A0, s_A1, r_B0, s_B1) == 0)) {
      # Skip this iteration
      next
      
    }
    
    # Increment the successful simulation counter
    successful_simulations <- successful_simulations + 1
    
    shape1AC <- c_AC + X_AC00 + (Y_AC01/alpha2C) + (Z_AC10/alpha1A) + (W_AC11/(alpha1A*alpha2C)) 
    shape2AC <- d_AC + (n_AC00 - X_AC00) + (m_AC01 - (Y_AC01/alpha2C)) + (l_AC10 - (Z_AC10/alpha1A)) + (h_AC11 - (W_AC11/(alpha1A*alpha2C)))
    
    shape1AD <- c_AD + X_AD00 + (Y_AD01/alpha2D) + (Z_AD10/alpha1A) + (W_AD11/(alpha1A*alpha2D))
    shape2AD <- d_AD + (n_AD00 - X_AD00) + (m_AD01 - (Y_AD01/alpha2D)) + (l_AD10 - (Z_AD10/alpha1A)) + (h_AD11 - (W_AD11/(alpha1A*alpha2D)))
    
    shape1BC <- c_BC + X_BC00 + (Y_BC01/alpha2C) + (Z_BC10/alpha1B) + (W_BC11/(alpha1B*alpha2C)) 
    shape2BC <- d_BC + (n_BC00 - X_BC00) + (m_BC01 - (Y_BC01/alpha2C)) + (l_BC10 - (Z_BC10/alpha1B)) + (h_BC11 - (W_BC11/(alpha1B*alpha2C)))
    
    shape1BD <- c_BD + X_BD00 + (Y_BD01/alpha2D) + (Z_BD10/alpha1B) + (W_BD11/(alpha1B*alpha2D)) 
    shape2BD <- d_BD + (n_BD00 - X_BD00) + (m_BD01 - (Y_BD01/alpha2D)) + (l_BD10 - (Z_BD10/alpha1B)) + (h_BD11 - (W_BD11/(alpha1B*alpha2D)))
    
    
    shape1ThetaA <- a_A + X_AA0 + (Y_AA1/alpha1A)
    shape2ThetaA <- b_A + (r_A0 - X_AA0) + (s_A1 - (Y_AA1/alpha1A))
    
    shape1ThetaB <- a_B + X_BB0 + (Y_BB1/alpha1B)
    shape2ThetaB <- b_B + (r_B0 - X_BB0) + (s_B1 - (Y_BB1/alpha1B))
    
    for (j in 1:10){ # number of parameter draws for Monte Carlo integration for the parameter integral 
      
      # Draw 1000 draws from the posterior of the probability of outcome at the end of the study for each of the 6 indifference embedded treatment sequences.
      # A0A (thetaA), B0B (thetaB), AC, AD, BC, BD
      
      # posterior pi_jk
      piAC_draws <- rbeta(n = 1000, shape1 = shape1AC, shape2 = shape2AC)
      piAD_draws <- rbeta(n = 1000, shape1 = shape1AD, shape2 = shape2AD)
      piBC_draws <- rbeta(n = 1000, shape1 = shape1BC, shape2 = shape2BC)
      piBD_draws <- rbeta(n = 1000, shape1 = shape1BD, shape2 = shape2BD)
      
      # posterior theta_j
      thetaA_draws <- rbeta(n = 1000, shape1 = shape1ThetaA, shape2 = shape2ThetaA)
      thetaB_draws <- rbeta(n = 1000, shape1 = shape1ThetaB, shape2 = shape2ThetaB)
      
      # Draw 1000 draws from the posterior of the probability of outcome at end of stage 1 for each of the indifference stage 1 treatment sequences 
      # posterior pi_j
      A0_draws <- rbeta(n = 1000, shape1 = r_A0 + a_A, shape2 = T_A0 - r_A0 + b_A)
      B0_draws <- rbeta(n = 1000, shape1 = r_B0 + a_B, shape2 = T_B0 - r_B0 + b_B)
      

      # Compute embedded DTR response probabilities using Robins' G-computation method
      thetadraws <- cbind(thetaA_draws * A0_draws + piAC_draws * (1 - A0_draws), # AAC00
                          thetaA_draws * A0_draws + piAD_draws * (1 - A0_draws), # AAD00
                          thetaB_draws * B0_draws + piBC_draws * (1 - B0_draws), # BBC00
                          thetaB_draws * B0_draws + piBD_draws * (1 - B0_draws)) # BBD00
      
      
      
      # Perform Bayesian MCB
      logORThreshold <- LogOR(response_prob = c(pA0A,pB0B,pi_AC,pi_AD,pi_BC,pi_BD),
                              stage_one_trt_one_response_prob = pi_A,
                              stage_one_trt_two_response_prob = pi_B)
      thetadraws_log_odds <- log(thetadraws / (1 - thetadraws))
      max_odds_ind <- which.max(colMeans(thetadraws_log_odds))
      
      log_OR_matrix <- thetadraws_log_odds - matrix(thetadraws_log_odds[, max_odds_ind], nrow = 1000, ncol = 4)
      rank_matrix <- apply(log_OR_matrix[, -max_odds_ind], 2, rank, ties.method = 'random')
      rank_max <- apply(rank_matrix, 1, max)
      
      new_dat <- apply(log_OR_matrix[,], 2, sort)
      ranks_quantile <- floor(quantile(rank_max, probs = 1 - alpha_type1))
      
      upper_limit[successful_simulations,j,] <- new_dat[ranks_quantile,]
      
      #CI_upper <- apply(new_dat[1:ranks_quantile,], 2, max) # CI upper limit
      
      #upper_limit[i, j, ] <- CI_upper # update array storing CI upper limit for indifference DTRs
      
      
      
    }
    
    # Update progress bar
    setTxtProgressBar(pb, successful_simulations)
  }
  
  
  close(pb)
  
  rejection_indices <- which(abs(logORThreshold)>=threshold)
  
  
  if (length(rejection_indices)==1) {
    power <- mean(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0)
    
  } else {
    power <- mean(apply(apply(upper_limit, 3, function(x) x)[, rejection_indices] < 0, 1, prod))
  }
  
  # check if set of best is correct 
  set_of_best_indicies <- which(abs(logORThreshold)<threshold)
  
  if (length(set_of_best_indicies)==1) {
    setofbestcorrect <- mean(apply(upper_limit, 3, function(x) x)[, set_of_best_indicies] >= 0)
    
  } else {
    setofbestcorrect <- mean(apply(apply(upper_limit, 3, function(x) x)[, set_of_best_indicies] >= 0, 1, prod))
  }
  
  # best index
  best_index <- which(abs(logORThreshold)==0)
  bestinset <- mean(apply(upper_limit, 3, function(x) x)[, best_index] >= 0)
  
  
  
  return(list(power, setofbestcorrect, bestinset))
  
  
  
}