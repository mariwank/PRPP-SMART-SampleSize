# Sub-scenario linkage parameters not equal to 1

# Data gen parameters 
beta1_link <- c(1,1)                                 # Beta1A, Beta1B 
alphaP_link <- c(1.1, 1.02, 1.05, 1.2)               # alpha1A, alpha1B, alpha2C, alpha2D

# First stage treatment outcome rates
pi_A <- 0.6                                          # stage 1 response rate to randomize A: Pr(Y1=1|T1=A,P1=0)
pi_B <- 0.45                                         # stage 1 response rate to randomize B: Pr(Y1=1|T1=B,P1=0)
pi_A1 <- pi_A*alphaP_link[1]                         # stage 1 response rate to prefer A: Pr(Y1=1|T1=A,P1=1)
pi_B1 <- pi_B*alphaP_link[2]                         # stage 1 response rate to prefer B: Pr(Y1=1|T1=B,P1=1)

# Second stage treatment outcome rates
pi_AC <- 0.5                                         # Second stage response rate of non-responders to randomized A who receive randomized C in the second stage: Pr(Y2=1|T1=A,P1=0,NR,P2=0,T2=C)        
pi_AD <- 0.4                                         # Second stage response rate of non-responders to randomized A who receive randomized D in the second stage: Pr(Y2=1|T1=A,P1=0,NR,P2=0,T2=D)
pi_BC <- 0.3                                         # Second stage response rate of non-responders to randomized B who receive randomized C in the second stage: Pr(Y2=1|T1=B,P1=0,NR,P2=0,T2=C)
pi_BD <- 0.2                                         # Second stage response rate of non-responders to randomized B who receive randomized D in the second stage: Pr(Y2=1|T1=B,P1=0,NR,P2=0,T2=D)
pA0A <- pi_A * beta1_link[1]                         # Second stage response rate of responders to randomized A
pB0B <- pi_B * beta1_link[2]                         # Second stage response rate of responders to randomized B
pA1A <- pi_A * alphaP_link[1] * beta1_link[1]        # Second stage response rate of responders to preferred A
pB1B <- pi_B * alphaP_link[2] * beta1_link[2]        # Second stage response rate of responders to preferred B
pA0C1 <- alphaP_link[3] * pi_AC                      # Second stage response rate of non-responders to randomized A who receive preferred C in the second stage
pA0D1 <- alphaP_link[4] * pi_AD                      # Second stage response rate of non-responders to randomized A who receive preferred D in the second stage
pA1C0 <- alphaP_link[1] * pi_AC                      # Second stage response rate of non-responders to preferred A who receive randomized C in the second stage
pA1D0 <- alphaP_link[1] * pi_AD                      # Second stage response rate of non-responders to preferred A who receive randomized D in the second stage
pA1C1 <- alphaP_link[1] * alphaP_link[3] * pi_AC     # Second stage response rate of non-responders to preferred A who receive preferred C in the second stage
pA1D1 <- alphaP_link[1] * alphaP_link[4] * pi_AD     # Second stage response rate of non-responders to preferred A who receive preferred D in the second stage
pB0C1 <- alphaP_link[3] * pi_BC                      # Second stage response rate of non-responders to randomized B who receive preferred C in the second stage
pB0D1 <- alphaP_link[4] * pi_BD                      # Second stage response rate of non-responders to randomized B who receive preferred D in the second stage
pB1C0 <- alphaP_link[2] * pi_BC                      # Second stage response rate of non-responders to preferred B who receive randomized C in the second stage
pB1D0 <- alphaP_link[2] * pi_BD                      # Second stage response rate of non-responders to preferred B who receive randomized D in the second stage
pB1C1 <- alphaP_link[2] * alphaP_link[3] * pi_BC     # Second stage response rate of non-responders to preferred B who receive preferred C in the second stage
pB1D1 <- alphaP_link[2] * alphaP_link[4] * pi_BD     # Second stage response rate of non-responders to preferred B who receive preferred D in the second stage
