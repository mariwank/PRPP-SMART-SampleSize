# PRPP-SMART-SampleSize
Companion code for "Bayesian Sample Size Determination in a Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial with Binary Outcomes”. Where applicable R code has been adapted from Artman et al. (2022): Bayesian set of best dynamic treatment regimes: Construction and sample size calculation for SMARTs with binary outcomes, [https://doi.org/10.1111/insr.12376](https://doi.org/10.1002/sim.9323).

## File Descriptions
- [PRPPSMART_DataGen.R](PRPPSMART_DataGen.R): Function used to create exemplarly PRPP-SMART data with binary end-of-stage outcomes to be used in the PRPP-SMART sample size calculation.
- [ComputePower_MCMC.R](ComputePower_MCMC.R): Code to calculate required sample size for a PRPP-SMART using MCMC-based posterior distributions. Note, calculation is reccomened to be run on a high-performance cluster rather than a laptop if setting $I>100$ and not implementing a stopping rule.
- [ComputePower_stopping.R](ComputePower_MCMC.R): Code to calculate required sample size for a PRPP-SMART using MCMC-based posterior distributions and a stopping rule to end the caculation once desired power is achieved.
- [PowerFunction_Approx.R](PowerFunction_Approx.R): Function to calculate the power in a PRPP-SMART using approximate closed-form posterior distributions. Called in [ComputePower_Approx.R](ComputePower_Approx.R).
- [ComputePower_Approx.R](ComputePower_Approx.R): Code to implement [PowerFunction_Approx.R](PowerFunction_Approx.R) to determine the required sample size for a PRPP-SMART. This code implements parallel processing of the sample size calculation.
- [Artman_SampleSize_calc.R](Artman_SampleSize_calc.R): Example R code to calculate $n_{min}$ in the PRPP-SMART sample size calculation using the method outlined in Artman et al. (2022) Bayesian set of best dynamic treatment regimes: Construction and sample size calculation for SMARTs with binary outcomes, [https://doi.org/10.1111/insr.12376](https://doi.org/10.1002/sim.9323).
- [LogOR_Function.R](LogOR_Function.R): Function to calculate the log-OR between each indfference DTR and the best indifference DTR given the sample size input parameters.  
- [nsim_toget_500.R](nsim_toget_500.R): Functions used to determine the number of simulations needed to run per scenario/sub-scenario in order to achieve 500 total simulations. Used for Empirical power calculations to ensure positivity assumption of BJSM model is met. 


## Folder Descriptions
- The [Empirical](Empirical) folder contains the code used to calculate empirical power in simulated PRPP-SMART datasets using the BJSM method. 

