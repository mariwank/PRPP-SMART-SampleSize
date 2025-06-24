LogOR <- function(response_prob,
                  stage_one_trt_one_response_prob,
                  stage_one_trt_two_response_prob) {
  
  # Arguments:
  # response_prob: probability of stage 2 outcom for each of 6 embedded treatment sequences in the following order (AA, BB, AC, AD, BC, BD)
  # stage_one_trt_one_response_prob: probability of stage 1 outcome to stage 1 treatment one
  # stage_one_trt_two_response_prob: probability of stage 1 outcome to stage 1 treatment two

  
  #Compute mean embedded dynamic treatment regime outcomes
  iDTRs <- c(
    (response_prob[1]) * stage_one_trt_one_response_prob + response_prob[3] * (1 - stage_one_trt_one_response_prob), # AAC00
    (response_prob[1]) * stage_one_trt_one_response_prob + response_prob[4] * (1 - stage_one_trt_one_response_prob), # AAD00
    (response_prob[2]) * stage_one_trt_two_response_prob + response_prob[5] * (1 - stage_one_trt_two_response_prob), # BBC00
    (response_prob[2]) * stage_one_trt_two_response_prob + response_prob[6] * (1 - stage_one_trt_two_response_prob)  # BBD00
  )
  
  
  
  
  
  # Compute log-OR
  dtrdraws_log_odds <- log(iDTRs / (1 - iDTRs))
  
  # Compute index of best EDTR
  max_odds_ind <- (which.max((dtrdraws_log_odds)))
  
  # Compute log-odds ratios between each EDTR and best
  Log_OR_output <- matrix((dtrdraws_log_odds - dtrdraws_log_odds[max_odds_ind]), nrow = 1, ncol = length(dtrdraws_log_odds))
  
  
  
  colnames(Log_OR_output) <- c(
    "iDTR 1",
    "iDTR 2",
    "iDTR 3",
    "iDTR 4"
  )
  
  return(Log_OR_output)
}
