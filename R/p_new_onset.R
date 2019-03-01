#' Computes the probability of new symptoms on a particular day for a particular contact
#'
#' @param incubation_period a vector of probabilities. Has to sum to 1.
#' @param date_analysis the date for which to compute the probability of new symptoms
#' @param exposure a vector containing one or several possible dates of exposure. They are equally weighted.
#' @return a probability
#'
#' @export
p_new_onset <- function(incubation_period, date_analysis, exposure) {

  if (sum(incubation_period) != 1) {
    stop("incubation_period doesn't sum to 1.")
  }

  if (any(incubation_period < 0)) {
    stop("incubation_period contains negative elements.")
  }

  idx <- as.integer(date_analysis - exposure) + 1
  incubation_period <- incubation_period[ 1 : max( which( incubation_period != 0 )) ] #remove trailing 0s
  idx <- idx[idx <= length(incubation_period)] #remove index larger than maximal incubation period

  # if (any(idx < 1)) {
  #   warning("Some exposure dates are after analysis dates.")
  # }

  idx <- idx[idx > 0] #remove index for exposure dates after analysis dates

  if (length(idx) == 0) {
    p <- 0 #all exposures are further in the past than the maximal incubation period or after the analysis date
  } else {
    p <- sum(incubation_period[idx])/length(exposure)
  }

  return(p)
}

