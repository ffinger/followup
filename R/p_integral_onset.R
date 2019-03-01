#' Computes the probability of new symptoms within a range of dates for a particular contact, given the incubation period.
#'
#' The lower but not the upper limit is included in the interval such that the probability returned is the sum over the probability for each date_analysis that fulfils the following: date_lower <= date_analysis < date_upper.
#'
#' @param incubation_period a vector of probabilities. Has to sum to 1.
#' @param date_lower the lower limit for which to compute the probability of new symptoms.
#' @param date_upper the upper limit for which to compute the probability of new symptoms.
#' @param exposure a vector containing one or several possible dates of exposure. They are equally weighted.
#' @return a probability
#'
#' @export
p_integral_onset <- function(incubation_period, date_lower, date_upper, exposure) {

  fct <- function(x) p_new_onset(incubation_period, x, exposure)

  if (date_lower > date_upper) {
    stop("date_upper is before date_lower.")
  } else if (date_lower == date_upper) {
    p <- 0
  } else {
    p <- sum(vapply(seq(from = date_lower, to = date_upper - 1, by = 1), fct, 1))
  }

  return(p)
}
