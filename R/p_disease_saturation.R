#' Helper function to compute the probability p_disease of getting disease after
#' a certain number of days of exposure based on a saturating exponential growth
#' model.
#'
#' @param d the number of exposure days, larger than 0. Can be a vector.
#' @param p_disease_max the maximum probability of getting disease, reached for
#'   high d only.
#' @param d50 the number of exposure days after which p_disease =
#'   0.5*p_disease_max is reached. Can be float.
#'
#' @return the proabbility p_disease, same length as d.
#'
#' @export
p_disease_saturation <- function(d, p_disease_max = 1, d50 = 1) {

  if ((p_disease_max < 0) | (p_disease_max > 1) | (length(p_disease_max) > 1)) {
    stop("p_disease_max has to be between 0 and 1 and of length 1.")
  }

  if (any(d < 0)) {
    stop("d has to be greater than or equal to 0.")
  }

  if (d50 < 0 | (length(d50) > 1)) {
    stop("d50 has to be greater than or equal to 0 and of length 1.")
  }

  k <- log(2)/d50
  return( p_disease_max * (1-exp(-k*d)) )
}
