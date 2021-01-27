#' @param incubation_period a vector of probabilities. Has to sum to 1.
#' @param date_low the lower date (integer or Date object)
#' @param date_high the upper date (integer or Date object)
#' 
#' @export
integrate_incubation <- function(incubation_period, date_low, date_high) {

    if (!any(
            inherits(date_low, "Date"), #check if date object
            (is.numeric(date_low) & as.numeric(date_low)%%1 == 0) #check if integer
        ) |
        length(date_low) != 1
    ) {
        stop("date_low needs to be integer or Date objecs of length 1")
    }

    if (!any(
            inherits(date_high, "Date"), #check if date object
            (is.numeric(date_high) & as.numeric(date_high)%%1 == 0) #check if integer
        ) |
        length(date_high) != 1
    ) {
        stop("date_high needs to be integer or Date object of length 1")
    }

    if (date_high < date_low) {
        stop("date_high needs to be greater or equal than date_low.")
    }

    if (inherits(date_high, "Date") & !inherits(date_low, "Date")) {
        stop("dates must be the same class (numeric or Dates)")
    }

    idx <- seq(0, as.integer(date_high - date_low)) + 1
    idx <- idx[idx > 0 & idx <= length(incubation_period)] #remove index larger than maximal incubation period and smaller than 1

    if (length(idx) == 0) {
        s_incubation <- 0 #all exposures are further in the past than the maximal incubation period or after the analysis date
    } else {
        s_incubation <- sum(incubation_period[idx])
    }

    return(s_incubation)
}



#' @param incubation_period a vector of probabilities. Has to sum to 1.
#' @param date_exposure the exposure date (integer or Date object)
#' @param date_analysis the analysis date (integer or Date object)
#' @param date_last_followup date of last follow-up  (integer or Date object, NA if not yet followed-up)
#' 
#' @export
integrate_incubation_followup <- function(incubation_period, date_exposure, date_analysis, date_last_followup = NA) {

    if (!any(
            inherits(date_last_followup, "Date"), #check if date object
            (is.numeric(date_last_followup) & as.numeric(date_last_followup)%%1 == 0), #check if integer
            is.na(date_last_followup)
        ) |
        (length(date_last_followup) != 1)
    ) {
        stop("date_last_followup needs to be integer or Date objecs of length 1")
    }

    if (!is.na(date_last_followup) & inherits(date_last_followup, "Date") & !inherits(date_exposure, "Date")) {
        stop("dates must be the same class (numeric or Dates)")
    }

    if (!is.na(date_last_followup) & date_last_followup > date_analysis) {
        stop("analysis date is before follow-up date")
    }

    if (is.na(date_last_followup) | (date_last_followup < date_exposure)) {
        s_incubation <- integrate_incubation(incubation_period, date_exposure, date_analysis)
    } else {
        up_to_follow_up <- integrate_incubation(incubation_period, date_exposure, date_last_followup)
        up_to_analysis <- integrate_incubation(incubation_period, date_exposure, date_analysis)
        
        if (up_to_follow_up >= 1) {
            s_incubation <- 0 #in this case the last follow-up was after the possible incubation period and showed that no infection developped
            return(s_incubation)
        }

        num <- up_to_analysis - up_to_follow_up
        denom <- 1 - up_to_follow_up
        s_incubation <- num / denom
    }

    return(s_incubation)
}



#' @param incubation_period a vector of probabilities. Has to sum to 1.
#' @param date_exposure the exposure date (integer or Date object). can be a vector.
#' @param date_analysis the analysis date (integer or Date object). can be a vector.
#' @param date_last_followup date of last follow-up  (integer or Date object, NA if not yet followed-up).  can be a vector.
#' 
#' If any of the dates is a vector the output is a vector of the same length
#' 
#' @export
integrate_incubation_followup <- Vectorize(integrate_incubation_followup, vectorize.args = c("date_exposure", "date_analysis", "date_last_followup"))
