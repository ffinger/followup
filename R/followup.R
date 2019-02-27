#' Compute the followup priorities for a list of contacts
#'
#' Outputs the probability that the symptoms onset of a contact was between the last follow up (day of last follow up included in interval by default, see parameter include_last_follow_up) and the analysis date, which defaults to day before current system date and is included in interval.
#'
#' @param contact_list A data.frame with one row per contact, containing at least a list column with possible dates of exposure.
#' @param incubation_period The incubation period distribution. Can be a distcrete distribution, an empirical incubation period returned by epitrix::empirical_incubation_dist() (of type data.frame) or an atomic vector whose elements 1:n correspond to the probability of the incubation period being 0:(n-1).
#' @param exposure The name of the column of contact_list containing the dates of exposure (bare variable name or in quotes). Can be a list column containing vectors with several possible exposure dates per contact.
#' @param exposure_end the name of a column containing dates representing the
#'   end of the exposure period. This is `NULL` by default, indicating
#'   all exposures are known and in the `exposure` column.
#' @param last_followup The name of the column of contact_list containing the last follow up date for each contact (bare variable name or in quotes). Should contain NA if unknown or contact has never been followed.
#' @param p_disease The overall probability that exposure leads to disease. Can either be a scalar applying to all contacts, or the name of a column (bare variable name or in quotes) in contact_list containing a different probability per contact. Defaults to 1.
#' @param date_analysis the date on which the prioritization should be done. The probability of symptoms onset is computed until the end of the previous day. Defaults to the system date (i.e. today).
#' @param include_last_follow_up TRUE if the date of last follow up should be included in the date range and thus the probability of symptoms starting on that date included in the result, FALSE if not (default TRUE).
#' @param sort If TRUE (default) the result is sorted by followup priority.
#'
#' @return A data.frame containing all columns in contact_list and the following additional columns:
#' @return   * p_onset, the probability of disease has onset by date_analysis given that the exposure results in disease
#' @return   * p_symptoms, the overall probability of disease has onset by date_analysis
#' @return   * followup_priority, an index ranging from 1 for the highest priority to the number of patients
#'
#' @export
#' @importFrom rlang enquo "!!" get_expr
#' @importFrom dplyr pull if_else
#' @importFrom tidyr complete full_seq
#' @importFrom purrr map2_dbl pmap_dbl
#' @importFrom checkmate assert_data_frame
followup_priorities <- function(contact_list, exposure, exposure_end = NULL,
  last_followup = NULL, p_disease = 1, incubation_period = NULL,
  date_analysis = Sys.Date(), include_last_follow_up = TRUE, sort = TRUE) {

  #-------------------------------------------------------------
  #------------- check inputs and transform --------------------
  #-------------------------------------------------------------

  #is contact list a df and has at least one row and col
  checkmate::assert_data_frame(contact_list, min.rows = 1, min.cols = 1)

  #get dates exposure
  exposure <- rlang::enquo(exposure)
  exposure <- dplyr::pull(contact_list, !!exposure)
  exposure_end  <- rlang::enquo(exposure_end)
  end_is_here   <- !is.null(rlang::get_expr(exposure_end))

  if (end_is_here) {
    # We need to create the list for each date
    if (is.list(exposure) || !inherits(exposure, "Date")) {
      stop("if exposure_end is specified, then exposure must be a vector of Dates")
    }
    e     <- exposure
    ee    <- dplyr::pull(contact_list, !! exposure_end)
    exposure <- vector(mode = "list", length = length(e))
    for (i in seq(exposure)) {
      exposure[[i]] <- seq(from = e[i], to = ee[i], by = "1 day")
    }
  }

  #get last_followup
  last_followup <- rlang::enquo(last_followup)

  if (is.null(rlang::get_expr(last_followup))) {
    #is last followup not given --> put NA.
    last_followup <- rep(as.Date(NA), nrow(contact_list))

  } else if (rlang::quo_text(last_followup) %in% names(contact_list)) {
    #is a quosure

    last_followup <- dplyr::pull(contact_list, !!last_followup)

    if (any(last_followup > date_analysis, na.rm = TRUE) & include_last_follow_up) {
      warning("Some followup dates are after the analysis date. Ignoring them.")
      last_followup[last_followup > date_analysis] <- as.Date(NA)
    } else if (any(last_followup >= date_analysis, na.rm = TRUE) & !include_last_follow_up) {
      warning("Some followup dates are equal or after the analysis date. Ignoring them.")
      last_followup[last_followup >= date_analysis] <- as.Date(NA)
    }
  } else {
    stop("last_followup is not a column of contact_list")
  }

  #get p_disease
  if (rlang::quo_text(rlang::enquo(p_disease)) %in% names(contact_list)) { #is a column
    p_disease <- rlang::enquo(p_disease)
    p_disease <- dplyr::pull(contact_list, !!p_disease)
  } else if ((inherits(p_disease, "numeric") & (length(p_disease) == 1)) ) { #is a single numeric
    contact_list$p_disease <- p_disease
  } else {
    stop("p_disease is not a scalar and not a column of contact_list")
  }

  if (any(p_disease < 0) | any(p_disease > 1, na.rm = TRUE)) {
    stop("p_disease should be between 0 and 1")
  }

  #incubation_period
  if (is.null(incubation_period)) {
    stop("Required argument incubation_period not given.")
  }

  # make a simple vector for incubation period distribution
  if (inherits(incubation_period, "distcrete")) {
    #from distcrete

    max_inc <- as.integer(
      max(date_analysis - as.Date(unlist(exposure)), origin = as.Date("1970-01-01"))
    )

    incubation_period <- incubation_period$d(0:max_inc)

  } else if (inherits(incubation_period, "data.frame")) {
    #from df given by epitrix::empirical_incubation_dist()

    if (any(incubation_period$incubation_period < 0)) {
      stop("Incubation periods can't be negative")
    }

    incubation_period <- tidyr::complete(incubation_period,
        incubation_period = tidyr::full_seq(c(0, incubation_period), 1),
        fill = list(relative_frequency = 0)
    )

    incubation_period <- incubation_period$incubation_period
  }

  #remove trailing 0s
  incubation_period <- incubation_period[ 1 : max( which(incubation_period != 0 )) ]


  # check if incubation_period distribution sums to 1, otherwise force it
  if (sum(incubation_period) != 1) {
    warning(paste0("Incubation period probabilities don't sum to 1 but to ", sum(incubation_period), " Automatically adjusted."))
    incubation_period <- incubation_period/sum(incubation_period)
  }

  if (any(is.na(incubation_period))) {
    stop("incubation_period contains NA.")
  }

  if (date_analysis < min(unlist(exposure))) {
    stop("date_analysis before first exposure date.")
  }


  #----------------------------------------------
  #------------- do the work --------------------
  #----------------------------------------------

  #the first exposure date for each contact
  min_exp <- as.Date(vapply(exposure, min, 1), origin = as.Date("1970-01-01"))

  increment <- 0
  if (!include_last_follow_up) {
    increment <- 1
  }

  #lower lim of integration
  date_lower <- dplyr::if_else(
    is.na(last_followup),
    min_exp,
    last_followup + increment #add increment of 1 if last followup not included in integration range
  )

  #compute probs
  p_onset <- purrr::map2_dbl(
    date_lower,
    exposure,
    function(x,y) p_integral_onset(incubation_period, x, date_analysis, y)
  )

  #put into df
  contact_list$p_onset <- p_onset

  #apply base probability of disease being transmitted
  contact_list$p_symptoms <- contact_list$p_onset * p_disease

  #sort and add priorities
  contact_list$followup_priority[order(contact_list$p_symptoms, decreasing = TRUE)] <- seq_len(nrow(contact_list))
  if (sort) {
    contact_list <- contact_list[order(contact_list$followup_priority),]
  }

  return(contact_list)
}


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



#' Helper function to compute the probability p_disease of getting disease after a certain number of days of exposure based on a saturating exponential growth model.
#'
#' @param d the number of exposure days, larger than 0. Can be a vector.
#' @param p_disease_max the maximum probability of getting disease, reached for high d only.
#' @param d50 the number of exposure days after which p_disease = 0.5*p_disease_max is reached. Can be float.
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


#' Helper function print contact list as html table.
#'
#' @param contact_list a data frame
#' @param cols_round a list of quoted column names to be rounded
#' @param date_cols a list of quoted column names containing dates
#' @param ndigit number of digits for rouning (default: 2)
#'
#' @importFrom dplyr mutate
#' @importFrom htmlTable txtRound
#' @importFrom kableExtra kable kable_styling
#' @importFrom magrittr "%>%"
#'
#' @export
to_html <- function(contact_list, cols_round = list(), date_cols = list(), ndigit = 2) {

  # date_cols <- dplyr::quos(unlist(date_cols))

  if (length(cols_round) > 0) {
    for (i in seq_along(cols_round)) {
      col <- cols_round[[i]]
      contact_list[col] <- htmlTable::txtRound(contact_list[col], ndigit)
    }
  }

  if (length(date_cols) > 0) {
    for (i in seq_along(date_cols)) {
      col <- date_cols[[i]]
      if (inherits(contact_list[[col]], "Date")) {
          contact_list[col] <- lapply(contact_list[col], as.character)
      } else if (inherits(contact_list[[col]], "list")) {
          contact_list[col] <- lapply(contact_list[col], datevec_to_txt)
      } else {
        stop(paste("unexpected column type", class(contact_list[[col]])))
      }
    }
  }

  kableExtra::kable(contact_list) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)
}

datevec_to_txt <- function(datevec) {
    datevec <- lapply(datevec, as.Date)
    datevec <- lapply(datevec, as.character)
    return(paste(datevec))
}
