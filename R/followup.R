#' Compute the followup priorities for a list of contacts
#'
#' Outputs the probability that the symptoms onset of a contact was between the
#' last follow up (day of last follow up not included in interval)
#' and the analysis date, which defaults to
#' day before current system date and is included in interval.
#'
#' @param contact_list A data.frame with one row per contact, containing at
#'   least a column with possible dates of exposure.
#'
#' @param incubation_period The incubation period distribution. Can be a
#'   distcrete distribution, an empirical incubation period returned by
#'   epitrix::empirical_incubation_dist() (of type data.frame) or an atomic
#'   vector whose elements 1:n correspond to the probability of the incubation
#'   period being 0:(n-1).
#'
#' @param exposure The name of the column of contact_list containing the dates
#'   of exposure (variable name or in quotes). Can be a list column
#'   containing vectors with several possible exposure dates per contact.
#'
#' @param exposure_end The name of a column containing dates representing the
#'   end of the exposure period. This is `NULL` by default, indicating all
#'   exposures are in the `exposure` column.
#'
#' @param rate_infectious_contact a vector containing rates of infectious contact corresponding
#'   to each exposure date in exposure, or between exposure and exposure_end.
#'   same length as exposure. Defaults to rep(1, length(exposure)).
#'
#' @param last_followup The name of the column of contact_list containing the
#'   last follow up date for each contact (variable name or in quotes).
#'   Should contain NA if unknown or contact has never been followed.
#' 
#' @param date_analysis the date on which the prioritization should be done. The
#'   probability of symptoms onset is computed until the end of the previous
#'   day. Defaults to the system date (i.e. today).
#'
#' @param sort If TRUE (default) the result is sorted by followup priority.
#'
#' @return A data.frame containing all columns in contact_list and the following
#'   additional columns:
#' - rate_onset the overall rate of disease onset by
#'   date_analysis
#' - p_onset the overall probability of disease onset by
#'   date_analysis
#' - followup_priority, an index ranging from 1 for the highest
#'   priority to the number of patients
#'
#' @export
#' @importFrom rlang enquo "!!" get_expr
#' @importFrom dplyr pull if_else
#' @importFrom tidyr complete full_seq
#' @importFrom purrr map2_dbl pmap_dbl
followup_priorities <- function(contact_list, exposure, exposure_end = NULL,
                                rate_infectious_contact = NULL,
                                last_followup = NULL,
                                incubation_period = NULL,
                                date_analysis = Sys.Date(),
                                sort = TRUE) {
  #-------------------------------------------------------------
  #------------- check inputs and transform --------------------
  #-------------------------------------------------------------

  #is contact list a df and has at least one row and col
  stopifnot(is.data.frame(contact_list), nrow(contact_list) > 0, ncol(contact_list) > 0)

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

    if (any(last_followup > date_analysis, na.rm = TRUE)) {
      warning("Some followup dates are after the analysis date. Ignoring them.")
      last_followup[last_followup > date_analysis] <- as.Date(NA)
    }
  } else {
    stop("last_followup is not a column of contact_list")
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

  #### THIS PART COULD BE DE-TIDYVERSIFIED IF NEEDED

  #prepare df to work with
  work_df <- tibble::tibble(
    exposure,
    last_followup,
    rate_infectious_contact
  )

  work_df <- tibble::rowid_to_column(work_df)

  #one row per exposure date
  work_df <- tidyr::unnest(work_df, c(exposure, rate_infectious_contact))
  
  #compute rates
  inc <- integrate_incubation_followup(
    incubation_period = incubation_period,
    date_exposure = work_df$exposure,
    date_analysis = date_analysis,
    date_last_followup = work_df$last_followup
    )

  work_df$rate_onset <- work_df$rate_infectious_contact * inc
  
  #sum rates of several for same contact exposures and put back to original df
  work_df <- summarize(group_by(work_df, rowid), rate_onset = sum(rate_onset), .groups = "drop")
  contact_list$rate_onset <- work_df$rate_onset

  #compute probability of onset from rate
  contact_list$p_onset <- 1 - exp(-contact_list$rate_onset)

  #sort and add priorities
  contact_list$followup_priority[order(contact_list$p_onset, decreasing = TRUE)] <- seq_len(nrow(contact_list))
  if (sort) {
    contact_list <- contact_list[order(contact_list$followup_priority),]
  }

  return(contact_list)
}
