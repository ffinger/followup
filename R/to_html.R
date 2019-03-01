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
