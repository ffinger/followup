context("Testing followup_priorities function")

myrank <- function(x) {
  res <- NA
  res[order(x, decreasing = TRUE)] <- 1:length(x)
  return(res)
}

inc <- c(0.1, 0.15, 0.2, 0.25, 0.3)

date_analysis <- as.Date("2019-02-15")

contact_list <- data.frame(
  date_last_followup = as.Date(c(NA, "2019-02-11", "2019-02-13", "2019-02-15"))
)

contact_list$dates_exposure <- list(
  as.Date(c("2019-02-11", "2019-02-12", "2019-02-13")),
  as.Date(c("2019-02-11", "2019-02-12", "2019-02-13")),
  as.Date(c("2019-02-11", "2019-02-12", "2019-02-13")),
  as.Date(c("2019-02-11", "2019-02-12", "2019-02-13"))
)

p_d <- 0.3

contact_list_2 <- contact_list
contact_list_2$p_disease <- c(0.2, 0.3, 0.1, 0.4)


expected_p_onset <- c(
  (0.1/3 + (0.1+0.15)/3 + (0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),
  (0.1/3 + (0.1+0.15)/3 + (0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),
  ((0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),
  0
)

expected_p_symptoms <- expected_p_onset * p_d
expected_p_symptoms_2 <- expected_p_onset * contact_list_2$p_disease

expected_p_onset_3 <- rep((0.1/3 + (0.1+0.15)/3 + (0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),4)
expected_p_symptoms_3 <- expected_p_onset_3 * contact_list_2$p_disease

expected_p_onset_4 <- c(
  (0.1/3 + (0.1+0.15)/3 + (0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),
  ((0.1+0.15)/3 + (0.1+0.15+0.2)/3 + (0.15+0.2+0.25)/3),
  ((0.15+0.2+0.25)/3)
)
expected_p_symptoms_4 <- expected_p_onset_4 * contact_list_2$p_disease[1:3]

test_that("output is right for constant p_disease", {
  fp <- followup_priorities(contact_list, dates_exposure, last_followup = date_last_followup, p_disease = p_d, incubation_period = inc, date_analysis = date_analysis, include_last_follow_up = TRUE, sort = FALSE)

  expect_equal(fp$p_onset, expected_p_onset)
  expect_equal(fp$p_symptoms, expected_p_symptoms)
  expect_equal(fp$followup_priority, myrank(expected_p_symptoms))
})

test_that("output is right for varying p_disease", {
  fp <- followup_priorities(contact_list_2, dates_exposure, last_followup = date_last_followup, p_disease = p_disease, incubation_period = inc, date_analysis = date_analysis, include_last_follow_up = TRUE, sort = FALSE)

  expect_equal(fp$p_onset, expected_p_onset)
  expect_equal(fp$p_symptoms, expected_p_symptoms_2)
  expect_equal(fp$followup_priority, myrank(expected_p_symptoms_2))
})

test_that("output is right for varying p_disease with sort = TRUE", {
  fp <- followup_priorities(contact_list_2, dates_exposure, last_followup = date_last_followup, p_disease = p_disease, incubation_period = inc, date_analysis = date_analysis, include_last_follow_up = TRUE, sort = TRUE)

  expect_equal(fp$p_onset, expected_p_onset[order(expected_p_symptoms_2, decreasing = TRUE)])
  expect_equal(fp$p_symptoms, expected_p_symptoms_2[order(expected_p_symptoms_2, decreasing = TRUE)])
  expect_equal(fp$followup_priority, 1:4)
})

test_that("output is right for varying p_disease with date last follow up null", {
  fp <- followup_priorities(contact_list_2, dates_exposure, last_followup = NULL, p_disease = p_disease, incubation_period = inc, date_analysis = date_analysis, include_last_follow_up = TRUE, sort = FALSE)

  expect_equal(fp$p_onset, expected_p_onset_3)
  expect_equal(fp$p_symptoms, expected_p_symptoms_3)
  expect_equal(fp$followup_priority, myrank(expected_p_symptoms_3))
})


test_that("output is right for varying p_disease with date last follow up null and include_last_follow_up = FALSE", {
  fp <- followup_priorities(contact_list_2[1:3,], dates_exposure, last_followup = date_last_followup, p_disease = p_disease, incubation_period = inc, date_analysis = date_analysis, include_last_follow_up = FALSE, sort = FALSE)

  expect_equal(fp$p_onset, expected_p_onset_4)
  expect_equal(fp$p_symptoms, expected_p_symptoms_4)
  expect_equal(fp$followup_priority, myrank(expected_p_symptoms_4))
})


test_that("error if analysis date before first exposure date past", {
  expect_error(
    followup_priorities(contact_list_2, dates_exposure, last_followup = NULL, p_disease = p_disease, incubation_period = inc, date_analysis = as.Date("2019-01-01"), include_last_follow_up = TRUE, sort = FALSE),
    "date_analysis before first exposure date."
  )
})

test_that("no error if analysis date >= first exposure date and before all followup dates", {
    expect_warning(
      fp <- followup_priorities(contact_list_2, dates_exposure, last_followup = date_last_followup, p_disease = p_disease, incubation_period = inc, date_analysis = as.Date("2019-02-11"), include_last_follow_up = TRUE, sort = FALSE),
      "Some followup dates are after the analysis date. Ignoring them."
    )

    expect_equal(fp$p_onset, rep(0,4))
    expect_equal(fp$p_symptoms, rep(0,4))
})

test_that("output is right for varying p_disease and analysis date far in future", {
  fp <- followup_priorities(contact_list_2, dates_exposure, last_followup = NULL, p_disease = p_disease, incubation_period = inc, date_analysis = as.Date("2019-04-01"), include_last_follow_up = TRUE, sort = FALSE)

  expect_equal(fp$p_onset, rep(1,4))
  expect_equal(fp$p_symptoms, c(0.2, 0.3, 0.1, 0.4))
})


#TODO
#testing argument types and presence
#include_last_follow_up = FALSE
