context("Testing p_new_onset function")

test_that("an error is thrown if incubation_period doesn't sum to 1", {
  expect_error(p_new_onset(c(.3, .5, .1), as.Date("2019-01-20"), as.Date("2019-01-01")), "incubation_period doesn't sum to 1.")
})

test_that("an error is thrown if incubation_period contains negatives", {
  expect_error(p_new_onset(c(.2, .6, -.2, .4), as.Date("2019-01-20"), as.Date("2019-01-01")), "incubation_period contains negative elements.")
})

inc <- c(0.1, 0.15, 0.2, 0.25, 0.3)

test_that("p is right for single exposure date", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), as.Date("2019-01-03")), 0.25)
})

test_that("p is 0 if exposure date is after analysis date", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), as.Date("2019-01-07")), 0)
})

test_that("p is 0 if exposure date is too long in the past", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), as.Date("2019-01-06") - length(inc)), 0)
})

test_that("p is right if incubation period is 0", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), as.Date("2019-01-06")), 0.1)
})

test_that("p is right for two exposure dates, one of them outside range", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), c(as.Date("2019-01-01"), as.Date("2019-01-03"))), .25/2)
})

test_that("p is right for several exposure dates, all within range", {
  expect_identical(p_new_onset(inc, as.Date("2019-01-06"), c(as.Date("2019-01-02"), as.Date("2019-01-03"), as.Date("2019-01-05"))), .7/3)
})
