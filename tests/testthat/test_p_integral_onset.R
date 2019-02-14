context("Testing p_integral_onset function")

inc <- c(0.1, 0.15, 0.2, 0.25, 0.3)

test_that("p is 0 if date_lower == date_upper", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-07"), as.Date("2019-01-07"), c(as.Date("2019-01-02"), as.Date("2019-01-04"))), 0)
})

test_that("p is 0 if all exposure further in the past than incubation period", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-02-05"), c(as.Date("2019-01-02"), as.Date("2019-01-04"))), 0)
})

test_that("p is 0 if all exposures after analysis", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-02-05"), c(as.Date("2019-02-10"), as.Date("2019-03-04"))), 0)
})

test_that("p is right for single exposure single analysis date", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-01-31"), c(as.Date("2019-01-28"))), inc[3])
})

test_that("p is right for multiple exposure single analysis date", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-01-31"), c(as.Date("2019-01-26"), as.Date("2019-01-28"))), (inc[3] + inc[5])/2)
})

test_that("p is right for multiple exposure (one further back than max incubation duration) and single analysis date", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-01-31"), c(as.Date("2019-01-25"), as.Date("2019-01-28"))), inc[3])
})

test_that("p is right for multiple exposure (one after analysis) and single analysis date", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-30"), as.Date("2019-01-31"), c(as.Date("2019-02-01"), as.Date("2019-01-28"))), inc[3])
})

test_that("p is right multiple analysis dates and single exposure date previous to analysis dates", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-25"), as.Date("2019-01-28"), c(as.Date("2019-01-24"))), sum(inc[2:4]))
})

test_that("p is right multiple analysis dates and single exposure date in the middle of analysis dates", {
  expect_identical(p_integral_onset(inc, as.Date("2019-01-25"), as.Date("2019-01-29"), c(as.Date("2019-01-26"))), sum(inc[1:3]))
})

test_that("p is right multiple analysis dates and multiple analysis and multiple exposure dates", {
  expect_identical(
    p_integral_onset(
      inc,
      as.Date("2019-01-25"),
      as.Date("2019-01-29"),
      c(as.Date("2019-01-23"), as.Date("2019-01-24"), as.Date("2019-01-26"))
    ),
    sum(inc[c(2,3)])/2 + sum(inc[c(1,3,4)])/3 + sum(inc[c(2,4,5)])/3 + sum(inc[c(3,5)])/2
  )
})
