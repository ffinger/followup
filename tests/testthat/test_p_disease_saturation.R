context("Testing p_disease_saturation function")

test_that("an error is thrown if p_disease_max is not between 0 and 1", {
  expect_error(
    p_disease_saturation(d = 5, p_disease_max = -0.1, d50 = 1),
    "p_disease_max has to be between 0 and 1 and of length 1."
  )
  expect_error(
    p_disease_saturation(d = 5, p_disease_max = 1.1, d50 = 1),
    "p_disease_max has to be between 0 and 1 and of length 1."
  )
})

test_that("an error is thrown if d is below 0", {
  expect_error(
    p_disease_saturation(d = -1, p_disease_max = 1, d50 = 1),
    "d has to be greater than or equal to 0."
  )
})

test_that("an error is thrown if d50 is below 0", {
  expect_error(
    p_disease_saturation(d = 3, p_disease_max = 1, d50 = -1),
    "d50 has to be greater than or equal to 0 and of length 1."
  )
})

test_that("p_disease is 0.5*p_disease_max if d=d50", {
  expect_identical(
    p_disease_saturation(d = 1, p_disease_max = .4, d50 = 1),
    .2
  )
})

test_that("p_disease is correct for test value", {
  expect_identical(
    p_disease_saturation(d = 0:5, p_disease_max = .4, d50 = 1),
    .4*(1-exp(-log(2)*0:5))
  )
})
