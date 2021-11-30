context("Hodges_Lehmann")

test_that("Output of scale_stat has the correct format", 
{
  x <- 1:5
  y <- shift_stat(x)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(shift_stat(x, version = "a"))
  expect_error(shift_stat(1))
})


test_that("scale_stat is computed correctly", 
{
  x <- c(92, 49, 48, 57, 27)
  b_n <- 2
  
  # empVar:
  y <-  5.44 / sqrt(5) / sqrt(lrv(x, "kernel", control = list(b_n = b_n)))
  z <- scale_stat(x, version = "empVar", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # MD:
  y <-  5.44 / sqrt(5) / sqrt(lrv(x, "kernel", control = list(b_n = b_n)))
  z <- scale_stat(x, version = "MD", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # GMD:
  y <-  5.44 / sqrt(5) / sqrt(lrv(x, "kernel", control = list(b_n = b_n)))
  z <- scale_stat(x, version = "GMD", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # MAD:
  y <-  5.44 / sqrt(5) / sqrt(lrv(x, "kernel", control = list(b_n = b_n)))
  z <- scale_stat(x, version = "MAD", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # QBeta:
  y <-  5.44 / sqrt(5) / sqrt(lrv(x, "kernel", control = list(b_n = b_n)))
  z <- scale_stat(x, version = "QBeta", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  x <- c(58, 2, 59, 26, 20, 88)
  b <- 4
  l <- 3
  y <- sqrt(6) * 310 / 36 * u_hat(x - c(rep(0, 5), 62), b) /
    sqrt(lrv(x, "subs", control = list(l = l, overlapping = TRUE, distr = TRUE)))
  z <- HodgesLehmann(x, b_u = b, control = list(l = l, overlapping = TRUE, distr = TRUE))
  attributes(z) <- NULL
  
  expect_equal(z, y)
})


test_that("Output of scale_cusum has the correct format", 
{
  
})


test_that("CUSUM test for changes in the scale is performed correctly", 
{
  
})