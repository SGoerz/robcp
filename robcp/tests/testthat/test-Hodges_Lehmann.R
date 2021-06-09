context("Hodges_Lehmann")

test_that("output has the correct format", 
{
  x <- 1:5
  y <- HodgesLehmann(x, b = 1)
  
  # X <- matrix(1:9, ncol = 3)
  # Y <- HodgesLehmann(X)
  
  expect_true(is.numeric(y))
  # expect_true(is.numeric(Y))
  expect_equal(length(y), 1)
  # expect_equal(length(Y), 1)
  expect_error(HodgesLehmann(x, b = 0))
})

test_that("u_hat computes the correct value",
{
  x <- rnorm(5)
  expect_error(u_hat(x, -2, kFun = "FT"))
  
  b <- 2
  y <- 1 - abs(apply(combn(x, 2), 2, diff) / b)
  y <- sum(ifelse(y < 0, 0, y))
  expect_equal(2 / (5 * 4 * b) * y, u_hat(x, b))
})

test_that("HodgesLehmann computes the correct value",
{
  x <- c(14, 49, 50, 47, 28)
  b <- 3
  l <- 2
  y <- sqrt(5) * 5.44 * u_hat(x - c(0, rep(34, 4)), b) / 
    sqrt(lrv(x, "subsampling", control = list(l = l, overlapping = TRUE, distr = TRUE)))
  z <- HodgesLehmann(x, control = list(l = l, overlapping = TRUE, distr = TRUE),
                     b = b)
  attributes(z) <- NULL
  
  expect_equal(z, y)
  
  x <- c(58, 2, 59, 26, 20, 88)
  b <- 4
  l <- 3
  y <- sqrt(6) * 310 / 36 * u_hat(x - c(rep(0, 5), 62), b) /
    sqrt(lrv(x, "subs", control = list(l = l, overlapping = TRUE, distr = TRUE)))
  z <- HodgesLehmann(x, control = list(l = l, overlapping = TRUE, distr = TRUE),
                                       b = b)
  attributes(z) <- NULL
  
  expect_equal(z, y)
})