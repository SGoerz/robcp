context("CUSUM")
require(MASS)

test_that("output has the correct format", 
{
  x <- 1:5
  y <- CUSUM(x)
  
  X <- matrix(1:9, ncol = 3)
  Y <- CUSUM(X)
  
  expect_true(is.numeric(y))
  expect_true(is.numeric(Y))
  expect_equal(length(y), 1)
  expect_equal(length(Y), 1)
  expect_error(CUSUM(x, b_n = 0))
})

test_that("CUSUM test statistic is computed correctly", 
{
  x <- 1:5
  n <- length(x)
  y <- psi(x)

  CUSUM2 <- function(X)
  {
    N <- length(X)
    numerator <- cumsum(X) - (1:N) * mean(X)
    return(numerator)
  }
  
  expect_equal(CUSUM(x), max(abs(CUSUM2(x))) / sqrt(lrv(x) * n))
  expect_equal(CUSUM(y), max(abs(CUSUM2(y))) / sqrt(lrv(y) * n))
  
  m <- 3
  X <- matrix(rnorm(9), ncol = m)
  Y <- psi(X)
  sigma <- lrv(Y)
  
  teststat <- apply(Y, 2, CUSUM2)
  
  mchol <- modifChol(sigma)
  swaps <- attr(mchol, "swaps")
  mchol.inv <- chol2inv(modifChol(sigma))
  ## swap-function
  sapply(m:1, function(i)
  {
    if(i != swaps[i] + 1)
    {
      ## rows
      temp <- mchol.inv[i, ]
      mchol.inv[i, ] <<- mchol.inv[swaps[i] + 1, ]
      mchol.inv[swaps[i] + 1, ] <<- temp
      ## columns
      temp <- mchol.inv[, i]
      mchol.inv[, i] <<- mchol.inv[, swaps[i] + 1]
      mchol.inv[, swaps[i] + 1] <<- temp
    }
  })
  
  g.inv <- ginv(sigma)
  
  res1 <- max(apply(teststat, 1, function(x) t(x) %*% mchol.inv %*% x)) / nrow(Y)
  res2 <- max(apply(teststat, 1, function(x) t(x) %*% g.inv %*% x)) / nrow(Y)
  
  expect_equal(res1, CUSUM(Y, inverse = "Cholesky"), tolerance = 1e-5)
  expect_equal(res2, CUSUM(Y, inverse = "generalized"), tolerance = 1e-5)
})

