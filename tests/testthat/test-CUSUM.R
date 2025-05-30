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
  expect_error(CUSUM(x, control = list(b_n = 0)))
})

test_that("CUSUM test statistic is computed correctly", 
{
  require(MASS)
  
  x <- 1:5
  n <- length(x)
  y <- psi(x)

  CUSUM2 <- function(X)
  {
    N <- length(X)
    numerator <- cumsum(X) - (1:N) * mean(X)
    return(numerator)
  }
  
  cx <- CUSUM(x, control = list(b_n = 2))
  attributes(cx) <- NULL
  cy <- CUSUM(y, control = list(b_n = 2))
  attributes(cy) <- NULL
  
  expect_equal(cx, max(abs(CUSUM2(x))) / 
                 sqrt(lrv(x, control = list(kFun = "TH", b_n = 2)) * n))
  expect_equal(cy, max(abs(CUSUM2(y))) / 
                 sqrt(lrv(y, control = list(kFun = "TH", b_n = 2)) * n))
  
  m <- 12
  X <- matrix(rnorm(m^2), ncol = 3)
  Y <- psi(X)
  sigma <- lrv(Y)
  
  teststat <- apply(Y, 2, CUSUM2)
  
  mchol <- modifChol(sigma)
  mchol.inv <- chol2inv(mchol)
  
  res <- svd(sigma)
  svd1inv <- res$u %*% diag(1 / res$d) %*% t(res$v)
  
  gen.inv <- ginv(sigma)
  
  # swaps <- attr(mchol, "swaps")
  # mchol.inv <- chol2inv(modifChol(sigma))
  # ## swap-function
  # sapply(m:1, function(i)
  # {
  #   if(i != swaps[i] + 1)
  #   {
  #     ## rows
  #     temp <- mchol.inv[i, ]
  #     mchol.inv[i, ] <<- mchol.inv[swaps[i] + 1, ]
  #     mchol.inv[swaps[i] + 1, ] <<- temp
  #     ## columns
  #     temp <- mchol.inv[, i]
  #     mchol.inv[, i] <<- mchol.inv[, swaps[i] + 1]
  #     mchol.inv[, swaps[i] + 1] <<- temp
  #   }
  # })
  # 
  # g.inv <- ginv(sigma)
  
  res1 <- max(apply(teststat, 1, function(x) t(x) %*% mchol.inv %*% x)) / nrow(Y)
  res2 <- max(apply(teststat, 1, function(x) t(x) %*% svd1inv %*% x)) / nrow(Y)
  res3 <- max(apply(teststat, 1, function(x) t(x) %*% gen.inv %*% x)) / nrow(Y)
  
  Ychol <- CUSUM(Y, inverse = "Cholesky")
  attributes(Ychol) <- NULL
  Ysvd <- CUSUM(Y, inverse = "svd")
  attributes(Ysvd) <- NULL
  Ygen <- CUSUM(Y, inverse = "generalized")
  attributes(Ygen) <- NULL
  
  expect_equal(res1, Ychol, tolerance = 1e-5)
  expect_equal(res2, Ysvd, tolerance = 1e-5)
  expect_equal(res3, Ygen, tolerance = 1e-5)
  
  # correct change point location
  x <- rnorm(100)
  x[50:100] <- x[50:100] + 10
  expect_equal(attr(CUSUM(x), "cp-location"), 50, tolerance = 1)
})

