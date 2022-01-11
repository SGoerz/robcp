context("rho test")

test_that("Output of rho_stat has the correct format", 
{
  x <- matrix(1:12, ncol = 2)
  y <- rho_stat(x)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(rho_stat(1))
})

test_that("rho_stat is computed correctly", 
{
  x <- matrix(c(80, 23, 74, 93, 39, 48, 26, 77), ncol = 2)
  b_n <- 2
  
  y <- 33/16 / sqrt(4 * lrv(x, "kernel", control = list(b_n = b_n, version = "rho")))
  z <- rho_stat(x, control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
})

test_that("Output of rho_cusum has the correct format", 
{
  x <- matrix(rnorm(10), ncol = 2)
  res <- suppressWarnings(rho_cusum(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "CUSUM test for changes in Spearman's rho")
})

test_that("CUSUM test for changes in the scale is performed correctly", 
{
  require(mvtnorm)
  
  n <- 500
  m <- 100
  N <- n + m
  
  rho <- 0.4
  
  theta1 <- 0
  theta2 <- 0
  theta <- cbind(c(theta1, 0), c(0, theta2))
  q <- rho * sqrt( (theta1^2 + 1) * (theta2^2 + 1) / (theta1 * theta2 + 1))
  S <- cbind(c(1, q), c(q, 1))
  
  y <- replicate(1000, 
  {
    e <- rmvt(N, S)
    x <- matrix(numeric(N * 2), ncol = 2)
    x[1, ] <- e[1, ]
    t(sapply(2:N, function(i) x[i, ] <<- e[i, ] + theta %*% e[i-1, ]))
    
    x <- x[-(1:m), ]
    rho_cusum(x)$p.value
  })
  
  hist(y)
  mean(y < 0.05)

  
  # correct change point location
  x1 <- rnorm(100)
  x2 <- arima.sim(list(ar = 0.8), 100)
  x <- matrix(c(x1, x2), byrow = TRUE, ncol = 2)
  expect_equal(attr(rho_cusum(x)$statistic, "cp-location"), 50, tolerance = 1)
  
  ## maybe some more tests
  ## best to be checked graphically:
  ## hist(p)
})