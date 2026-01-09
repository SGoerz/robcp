 context("cor test")

test_that("Output of cor_stat has the correct format", 
{
  x <- matrix(1:12, ncol = 2)
  y <- cor_stat(x)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  y <- cor_stat(x, "tau")
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(cor_stat(1))
  expect_error(cor_stat(1, "tau"))
})

test_that("cor_stat is computed correctly", 
{
  ## rho
  x <- matrix(c(80, 23, 74, 93, 39, 48, 26, 77), ncol = 2)
  b_n <- 2
  
  y <- 33/16 / sqrt(4 * lrv(x, "kernel", control = list(b_n = b_n, version = "rho")))
  z <- cor_stat(x, "rho", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  ## tau
  x <- matrix(c(82, 43, 71, 54, 90, 40, 12, 75, 91, 49), ncol = 2)
  b_n <- 2
  
  y <- 2 / sqrt(5 * lrv(x, "kernel", control = list(b_n = b_n, version = "tau",
                                                    kFun = "quadratic")))
  z <- cor_stat(x, "tau", control = list(b_n = b_n, kFun = "quadratic"))
  attributes(z) <- NULL
  expect_equal(z, y)
})

test_that("Output of cor_cusum has the correct format", 
{
  x <- matrix(rnorm(10), ncol = 2)
  res <- suppressWarnings(cor_cusum(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "CUSUM test for changes in the correlation")
  
  x <- matrix(rnorm(100), ncol = 2)
  y <- cor_cusum(x)
  
  expect_true(is.list(y))
  expect_true(is(y$lrv, "list"))
  expect_equal(class(y), "htest")
  
  expect_true(all(c("alternative", "method", "data.name", "statistic", "p.value",
                    "cp.location", "lrv") %in% names(y)))
  expect_true(all(c("method", "param", "value") %in% names(y$lrv)))
})

test_that("CUSUM test for changes in the scale is performed correctly", 
{
  suppressWarnings(require(mvtnorm))
  
  n <- 500
  m <- 100
  N <- n + m
  k <- m + floor(n * 0.5)
  n1 <- N - k
  
  rho <- c(0.4, -0.9)
  
  ## rho
  theta1 <- 0
  theta2 <- 0
  theta <- cbind(c(theta1, 0), c(0, theta2))
  q <- rho * sqrt( (theta1^2 + 1) * (theta2^2 + 1) / (theta1 * theta2 + 1))
  S0 <- cbind(c(1, q[1]), c(q[1], 1))
  S1 <- cbind(c(1, q[2]), c(q[2], 1))
  
  p <- replicate(1000, 
  {
    e0 <- rmvt(k, S0, 5)
    e1 <- rmvt(n1, S1, 5)
    e <- rbind(e0, e1)
    x <- matrix(numeric(N * 2), ncol = 2)
    x[1, ] <- e[1, ]
    invisible(sapply(2:N, function(i) x[i, ] <<- e[i, ] + theta %*% e[i-1, ]))
    
    x <- x[-(1:m), ]
    cor_cusum(x, "rho")$p.value
  })
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  ## tau
  S0 <- cbind(c(1, rho[1]), c(rho[1], 1))
  S1 <- cbind(c(1, rho[2]), c(rho[2], 1))
  
  p <- replicate(1000, 
  {
    e0 <- rmvt(k, S0, 5)
    e1 <- rmvt(n1, S1, 5)
    e <- rbind(e0, e1)
    x <- matrix(numeric(N * 2), ncol = 2)
    x[1, ] <- e[1, ]
    # AR(1):
    invisible(sapply(2:N, function(i) x[i, ] <<- 0.8 * x[i-1, ] + e[i, ]))
    
    x <- x[-(1:m), ]
    cor_cusum(x, version = "tau")$p.value
  })
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  
  # correct change point location
  ## rho
  set.seed(1895)
  rho <- c(0.9, -0.9)
  q <- rho * sqrt( (theta1^2 + 1) * (theta2^2 + 1) / (theta1 * theta2 + 1))
  S0 <- cbind(c(1, q[1]), c(q[1], 1))
  S1 <- cbind(c(1, q[2]), c(q[2], 1))
  
  e0 <- rmvt(k, S0, 5)
  e1 <- rmvt(n1, S1, 5)
  e <- rbind(e0, e1)
  x <- matrix(numeric(N * 2), ncol = 2)
  x[1, ] <- e[1, ]
  invisible(t(sapply(2:N, function(i) x[i, ] <<- e[i, ] + theta %*% e[i-1, ])))
  
  x <- x[-(1:m), ]
  expect_equal(attr(cor_cusum(x, "rho")$statistic, "cp-location"), 250, tolerance = 0.2)
  
  ## tau
  expect_equal(attr(cor_cusum(x, "tau")$statistic, "cp-location"), 250, tolerance = 0.2)
  
  
  ## maybe some more tests
  ## best to be checked graphically:
  ## hist(p)
})