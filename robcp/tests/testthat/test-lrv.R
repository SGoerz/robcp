context("lrv")

test_that("output has the correct format",
{
  x <- rnorm(100)
  y <- lrv(x, control = list(b_n = 10))
  expect_equal(class(y), "numeric")
  expect_equal(length(y), 1)
  
  expect_error(lrv(x, control = list(b_n = 0)))
  expect_warning(lrv(x, control = list(kFun = "abc")))
  
  x <- ts(x)
  y <- lrv(x, control = list(b_n = 10))
  expect_equal(class(y), "numeric")
  expect_equal(length(y), 1)
  
  X <- matrix(rnorm(100), ncol = 10)
  Y <- lrv(X, control = list(b_n = 5))
  expect_is(Y, "matrix")
  expect_true(is.numeric(Y))
  expect_equal(dim(Y), c(10, 10))
  expect_equal(Y, t(Y))
  
  X <- ts(X)
  Y <- lrv(X, control = list(b_n = 5))
  expect_is(Y, "matrix")
  expect_true(is.numeric(Y))
  expect_equal(dim(Y), c(10, 10))
  
  expect_error(lrv(X, control = list(b_n = 11)))
  expect_error(lrv(X, control = list(b_n = 0)))
})

test_that("correct warnings and errors", 
{
  x <- c(1, -1, 1, -1)
  expect_warning(lrv(x, control = list(kFun = "FT", b_u = 1)))
  
  ## other names than in control
  expect_warning(lrv(x, control = list(fun = "TH")))
  
  ## wrong parameters specified
  expect_error(lrv(x, control = list(b_n = 5)))
  expect_error(lrv(x, method = "subsampling", control = list(l = -1)))
  expect_error(lrv(x, method = "bootstrap", control = list(B = 0)))
  expect_error(lrv(x, method = "bootstrap", control = list(l = 5)))
  
  ## methods getting confused
  
  ## kernel functions for bootstrap
  expect_warning(lrv(x, method = "bootstrap", control = list(kFun = "FT", l = 1)))
})

test_that("kernel-based estimation is correctly computed", 
{
  x <- 1:5
  #y <- c(2, 0.8, -0.2, -0.8, -0.8)
  
  res1 <- 2
  res2 <- 2 + 2*0.8
  res3 <- 2 + 2*0.8 - 4/3*0.2
  res4 <- 2 + 2*0.8 - 2*0.2 - 0.8
  res5 <- 2 + 2*0.8 - 2*0.2 - 2*0.8*0.8 - 2*0.4*0.8
  
  expect_equal(lrv(x, control = list(b_n = 1, kFun = "FT")), res1)
  expect_equal(lrv(x, control = list(b_n = 2, kFun = "FT")), res2)
  expect_equal(lrv(x, control = list(b_n = 3, kFun = "FT")), res3)
  expect_equal(lrv(x, control = list(b_n = 4, kFun = "FT")), res4)
  expect_equal(lrv(x, control = list(b_n = 5, kFun = "FT")), res5)
  
  X <- cbind(x, seq(-10, 10, 5))

  expect_equal(lrv(X, control = list(b_n = 2, kFun = "FT")), matrix(c(3.6, 18, 18, 90), ncol = 2))
  expect_equal(lrv(X, control = list(b_n = 1, kFun = "FT")), matrix(c(2, 10, 10, 50), ncol = 2))
})

test_that("kernel-based estimation for the scale is correctly computed", 
{
  x <- arima.sim(list(ar = 0.5), 5)
  b_n <- 3
  
  ft <- 1:(b_n-1) / b_n
  ft <- ifelse(abs(ft) < 0.5, 1, ifelse(abs(ft) < 1, 2 - 2 * abs(ft), 0))
  
  lrvTest <- function(z, b_n)
  {
    ac <- acf(z, plot = FALSE, type = "covariance", demean = FALSE,
              lag.max = b_n - 1)$acf[, , 1]
    a <- sum(2 * ac[-1] * ft) + ac[1]
    if(a > 0) return(a)
    return(var(z) * (length(z) - 1) / length(z))
  }

  # empVar:
  m <- mean(x)
  y <- suppressWarnings({lrv(x, method = "kernel",
           control = list(version = "empVar", b_n = b_n, kFun = "FT"))})
  
  newx <- (x - m)^2
  expect_equal(y, lrvTest(newx - mean(newx), b_n))
  
  # MD:
  m <- median(x)
  y <- suppressWarnings({lrv(x, method = "kernel",
           control = list(version = "MD", b_n = b_n, kFun = "FT"))})
  
  newx <- abs(x - m)
  expect_equal(y, lrvTest(newx - mean(newx), b_n))
  
  # GMD:
  v <- sum(sapply(2:5, function(j) sum(abs(x[j] - x[1:(j-1)])))) / 10
  y <- suppressWarnings({lrv(x, method = "kernel",
                             control = list(version = "GMD", b_n = b_n,
                                            kFun = "FT"))})
  
  newx <- sapply(seq_along(x), function(i) mean(abs(x[-i] - x[i])))
  expect_equal(y, 4 * lrvTest(newx - mean(newx), b_n))
  
  
  x <- c(85, 89, 36, 12, 51, 39, 24)
  # # MAD:
  # m <- 39 
  # v <- 15
  # y <- lrv(x, "kernel", 
  #          control = list(version = "MAD", scale = 15, b_n = b_n, kFun = "FT"))
  # mad_f <- sum(3 / 4 * (1 - (c(-12, 12, -3, -15, 0) / 38 / 7^(-1/3))^2)) / 38 / 7^(2/3)
  # 
  # expect_equal(y, lrvTest(as.numeric(abs(x - m) <= v) - 0.5, b_n) / mad_f)
  # 
  # Qalpha:
  beta <- 0.5
  v <- 34
  y <- lrv(x, "kernel", control = list(alpha_Q = beta, version = "Qalpha",
                        b_n = b_n, kFun = "FT"))

  qalpha_u <- sum(sapply(2:7, function(j)
  {
    temp <- (abs(x[1:(j-1)] - x[j]) - v) / IQR(x) / 7^(-1/3)
    temp <- ifelse(abs(temp) < 1, 3 / 4 * (1 - temp^2), 0)
    return(sum(temp))
  })) / (21 * 7^(-1/3) * IQR(x))

  expect_equal(y, lrvTest(sapply(x, function(xi)
    mean(as.numeric(abs(x - xi) <= v))) - beta, b_n) * 4 / qalpha_u^2)
})

test_that("kernel-based estimation for the correlation is computed correctly", 
{
  ## rho
  x <- matrix(c(80, 23, 74, 93, 39, 48, 26, 77), ncol = 2)
  
  res1 <- 144 * 75 / 16^3
  res2 <- 144 * (75/16^3 - 50/16^3)
  res3 <- 144 * (75/16^3 - 50/16^3 - 4/3 * 73/16^3)
  res4 <- 144 * (75/16^3 - 50/16^3 - 146/16^3 - 121/16^3)
  
  expect_equal(lrv(x, control = list(b_n = 1, kFun = "FT", version = "rho")), res1)
  expect_equal(lrv(x, control = list(b_n = 2, kFun = "FT", version = "rho")), res2)
  expect_equal(lrv(x, control = list(b_n = 3, kFun = "FT", version = "rho")), res3)
  expect_equal(lrv(x, control = list(b_n = 4, kFun = "FT", version = "rho")), res4)
  
  
  ## tau
  x <- matrix(c(47, 82, 13, 53, 45, 15, 23, 19, 63, 88), ncol = 2)
  v <- 0
  
  res1 <- 4 * 8 / 5^3
  res3 <- 4 * (8 + 0 - 16/3) / 5^3
  
  expect_equal(lrv(x, control = list(b_n = 1, kFun = "FT", version = "tau")), res1)
  expect_equal(lrv(x, control = list(b_n = 3, kFun = "FT", version = "tau")), res3)
})

test_that("subsampling estimation is correctly computed", 
{
  n <- 100
  x <- arima.sim(model = list(ar = 0.5), n)
  y <- ecdf(x)(x)
  l <- 10
  
  ## overlapping
  res1 <- sapply(0:(n-l), function(i)
  {
    abs(sum(y[(i+1):(i+l)] - 0.5))
  })
  res1 <- (sum(res1) * sqrt(pi) / (sqrt(2*l) * (n - l + 1)))^2
  
  expect_equal(lrv(x, method = "subsampling", 
                   control = list(l = l, overlapping = TRUE, distr = TRUE)), 
               res1)
  
  ## non-overlapping 
  a <- floor(n / l) ### ??????????
  meanx <- mean(x)
  res2 <- sapply(1:a, function(i) 
  {
    (sum(x[((i-1)*l+1):(i*l)]) - l * meanx)^2
  })
  res2 <- sum(res2) / l / a
  
  expect_equal(lrv(x, method = "subsampling", 
                   control = list(l = l, overlapping = FALSE, distr = FALSE)), 
               res2)
  
  ## non-overlapping & distr
  meany <- mean(y)
  res3 <- sapply(1:a, function(i)
  {
    abs(sum(y[((i-1)*l+1):(i*l)]) - l * meany)
  })
  res3 <- (sum(res3) * sqrt(pi/2) / sqrt(l) / a)^2
  
  expect_equal(lrv(x, method = "subsampling", 
                   control = list(l = l, overlapping = FALSE, distr = TRUE)), 
               res3)
})

test_that("bootstrap estimation is correctly computed", 
{
  #skip_if_not_installed("mvtnorm", minimum_version = NULL)
  require("mvtnorm")
  x <- rnorm(5)
  sigma <- matrix(c(3, 2, 1, 0, 0, 
                    2, 3, 2, 1, 0,
                    1, 2, 3, 2, 1,
                    0, 1, 2, 3, 2,
                    0, 0, 1, 2, 3), ncol = 5) / 3
  set.seed(1895)
  dwb <- replicate(1000, 
  {
    eps <- rmvnorm(1, rep(0, 5), sigma)
    mean(mean(x) + (x - mean(x)) * eps)
  })
  
  expect_equal(var((dwb - mean(x)) * sqrt(5)), 
               lrv(x, method = "bootstrap", 
                   control = list(l = 3, B = 1000, seed = 1895)))

})