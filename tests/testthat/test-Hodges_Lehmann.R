context("Hodges_Lehmann")

test_that("output has the correct format", 
{
  skip_on_os(os = "solaris")
  
  x <- 1:5
  y <- HodgesLehmann(x, b_u = 1, control = list(b_n = 1, l_n = 1))
  
  # X <- matrix(1:9, ncol = 3)
  # Y <- HodgesLehmann(X)
  
  expect_true(is.numeric(y))
  # expect_true(is.numeric(Y))
  expect_equal(length(y), 1)
  # expect_equal(length(Y), 1)
  # expect_error(HodgesLehmann(x, b_u = 0))
  
  testStructure(hl_test, "kernel")
  testStructure(hl_test, "subsampling")
  # testStructure(hl_test, "bootstrap")
})

# test_that("u_hat computes the correct value",
# {
#   skip_on_os(os = "solaris")
#   
#   x <- rnorm(5)
#   expect_error(u_hat(x, -2, kFun = "FT"))
#   
#   b <- 2
#   y <- abs(apply(combn(x, 2), 2, diff) / b)
#   y <- ifelse(y > 1, 1, y)
#   y <- sum(1 - y)
#   expect_equal(2 / (5 * 4 * b) * y * 2 / 3, u_hat(x, b, kFun = "bartlett"))
# })

test_that("HodgesLehmann computes the correct value",
{
  skip_on_os("solaris")
  
  x <- c(14, 49, 50, 47, 28)
  l <- 2
  y <- sqrt(5) * 5.44 * u_hat(x - c(0, rep(34, 4))) / 
    sqrt(lrv(x, "kernel", control = list(l_n = l, overlapping = TRUE, distr = TRUE)))
  z <- HodgesLehmann(x, control = list(l_n = l, overlapping = TRUE, distr = TRUE))
  attributes(z) <- NULL
  
  expect_equal(z, y)
  
  x <- c(58, 2, 59, 26, 20, 88)
  l <- 3
  y <- sqrt(6) * 310 / 36 * u_hat(x - c(rep(0, 5), 62)) /
    sqrt(lrv(x, "subs", control = list(l_n = l, overlapping = TRUE, distr = TRUE)))
  z <- HodgesLehmann(x, method = "subsampling",
                     control = list(l_n = l, overlapping = TRUE, distr = TRUE))
  attributes(z) <- NULL
  
  expect_equal(z, y)
})

test_that("The output of hl_test has the correct format",
{
  skip_on_os(os = "solaris")
  
  x <- rnorm(10)
  res <- suppressWarnings(hl_test(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "Hodges-Lehmann change point test")
})

test_that("Hodges-Lehmann change point test is performed correctly", 
{
  skip_on_os(os = "solaris")
  skip_on_cran()
  
  ## simulation might run too long
  suppressWarnings({p <- replicate(200, 
  {
   x <- rnorm(200)
   x[101:200] <- x[101:200] + 1
   hl_test(x, b_u = 0.05, control = list(b_n = 10))$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.1)
  
  # correct change point location
  x <- rnorm(100)
  x[50:100] <- x[50:100] + 10
  expect_equal(attr(hl_test(x, b_u = 0.1)$statistic, "cp-location"), 50, tolerance = 1)
  
  ## maybe some more tests
  ## best to be checked graphically:
  ## hist(p)
})
