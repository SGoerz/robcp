context("asymptotic cumulative distribution of CUSUM")

test_that("The output of pKSdist and pbessel has the correct format",
{
  x <- runif(1, 0, 10)
  y <- pKSdist(x)
  Y <- pBessel(x, 2L)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  expect_true(is.numeric(Y))
  expect_equal(length(Y), 1)
})

test_that("pKSdist returns the correct value", 
{
  expect_equal(pKSdist(4), 1)
  expect_equal(pKSdist(-1), 0)
  
  skip_on_cran() ## simulation takes a few seconds
  n <- 1001
  times <- seq(0, 1, length.out = n)
  res <- replicate(10000, 
  {
   dW <- rnorm(n) / sqrt(n)
   W <- cumsum(dW)
   
   B <- W - times * W[n] 
   
   max(abs(B))
  })
  
  expect_equal(pKSdist(0.5), mean(res <= 0.5), tolerance = 0.05)
  expect_equal(pKSdist(1), mean(res <= 1), tolerance = 0.05)
  expect_equal(pKSdist(1.5), mean(res <= 1.5), tolerance = 0.05)
})

test_that("pBessel returns the correct value",
{
  expect_equal(pBessel(2.114, 2), 0.9, tolerance = 1e-4)
  expect_equal(pBessel(3.396, 2), 0.99, tolerance = 1e-5)
  expect_equal(pBessel(32.624, 100), 0.9, tolerance = 1e-4)
  expect_equal(pBessel(36.783, 100), 0.99, tolerance = 1e-5)
})



