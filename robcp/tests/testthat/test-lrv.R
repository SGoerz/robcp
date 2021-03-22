context("lrv")

test_that("output has the correct format",
{
  x <- rnorm(100)
  y <- lrv(x, 10)
  expect_equal(class(y), "numeric")
  expect_equal(length(y), 1)
  
  expect_error(lrv(x, 0))
  
  x <- ts(x)
  y <- lrv(x, 10)
  expect_equal(class(y), "numeric")
  expect_equal(length(y), 1)
  
  X <- matrix(rnorm(100), ncol = 10)
  Y <- lrv(X, 5)
  expect_is(Y, "matrix")
  expect_true(is.numeric(Y))
  expect_equal(dim(Y), c(10, 10))
  expect_equal(Y, t(Y))
  
  X <- ts(X)
  Y <- lrv(X, 5)
  expect_is(Y, "matrix")
  expect_true(is.numeric(Y))
  expect_equal(dim(Y), c(10, 10))
  
  expect_error(lrv(X, 11))
  expect_error(lrv(X, 0))
  
  x <- c(1, -1, 1, -1)
  expect_warning(lrv(x))
})

test_that("long run variance is correctly computed", 
{
  x <- 1:5
  y <- c(2, 0.8, -0.2, -0.8, -0.8)
  
  res1 <- 2
  res2 <- 2 + 2*0.8
  res3 <- 2 + 2*0.8 - 4/3*0.2
  res4 <- 2 + 2*0.8 - 2*0.2 - 0.8
  res5 <- 2 + 2*0.8 - 2*0.2 - 2*0.8*0.8 - 2*0.4*0.8
  
  expect_equal(lrv(x, 1), res1)
  expect_equal(lrv(x, 2), res2)
  expect_equal(lrv(x, 3), res3)
  expect_equal(lrv(x, 4), res4)
  expect_equal(lrv(x, 5), res5)
  
  X <- cbind(x, seq(-10, 10, 5))

  expect_equal(lrv(X, 2), matrix(c(3.6, 18, 18, 90), ncol = 2))
  expect_equal(lrv(X, 1), matrix(c(2, 10, 10, 50), ncol = 2))
})