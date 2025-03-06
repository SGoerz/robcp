context("modifChol")

test_that("output has the correct format", 
{
  X1 <- diag(4)
  Y1 <- modifChol(X1)
  X2 <- matrix(c(1, 1, 2, 2, 3, 3), ncol = 2)
  
  expect_is(Y1, "matrix")
  expect_true(is.numeric(Y1))
  expect_equal(dim(Y1), c(4, 4))
  expect_true(all(Y1[lower.tri(Y1)] == 0))
  
  expect_error(modifChol(X2))
})

test_that("Cholesky decomposition is correctly computed", 
{
  X1 <- diag(4)
  Y1 <- X1

  ## limit case
  X2 <- matrix(c(3, 2, 2, 1), ncol = 2)
  Y2 <- modifChol(X2)
  
  X3 <- matrix(c(4, 2, 2, 1), ncol = 2)
  Y3 <- modifChol(X3)
  
  X4 <- matrix(c(5, 2, 2, 1), ncol = 2)
  Y4 <- modifChol(X4)
  
  X5 <- diag(c(3, 2, 1)) + matrix(4, ncol = 3, nrow = 3)
  Y5 <- chol(X5)

  expect_equal(modifChol(X1), Y1)
  expect_equal(t(Y2) %*% Y2, X2, tolerance = 0.3)
  expect_equal(t(Y3) %*% Y3, X3)
  expect_equal(t(Y4) %*% Y4, X4)
  expect_equal(modifChol(X5), Y5)
})