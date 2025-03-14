context("psi function")

test_that("output is of the same format as input", 
{
  testFormat(psi)
  testFormat(function(y) psi(y, fun = "HCm"))
  
  matrixDim <- function(fun)
  {
    m <- 10
    X <- matrix(rnorm(1000), ncol = m)
    Y <- psi(X, fun)
    
    cols <- switch(fun, "HCm" = m * (m + 1) / 2, "HCg" = m * (m + 1) / 2, 
                   "SCm" = m * (m - 1) / 2, "SCg" = m * (m + 1) / 2 - 1, m)
    
    expect_equal(ncol(Y), cols)
    expect_equal(class(Y), class(X))
    expect_equal(nrow(Y), nrow(X))
    
    X <- ts(X)
    Y <- psi(X, fun)
    expect_equal(ncol(Y), cols)
    expect_equal(class(Y), class(X))
    expect_equal(nrow(Y), nrow(X))
  }
  
  matrixDim("HLm")
  matrixDim("SCm")
  matrixDim("SCg")
  
  expect_error(psi(1:5, fun = "xyz"))
})


test_that("Time series are correctly transformed", 
{
  X <- cbind(c(-1, 3, 6, 0, -1), c(1, 1, 2, 3, 3), c(-50, 1, 2, 1, 50))
  X.tilde <- cbind(X[, 1] / 1.4826, (X[, 2] - 2) / 1.4826, (X[, 3] - 1) / 1.4826)
  
  X.HLm <- X.tilde
  X.HLm[X.tilde > 1.5] <- 1.5
  X.HLm[X.tilde < -1.5] <- -1.5
  
  temp <- sqrt(rowSums(X.tilde^2))
  X.HLg <- t(sapply(1:5, function(i)
  {
    if(temp[i] == 0) return(0)
    if(temp[i] <= 1.5) return(X.tilde[i, ])
    else return(1.5 * X.tilde[i, ] / temp[i])
  }))
  
  X.SLm <- sign(X.tilde)
  
  X.SLg <- t(sapply(1:5, function(i)
  {
    if(temp[i] == 0) return(0) else return(X.tilde[i, ] / temp[i])
  }))
  
  expect_equal(psi(X, fun = "HLm", k = 1.5), X.HLm)
  expect_equal(psi(X, fun = "HLg", k = 1.5), X.HLg)
  expect_equal(psi(X, fun = "SLm", k = 1.5), X.SLm)
  expect_equal(psi(X, fun = "SLg", k = 1.5), X.SLg)
  
  
  X.HCm <- t(apply(X.HLm, 1, function(x) 
  {
    a <- x %*% t(x)
    a[lower.tri(a, diag = TRUE)]
  }))
  
  X.HCg <- t(apply(X.HLg, 1, function(x) 
  {
    a <- x %*% t(x)
    a[lower.tri(a, diag = TRUE)]
  }))
  
  X.SCm <- t(apply(X.SLm, 1, function(x) 
  {
    a <- x %*% t(x)
    a[lower.tri(a, diag = FALSE)]
  }))
  
  X.SCg <- t(apply(X.SLg, 1, function(x) 
  {
    a <- x %*% t(x)
    a <- a[lower.tri(a, diag = TRUE)]
    a <- a[-length(a)]
  }))
  
  expect_equal(psi(X, fun = "HCm", k = 1.5), X.HCm)
  expect_equal(psi(X, fun = "HCg", k = 1.5), X.HCg)
  expect_equal(psi(X, fun = "SCm", k = 1.5), X.SCm)
  expect_equal(psi(X, fun = "SCg", k = 1.5), X.SCg)
})