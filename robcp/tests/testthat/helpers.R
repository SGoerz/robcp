testFormat <- function(fun)
{
  x <- rnorm(100)
  y <- fun(x)
  expect_equal(class(y), class(x))
  expect_equal(length(y), length(x))
  
  x <- ts(x)
  y <- fun(x)
  expect_equal(class(y), class(x))
  expect_equal(length(y), length(x))
  expect_equal(attr(y, "tsp"), attr(x, "tsp"))
  
  X <- matrix(rnorm(1000), ncol = 10)
  Y <- fun(X)
  expect_equal(class(Y), class(X))
  expect_equal(nrow(Y), nrow(X))
  
  X <- ts(X)
  Y <- psi(X)
  expect_equal(class(Y), class(X))
  expect_equal(nrow(Y), nrow(X))
  expect_equal(attr(Y, "tsp"), attr(X, "tsp"))
}
