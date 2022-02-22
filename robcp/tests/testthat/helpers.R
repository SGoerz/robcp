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

testStructure <- function(fun, method)
{
  x <- rnorm(100)
  y <- fun(x, method = method)
  
  expect_true(is.list(y))
  expect_true(is(y$lrv, "list"))
  expect_equal(class(y), "htest")
  
  expect_true(all(c("alternative", "method", "data.name", "statistic", "p.value",
                    "cp.location", "lrv") %in% names(y)))
  expect_true(all(c("method", "param", "value") %in% names(y$lrv)))
}
