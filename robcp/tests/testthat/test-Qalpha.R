context("Qalpha")

test_that("Output of Qalpha has the correct format", 
{
  n <- 5
  x <- sample(1:50, n)
  y <- Qalpha(x, 0.5)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), n-1)
  
  expect_error(Qalpha(x, "a"))
  expect_error(Qalpha(x, -0.4))
  expect_error(Qalpha(x, 4))
})

test_that("Qalpha returns the correct values", 
{
  x <- sample(1:50, 5)
  diffs <- sort(abs(unlist(sapply(1:4, function(i) x[i] - x[(i+1):5]))))

  expect_equal(sapply(1:10, function(alpha) Qalpha(x, alpha/10)[4]), diffs)
})