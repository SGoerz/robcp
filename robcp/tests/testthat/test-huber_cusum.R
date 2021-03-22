context("huber_cusum")

test_that("The output of huber_cusum has the correct format",
{
  x <- rnorm(10)
  res <- huber_cusum(x)
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "Huberized CUSUM test")
})

test_that("Huberized CUSUM test is performed correctly", 
{
  p <- replicate(1000, 
  {
    x <- rnorm(1000)
    x[501:1000] <- x[501:1000] + 1
    huber_cusum(x)$p.value
  })
  
  expect_equal(mean(p < 0.05), 1)
  
  
  expect_true(FALSE)
})