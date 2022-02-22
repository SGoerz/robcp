context("huber_cusum")

test_that("The output of huber_cusum has the correct format",
{
  x <- rnorm(10)
  res <- suppressWarnings(huber_cusum(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "Huberized CUSUM test")
  
  testStructure(huber_cusum, "kernel")
  testStructure(huber_cusum, "subsampling")
  testStructure(huber_cusum, "bootstrap")
})

test_that("Huberized CUSUM test is performed correctly", 
{
  ## simulation might run too long
  skip_on_cran()
  suppressWarnings({p <- replicate(5000, 
  {
    x <- rnorm(1000)
    x[501:1000] <- x[501:1000] + 1
    huber_cusum(x)$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.001)
  
  ## maybe some more test
  ## best to be checked graphically:
  ## hist(p)
})