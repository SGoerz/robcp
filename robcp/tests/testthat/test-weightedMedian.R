context("weightedMedian")

test_that("Output of weightedMedian is correct",
{
  expect_equal(weightedMedian(c(1, 5, 9), c(5, 1, 1)), 1)
  expect_equal(weightedMedian(1:3, c(1, 1, 1)), 2)
  expect_equal(weightedMedian(1:2, c(1, 1)), 2)
  expect_equal(weightedMedian(1:4, c(2, 2, 0, 2)), 2)
  
   
  skip_on_os(os = "solaris") 
  expect_error(weightedMedian(1:2, 1))
  expect_error(weightedMedian(1:2, -1:0))
})