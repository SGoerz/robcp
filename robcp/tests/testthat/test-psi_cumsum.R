context("psi_cumsum")

test_that("output is of the same format as input", 
{
  testFormat(psi_cumsum)
})

test_that("Cumulative sums are computed corretly",
{
  x <- c(0, 0, 1, 50, 50)
  mad(x, constant = 1)
  psi(x, constant = 1)
  X <- cbind(1:5, c(0, 0, 1, 50, 50))
  mad(1:5, constant = 1)
  
  expect_equal(psi_cumsum(X, constant = 1), 
               Y <- cbind(c(-1.5, -2.5, -2.5, -1.5, 0), c(-1, -2, -2, -0.5, 1)))
})