context("Wilcoxon")

test_that("output has the correct format", 
{
  x <- 1:5
  y <- wilcox_stat(x)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(wilcox_stat(x, b = 3L))
})

test_that("wilcox_stat returns the correct value", 
{
  x <- c(37, 40, 9, 14, 37, 3)
  t1 <- 4 / sqrt(6)^3 / 
    sqrt(lrv(x, method = "subsampling", 
             control = list(overlapping = FALSE, distr = TRUE)))
  
  y <- wilcox_stat(x, 1, method = "subsampling", 
                   control = list(overlapping = FALSE, distr = TRUE))
  attributes(y) <- NULL
  expect_equal(y, t1)
  
  x <- c(27, 29, 42, 50, 39, 40, 17)
  t1 <- 3 / sqrt(7)^3 / 
    sqrt(lrv(x, method = "subsampling", 
             control = list(overlapping = FALSE, distr = TRUE)))
  y <- wilcox_stat(x, 1, method = "subsampling", 
                   control = list(overlapping = FALSE, distr = TRUE))
  attributes(y) <- NULL
  expect_equal(y, t1)
  
  x <- c(21, 39, 1, 36, 32, 23, 14)
  t2 <- 88 / sqrt(7)^3 / 
    sqrt(lrv(x, method = "subs", control = list(overlapping = FALSE, 
                                                distr = FALSE)))
  y <- wilcox_stat(x, 2, "subs", list(overlapping = FALSE, distr = FALSE))
  attributes(y) <- NULL
  expect_equal(y, t2)
  
  x <- c(15, 20, 13, 1, 16, 22)
  t2 <- 54   / sqrt(6)^3 / sqrt(lrv(x, "subs", list(overlapping = FALSE, 
                                                distr = FALSE)))
  y <- wilcox_stat(x, 2, "subs", list(overlapping = FALSE, distr = FALSE))
  attributes(y) <- NULL
  expect_equal(y, t2)
})
