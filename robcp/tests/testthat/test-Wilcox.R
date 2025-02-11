context("Wilcoxon-Mann-Whitney")

test_that("output of wilcox_stat has the correct format", 
{
  x <- 1:5
  # expect_error(wilcox_stat(x))
  # expect_error(wilcox_stat(x, method = "bootstrap"))
  
  y <- wilcox_stat(x, control = list(l_n = 2))
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(wilcox_stat(x, b = 3L))
  
  expect_warning(wilcox_stat(x, h = 1, control = list(distr = FALSE, l_n = 2)))
  
  expect_error(wilcox_stat(x, h = "a"))
  
  testStructure(wmw_test, "kernel")
  testStructure(wmw_test, "subsampling")
  testStructure(wmw_test, "bootstrap")
})

test_that("wilcox_stat returns the correct value", 
{
  x <- c(37, 40, 9, 14, 37, 3)
  t1 <- 3.5 / sqrt(6)^3 / 
    sqrt(lrv(x, method = "subsampling", 
             control = list(overlapping = FALSE, distr = TRUE, l_n = 2)))
  
  y <- wilcox_stat(x, 1, method = "subsampling", 
                   control = list(overlapping = FALSE, distr = TRUE, l_n = 2))
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
                                                distr = FALSE, l_n = 2)))
  y <- wilcox_stat(x, 2, "subs", list(overlapping = FALSE, distr = FALSE, l_n = 2))
  attributes(y) <- NULL
  expect_equal(y, t2)
  
  x <- c(15, 20, 13, 1, 16, 22)
  t2 <- 54 / sqrt(6)^3 / sqrt(lrv(x, "subs", list(overlapping = FALSE, 
                                                distr = FALSE)))
  y <- wilcox_stat(x, 2, "subs", list(overlapping = FALSE, distr = FALSE))
  attributes(y) <- NULL
  expect_equal(y, t2)
  
  # wmw_test and CUSUM test are equal for h = 2 and the kernel-based long run
  # variance estimation
  x <- rnorm(100)
  y <- wmw_test(x, h = 2, method = "kernel",
                control = list(kFun = "TH", b_n = 5))$statistic
  attr(y, "names") <- NULL
  expect_equal(y, CUSUM(x, method = "kernel", control = list(b_n = 5)))
  
  # wmw_tests are equal for h = 2 and h = function(x, y) x - y
  y1 <- wmw_test(x, h = 2, method = "kernel")
  y2 <- wmw_test(x, h = function(x, y) x - y, method = "kernel")
  expect_equal(y1, y2)
  
  # wmw_tests are equal for h = 1 and h = function(x, y) as.numeric(x < y) - 0.5
  y1 <- wmw_test(x, h = 1, method = "kernel", control = list(distr = TRUE, b_n = 2))
  y2 <- wmw_test(x, h = function(x, y) as.numeric(x < y) - 0.5, method = "kernel", 
                 control = list(distr = TRUE, b_n = 2))
  expect_equal(y1, y2)
})


test_that("The output of wmw_test has the correct format",
{
  x <- rnorm(10)
  res <- suppressWarnings(wmw_test(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "Wilcoxon-Mann-Whitney change point test")
})

test_that("Wilcoxon-Mann-Whitney change point test is performed correctly", 
{
  ## simulation might run too long
  skip_on_cran()
  suppressWarnings({p <- replicate(200, 
  {
   x <- rnorm(200)
   x[101:200] <- x[101:200] + 1
   wmw_test(x)$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.001)
  
  
  skip_on_cran()
  suppressWarnings({p <- replicate(200, 
  {
   x <- rnorm(200)
   x[101:200] <- x[101:200] + 1
   wmw_test(x, h = 2L)$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.001)
  
  # correct change point location
  x <- rnorm(100)
  x[50:100] <- x[50:100] + 10
  expect_equal(attr(wmw_test(x)$statistic, "cp-location"), 50, tolerance = 1)
  
  ## maybe some more tests
  ## best to be checked graphically:
  ## hist(p)
})
