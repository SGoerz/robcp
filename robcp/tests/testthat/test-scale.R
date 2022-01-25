context("Scale test")

test_that("Output of scale_stat has the correct format", 
{
  x <- 1:6
  y <- scale_stat(x)
  
  expect_true(is.numeric(y))
  expect_equal(length(y), 1)
  
  expect_error(scale_stat(x, version = "a"))
  expect_error(scale_stat(1))
})


test_that("scale_stat is computed correctly", 
{
  x <- c(92, 49, 48, 57, 27)
  b_n <- 2
  
  # empVar:
  y <-  5034.6 / sqrt(5 * lrv(x,control = list(b_n = b_n, mean = mean(x),
                                               var = var(x), kFun = "SFT", 
                                               version = "empVar")))
  z <- scale_stat(x, version = "empVar", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # MD:
  y <-  49 / sqrt(5 * lrv(x, control = list(b_n = b_n, mean = median(x), kFun = "SFT",
                                               var = mean(abs(x - median(x))) * 5 / 4, 
                                               version = "MD")))
  z <- scale_stat(x, version = "MD", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # GMD:
  y <-  30.4 / sqrt(5 * lrv(x, "kernel", 
                            control = list(b_n = b_n, mean = 0, kFun = "SFT",
                                           var = sum(sapply(2:5, function(j) sum(abs(x[j] - x[1:(j-1)])))) / 10, 
                                           version = "GMD")))
  z <- scale_stat(x, version = "GMD", control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # MAD:
  y <-  27 / sqrt(5 * lrv(x, "kernel", 
                          control = list(b_n = b_n, mean = median(x), 
                                         var = mad(x, constant = 1), version = "MAD")))
  z <- scale_stat(x, version = "MAD", constant = 1, control = list(b_n = b_n))
  attributes(z) <- NULL
  expect_equal(z, y)
  
  # QBeta:
  x <- c(6, 2, 1, 10, 5, 0)
  beta <- 0.9
  y <-  12 / sqrt(6 * lrv(x, "kernel", control = list(version = "QBeta",
                              kFun = "SFT", b_n = b_n, mean = beta, var = 9)))
  z <- scale_stat(x, version = "QBeta", control = list(b_n = b_n), beta = beta)
  attributes(z) <- NULL
  expect_equal(z, y)
})


test_that("Output of scale_cusum has the correct format", 
{
  x <- rnorm(10)
  res <- suppressWarnings(scale_cusum(x))
  
  expect_equal(class(res), "htest")
  expect_equal(res$alternative, "two-sided")
  expect_equal(res$method, "CUSUM test for scale changes")
})


test_that("CUSUM test for changes in the scale is performed correctly", 
{
  # empVar:
  suppressWarnings({p <- replicate(200, 
  {
    x <- rnorm(200)
    x[101:200] <- x[101:200] * 3
    scale_cusum(x, control = list(b_n = NA))$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  
  # MD:
  suppressWarnings({p <- replicate(200, 
  {
    x <- rnorm(200)
    x[101:200] <- x[101:200] * 3
    scale_cusum(x, version = "MD", control = list(b_n = NA))$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  
  # GMD:
  suppressWarnings({p <- replicate(200, 
  {
    x <- rnorm(200)
    x[101:200] <- x[101:200] * 10
    scale_cusum(x, version = "GMD")$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  
  # MAD:
  suppressWarnings({p <- replicate(200, 
  {
    x <- rnorm(200)
    x[101:200] <- x[101:200] * 3
    scale_cusum(x, version = "MAD")$p.value
  })})

  expect_equal(mean(p < 0.05), 1, tolerance = 0.01)
  
  # QBeta:
  suppressWarnings({p <- replicate(200, 
  {
    x <- rnorm(200)
    x[101:200] <- x[101:200] * 3
    scale_cusum(x, version = "QBeta", control = list(b_n = 10), beta = 0.9)$p.value
  })})
  
  expect_equal(mean(p < 0.05), 1, tolerance = 0.1)
  
  # correct change point location
  x <- rnorm(100)
  x[50:100] <- x[50:100] * 3
  expect_equal(attr(scale_cusum(x, version = "empVar")$statistic, "cp-location"), 50, tolerance = 1)
  
  ## maybe some more tests
  ## best to be checked graphically:
  ## hist(p)
})