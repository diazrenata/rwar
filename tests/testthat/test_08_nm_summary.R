library(BBSsize)


test_that("ses works correctly", {

  expect_equal(0, ses(50, 0:100))

  rvect <- runif(10000, 0, 1000000)

  expect_equal(0, ses(mean(rvect), rvect))
  expect_equal(1, ses(mean(rvect) + sd(rvect), rvect))
  expect_equal(-1, ses(mean(rvect) - sd(rvect), rvect))

  }
)

test_that("percentile works correctly", {

  expect_equal(.49, percentile_score(50, 1:100))
  expect_equal(0, percentile_score(1, 1:100))
  expect_equal(.99, percentile_score(100, 1:100))

  set.seed(1977)
  rvect <- runif(10, 0, 10000)

  expect_equal(.5, percentile_score(median(rvect), rvect))
  expect_equal(0, percentile_score(min(rvect), rvect))
  expect_equal(.9, percentile_score(max(rvect), rvect))

  tievect <- c(1, 2, 3, 4, 4, 4, 4, 4, 5, 6)

  expect_equal(.3, percentile_score(median(tievect), tievect))
  expect_equal(0, percentile_score(min(tievect), tievect))
  expect_equal(.9, percentile_score(max(tievect), tievect))

  }
)

test_that('summary works correctly', {


  h <- BBSsize::hartland
  ar_h <- all_core_analyses(h, 1994:1998, 2014:2018, 1977)
  n_h <- local_null_model_wrapper(h, 5, 1994:1998, 2014:2018)

  g <- BBSsize::granby
  ar_g <- all_core_analyses(g, 1988:1992, 2014:2018, 1977)
  n_g <- local_null_model_wrapper(g, 5, 1988:1992, 2014:2018)

  nulls <- dplyr::bind_rows(n_h, n_g)
  ars <- dplyr::bind_rows(ar_h, ar_g)

  ns <- summarize_null_results(nulls, ars)

  expect_true(all(ns[1, ] == summarize_null_results(n_h, ar_h)))
  expect_true(all(ns[2, ] == summarize_null_results(n_g, ar_g)))

  expect_true(ns$ses[1] == ses(ar_h$isd_turnover, n_h$isd_turnover))
  expect_true(ns$perc[1] == percentile_score(ar_h$isd_turnover, n_h$isd_turnover))

})
