library(BBSsize)

h <- BBSsize::granby
ranges <- rwar:::test_ranges

test_that("regional nm matches species and differs size", {

  nm_results <- regional_null_model(h, ranges, 1994, 1, 1988:1992, 2014:2018, 1977)

  actual_results <- all_core_analyses(h, 1988:1992, 2014:2018, 1977)

  expect_equal(nm_results$bcd, actual_results$bcd)
  expect_equal(nm_results$sp_turnover, actual_results$sp_turnover)
  expect_false(nm_results$isd_turnover == actual_results$isd_turnover)


  }
)


test_that("local shuffle matches manual shuffle", {

  null_shuffle <- shuffle_regional(h, ranges, 1994)

  set.seed(1994)
  manual_shuffle <- sample(colnames(h$abundance), size = length(colnames(h$abundance)), replace = F)
  manual_dat <- h
  colnames(manual_dat$abundance) <- manual_shuffle

  expect_true(all(colnames(null_shuffle$abundance) == manual_shuffle))

  expect_true(all(all_core_analyses(null_shuffle, 1988:1992, 2014:2018, 1977) == all_core_analyses(manual_dat, 1988:1992, 2014:2018, 1977)))

  }
)

test_that("local wrapper works as intended", {

  several_null <- local_null_model_wrapper(h, 3, 1988:1992, 2014:2018)

  expect_equal(length(unique(several_null$isd_seed)), 3)
  expect_equal(length(unique(several_null$null_mod_seed)), 3)
  expect_equal(length(unique(several_null$sp_turnover)), 1)
  expect_equal(length(unique(several_null$isd_turnover)), 3)

}
)

