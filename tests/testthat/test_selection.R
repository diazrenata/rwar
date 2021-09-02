library(testthat)
library(cvlt)

abund_dat <- get_rodents_annual()


cv_fits <- fit_ldats_crossval(abund_dat, k = 2, lda_seed = 2, cpts = 0, nit = 50)


cv_fits2 <- fit_ldats_crossval(abund_dat, k = 2, lda_seed = 2, cpts = 1, nit = 50)


test_that("toplevel works", {

  expect_true(is.data.frame(cv_fits))
  expect_false(any(is.na(cv_fits)))

})

test_that("selection works", {

  both_fits <- dplyr::bind_rows(cv_fits, cv_fits2)

  best <- select_cvlt(both_fits)

  expect_true(nrow(best) == 1)

  expect_true(best$cpts == 1)

})

test_that("run_best works", {

  both_fits <- dplyr::bind_rows(cv_fits, cv_fits2)

  best <- select_cvlt(both_fits)

  best_mod <- run_best_model(abund_dat, best, nit = 50)

  expect_true(length(best_mod) == 4)

  })
