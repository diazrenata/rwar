test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

g <- BBSsize::granby

gsims <- ssims_wrapper(g, "actual")

test_that("fit works", {

  alm <- fit_stanlm(gsims)

  expect_true(length(alm) == 4)

  call_chr <-  toString(alm$te_stanlms$te_stanlm_full$call)

  expect_true(call_chr == "rstanarm::stan_glm, total_energy ~ (timeperiod * source), some_sims, 8000, 4")

  comps <- compare_both_stanarms(alm)

  expect_true(nrow(comps == 6))
  expect_false(anyNA(comps))

  best <- loo_select(comps)

  expect_true(nrow(best) == 2)
  expect_true(all(best$model_complexity == c(1,3)))

  best_draws <- winner_draws(best, alm)

  expect_true(nrow(best_draws) == 8000)
  expect_true(all(unique(best_draws$modtype) == c("te_stanlm_notime", "tb_stanlm_full")))

  qis <- winner_qis(best_draws)

  expect_true(nrow(qis) == 4)
  expect_true(all(unique(qis$modtype) == c("tb_stanlm_full", "te_stanlm_notime")))

  diag <- extract_diagnostics(alm)

  expect_true(all(diag$divergent_sum == 0))
  expect_true(anyNA(diag$rhat_timeperiodend))
  expect_false(anyNA(diag$`rhat_(Intercept)`))
})



