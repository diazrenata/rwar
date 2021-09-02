library(BBSsize)

h <- BBSsize::hartland
h_isd <- BBSsize::simulate_isd_ts(h, isd_seed = 1977)
h_svs <- get_annual_svs(h_isd$isd)

test_that("lm on full timeseries works", {

  ts_lm <- fit_all_timeseries_lms(h_svs)

  abund_lm <- lm(abundance ~ year, data = h_svs)

  abund_p <- anova(abund_lm)$`Pr(>F)`[1]

  abund_summary <- summary(abund_lm)

  abund_r2 <- abund_summary$r.squared

  abund_coefs <- coef(abund_lm)

  abund_slope <- abund_coefs["year"]

  abund_fitted_ratio <- abund_lm$fitted.values[25] / abund_lm$fitted.values[1]

  expect_true(abund_p == ts_lm$p_ts_abundance)
  expect_true(abund_r2 == ts_lm$r2_ts_abundance)
  expect_true(abund_slope == ts_lm$slope_ts_abundance)
  expect_true(abund_fitted_ratio == ts_lm$fitted_ratio_ts_abundance)

  expect_false(anyNA(ts_lm))
  expect_false(any(ts_lm[,1:5] > 1) )
  expect_false(any(ts_lm[,1:5] < 0) )

  })

test_that("lm on caps works", {

  h_caps <- pull_caps(h_svs)

  caps_lm <- fit_all_caps_lms(h_caps)

  abund_lm <- lm(abundance ~ timeperiod, data = h_caps)

  abund_p <- anova(abund_lm)$`Pr(>F)`[1]

  abund_summary <- summary(abund_lm)

  abund_r2 <- abund_summary$r.squared

  abund_coefs <- coef(abund_lm)

  abund_slope <- abund_coefs["timeperiodend"]

  abund_fitted_ratio <- abund_lm$fitted.values[10] / abund_lm$fitted.values[1]

  expect_true(abund_p == caps_lm$p_caps_abundance)
  expect_true(abund_r2 == caps_lm$r2_caps_abundance)
  expect_true(abund_slope == caps_lm$slope_caps_abundance)
  expect_true(abund_fitted_ratio == caps_lm$fitted_ratio_caps_abundance)

  expect_false(anyNA(caps_lm))
  expect_false(any(caps_lm[,1:5] > 1) )
  expect_false(any(caps_lm[,1:5] < 0) )

})
