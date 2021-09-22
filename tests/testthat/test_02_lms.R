library(BBSsize)

h <- BBSsize::granby
h_isd <- BBSsize::simulate_isd_ts(h, isd_seed = 1977)
h_svs <- get_annual_svs(h_isd$isd)
h_caps_svs <- pull_caps(h_svs, 1988:1992, 2014:2018)

test_that("lm on full timeseries works", {

  ts_lm <- fit_all_timeseries_lms(h_svs)

  abund_lm <- lm(abundance ~ year, data = h_svs)

  abund_p <- anova(abund_lm)$`Pr(>F)`[1]

  abund_summary <- summary(abund_lm)

  abund_r2 <- abund_summary$r.squared

  abund_coefs <- coef(abund_lm)

  abund_slope <- abund_coefs["year"]

  abund_fitted_ratio <- abund_lm$fitted.values[18] / abund_lm$fitted.values[1]

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



test_that("rangescale works as intended", {

  rs <- rangescale(h_caps_svs$abundance)

  expect_true(min(rs) == 0)
  expect_true(max(rs) == 1)
  expect_false(anyNA(rs))

})


test_that("interaction lm works as intended", {

  expect_error(interaction_lms(h_caps_svs, scaling = "blah"))

  ilm <- interaction_lms(h_caps_svs)

  h_scaled <- h_caps_svs %>%
    dplyr::mutate(
      abundance = scale(sqrt(abundance)),
      biomass = scale(sqrt(biomass)),
      energy = scale(sqrt(energy))
    )

  expect_equivalent(mean(h_scaled$energy), 0)
  expect_equivalent(sd(h_scaled$energy), 1)

  h_scaled_long <- h_scaled %>%
    tidyr::pivot_longer(c(-year, -timeperiod), names_to = "currency", values_to = "val")%>%
    dplyr::filter(currency %in% c("energy", "abundance", "biomass"))

  h_scaled_energy <- dplyr::filter(h_scaled_long, currency == "energy")
  expect_true(all(h_scaled_energy$val == h_scaled$energy))

  h_ilm <- lm(val ~ timeperiod * currency, h_scaled_long)

  h_ilm_coef <- summary(h_ilm)$coefficients %>% as.data.frame()

  expect_true(all((h_ilm_coef$`Pr(>|t|)`) > .1))

  expect_true(h_ilm_coef$Estimate[1] == ilm$`Estimate_(Intercept)`)
  expect_true(h_ilm_coef$Estimate[5] == ilm$`Estimate_timeperiodend:currencybiomass`)

  expect_true(h_ilm_coef$`Pr(>|t|)`[1] == ilm$`Pr(>|t|)_(Intercept)`)
  expect_true(h_ilm_coef$`Pr(>|t|)`[5] == ilm$`Pr(>|t|)_timeperiodend:currencybiomass`)

  h_ilm_summary <- summary(h_ilm)

  expect_true(h_ilm_summary$r.squared == ilm$overall_r2)

  h_ilm_contrasts <- as.data.frame(pairs(emmeans::emmeans(h_ilm, ~ timeperiod | currency)))

  expect_true(h_ilm_contrasts$p.value[1] == ilm$abundance_contrastp.value)

})
