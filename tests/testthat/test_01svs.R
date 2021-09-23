test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

h <- BBSsize::hartland

h_isd <- BBSsize::simulate_isd_ts(h, isd_seed = 2021)


test_that("svs works", {

  h_svs <- get_annual_svs(h_isd$isd)

  expect_true(all(h_svs$year == 1994:2018))
  expect_true(nrow(h_svs) == 25)
  expect_true(ncol(h_svs) == 7)
  expect_false(anyNA(h_svs))
  expect_true(all(h_svs$abundance == rowSums(h$abundance)))
  expect_true((h_svs$biomass[1]) == sum((dplyr::filter(h_isd$isd, year == 1994)$mass)))
  expect_true((h_svs$energy[1]) == sum(BBSsize::estimate_b(dplyr::filter(h_isd$isd, year == 1994)$mass)))

  expect_true(all(floor(h_svs$mean_energy) == floor(h_svs$energy/h_svs$abundance)))
  expect_true(all(floor(h_svs$mean_biomass) == floor(h_svs$biomass/h_svs$abundance)))

  }
)

test_that("pull caps works no caps specified", {

  h_svs <- get_annual_svs(h_isd$isd)

  h_caps <- pull_caps(h_svs)

  expect_true(nrow(h_caps) == 10)

  expect_true(all(h_caps[1:5, 1:7] == h_svs[1:5, ]))
  expect_true(all(h_caps[6:10, 1:7] == h_svs[21:25, ]))
  expect_true(all(h_caps$timeperiod[1:5] == "begin"))
  expect_true(all(h_caps$timeperiod[6:10] == "end"))

  }
)


test_that("pull caps works with caps specified", {

  h_svs <- get_annual_svs(h_isd$isd)

  h_caps <- pull_caps(h_svs, begin_years = c(2000:2004), end_years = c(2010:2014))

  expect_true(nrow(h_caps) == 10)

  expect_true(all(h_caps[1:5, 1:7] == h_svs[7:11, ]))
  expect_true(all(h_caps[6:10, 1:7] == h_svs[17:21, ]))
  expect_true(all(h_caps$timeperiod[1:5] == "begin"))
  expect_true(all(h_caps$timeperiod[6:10] == "end"))

  expect_error(pull_caps(h_svs, begin_years = c(2000:2004), end_years = c(1994:1998)))

}
)


test_that("raw ratios works", {

  h_svs <- get_annual_svs(h_isd$isd)

  h_caps <- pull_caps(h_svs)

  raw_change <- compute_raw_sv_change(h_caps)

  expect_true(nrow(raw_change) == 1)
  expect_true(ncol(raw_change) == 5)

  expect_false(anyNA(raw_change))

  expect_true(raw_change$energy_raw_ratio == sum(h_caps$energy[6:10])/sum(h_caps$energy[1:5]))

  expect_true(raw_change$abundance_raw_ratio == sum(h_caps$abundance[6:10])/sum(h_caps$abundance[1:5]))

  expect_true(raw_change$biomass_raw_ratio == sum(h_caps$biomass[6:10])/sum(h_caps$biomass[1:5]))

  expect_true(raw_change$mean_energy_raw_ratio == (sum(h_caps$energy[6:10]) / sum(h_caps$abundance[6:10]))/(sum(h_caps$energy[1:5])/sum(h_caps$abundance[1:5])))

  expect_true(raw_change$mean_biomass_raw_ratio == (sum(h_caps$biomass[6:10]) / sum(h_caps$abundance[6:10]))/(sum(h_caps$biomass[1:5])/sum(h_caps$abundance[1:5])))

  }
)
