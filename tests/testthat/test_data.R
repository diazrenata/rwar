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
  expect_true(floor(h_svs$energy[1]) == 94302)
  expect_true(floor(h_svs$biomass[1]) == 36983)
  expect_true(all(floor(h_svs$mean_energy) == floor(h_svs$energy/h_svs$abundance)))
  expect_true(all(floor(h_svs$mean_biomass) == floor(h_svs$biomass/h_svs$abundance)))

  }
)

