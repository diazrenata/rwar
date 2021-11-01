test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

g <- BBSsize::granby

test_that("just_isd", {

  granby_isd <- just_isd(g, isd_seed = 1989)

  expect_true(nrow(granby_isd) == sum(g$abundance))
  expect_true(unique(granby_isd$isd_seed) == 1989)
  expect_true(all(granby_isd$mass > 0))

  granby_isd2 <- just_isd(g, isd_seed = 1977)
  granby_isd3 <- just_isd(g, isd_seed = 1977)

  expect_false(all(granby_isd$mass == granby_isd2$mass))

  expect_true(all(granby_isd2 == granby_isd3))

})

test_that("sampling_gmms", {

  sgmms <- construct_sampling_gmm(g, initial_isd_seed = 1989)
  sgmms2 <- construct_sampling_gmm(g, initial_isd_seed = 1977)

  sgmms3 <- construct_sampling_gmm(g, initial_isd_seed = 1977)

  sgmms4 <-  construct_sampling_gmm(g)

  expect_true(sum(sgmms$begin$density) == 1)
  expect_true(sum(sgmms$end$density) == 1)

  expect_true(all(sgmms$begin$mass == sgmms2$begin$mass))
  expect_false(all(sgmms$begin$density == sgmms2$begin$density))

  expect_true(all(sgmms2$begin == sgmms3$begin))

  expect_false(all(sgmms2$begin$density == sgmms4$begin$density))

})

test_that("draw_communities", {

  sgmms <- construct_sampling_gmm(g, initial_isd_seed = 1989)

  draws1 <- draw_communities_wrapper(g, ndraws = 2, draw_seed = 1989, sampling_gmms = sgmms, raw_isd_seed = 1977)

  draws2 <- draw_communities_wrapper(g, ndraws = 2, draw_seed = 1989, sampling_gmms = sgmms, raw_isd_seed = 1977)


  draws3 <- draw_communities_wrapper(g, ndraws = 2, sampling_gmms = sgmms)

  draws4 <- draw_communities_wrapper(g, ndraws = 2, sampling_gmms = sgmms)

  draws5 <- draw_communities_wrapper(g, ndraws = 2)
  draws6 <- draw_communities_wrapper(g, ndraws = 2)

  expect_true(dplyr::all_equal(draws1, draws2))

  expect_false(all(draws1$total_biomass == draws3$total_biomass))

  expect_false(all(draws3$total_biomass == draws4$total_biomass))

  expect_false(all(draws5$total_biomass == draws1$total_biomass))

  expect_false(all(draws5$total_biomass == draws6$total_biomass))


  expect_true(all(sort(unique(draws3$isd_timeperiod)) == c("begin", "end", "raw")))

  expect_true(all(sort(unique(draws3$timeperiod)) == c("begin", "end")))

  expect_true(all(unique(draws3$year) == c(1988:1992, 2014:2018)))

  abund <- dplyr::filter(draws1, source == "abundance")

  expect_true(unique(abund$isd_timeperiod) == "begin")

  expect_true(all(unique(abund$timeperiod) == c("begin", "end")))

  curr <- dplyr::filter(draws1, source == "currency")

  expect_true(all(unique(curr$isd_timeperiod) == c("begin", "end")))


  expect_true(all(unique(curr$timeperiod) == c("begin", "end")))
  raw <- dplyr::filter(draws1, source == "raw")

  expect_true(all(unique(raw$isd_timeperiod) == "raw"))


  expect_true(all(unique(raw$timeperiod) == c("begin", "end")))

})
