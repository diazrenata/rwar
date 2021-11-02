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

  dat= granby

  dat_gmms <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
  dat_gmms2 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
  dat_gmms3 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 22)
  dat_gmms4 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2)

  expect_true(dplyr::all_equal(dat_gmms$begin, dat_gmms2$begin)) # expect yes
  expect_true(dplyr::all_equal(dat_gmms$end, dat_gmms2$end)) # expect yes

  expect_false(all(dat_gmms$begin$density == dat_gmms3$begin$density)) # expect no
  expect_false(all(dat_gmms$end$density == dat_gmms3$end$density)) # expect no
  expect_true(all(dat_gmms$end$mass == dat_gmms3$end$mass)) # expect yes


  expect_false(all(dat_gmms4$begin$density == dat_gmms3$begin$density)) # expect no
  expect_false(all(dat_gmms4$end$density == dat_gmms3$end$density)) # expect no
  expect_true(all(dat_gmms4$end$mass == dat_gmms3$end$mass)) # expect yes

})

test_that("draw_communities", {

  sgmms <- construct_sampling_gmm(g, initial_isd_seed = 1989)

  draws1 <- draw_communities_wrapper(g, ndraws = 2, initial_draw_seed = 1989, sampling_gmms = sgmms)
  draws2 <- draw_communities_wrapper(g, ndraws = 2, initial_draw_seed = 1989, sampling_gmms = sgmms)
  draws3 <- draw_communities_wrapper(g, ndraws = 2, initial_draw_seed = 1994, sampling_gmms = sgmms)


  expect_true(dplyr::all_equal(draws1, draws2))

  expect_false(all(draws1$total_biomass == draws3$total_biomass))

  expect_true(all(sort(unique(draws3$isd_timeperiod)) == c("begin", "end")))

  expect_true(all(sort(unique(draws3$timeperiod)) == c("begin", "end")))

  expect_true(all(unique(draws3$year) == c(1988:1992, 2014:2018)))

  abund <- dplyr::filter(draws1, source == "abundance")

  expect_true(unique(abund$isd_timeperiod) == "begin")

  expect_true(all(unique(abund$timeperiod) == c("begin", "end")))

  curr <- dplyr::filter(draws1, source == "currency")

  expect_true(all(unique(curr$isd_timeperiod) == c("begin", "end")))

  expect_true(all(unique(curr$timeperiod) == c("begin", "end")))

  multi_draws <- draw_communities_wrapper(g, ndraws = 5, sampling_gmms = sgmms)

  multi_draws_summ <- multi_draws %>%
    dplyr::group_by(source,
             timeperiod,
             isd_timeperiod) %>%
    dplyr::summarize(n_seeds = length(unique(sampling_seed)),
              seeds = toString(unique(sampling_seed)))

  expect_true(all(multi_draws_summ$n_seeds == 5))
  expect_true(length(unique(multi_draws$sampling_seed) == 20))


})



test_that("sims wrappers", {

  sgmms <- construct_sampling_gmm(g, initial_isd_seed = 1989)

  multi_draws <- draw_communities_wrapper(g, ndraws = 5, sampling_gmms = sgmms, initial_draw_seed = 1989)

  multi_actual_draws <- make_actual_sims(g, ndraws= 5, initial_draw_seed = 1989, initial_isd_seed_gmm = 1989)

  expect_true(all(multi_actual_draws[,1:22] == multi_draws))


  multi_actual_draws2 <- make_actual_sims(g, ndraws= 5, initial_draw_seed = 2000, initial_isd_seed_gmm = 1989)

  expect_false(all(multi_actual_draws2$total_biomass == multi_actual_draws$total_biomass))

  multi_nc_draws <- make_nochange_sims(g,n_isd_draws = 2, ndraws= 5)
  multi_nc_draws2 <- make_nochange_sims(g,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)

  multi_draws_summ <- multi_nc_draws %>%
    dplyr::group_by(source,
             timeperiod,
             isd_timeperiod) %>%
    dplyr::summarize(n_seeds = length(unique(sampling_seed)),
              seeds = toString(unique(sampling_seed)))

  expect_true(all(multi_draws_summ$n_seeds == 5))
  expect_true(length(unique(multi_draws$sampling_seed) == 20))


  multi_nsc_draws <- make_nosizechange_sims(granby, n_isd_draws = 2, ndraws = 5)
  multi_nsc_draws2 <- make_nosizechange_sims(granby,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)
  multi_nsc_draws3 <- make_nosizechange_sims(granby, n_isd_draws = 2, ndraws = 5, initial_draw_seed = 1989)

  expect_false(all(multi_nsc_draws$total_biomass == multi_nsc_draws2$total_biomass)) #xpect no
  expect_true(dplyr::all_equal(multi_nsc_draws, multi_nsc_draws3)) # expect y

  multi_draws_summ <- multi_nsc_draws %>%
    dplyr::group_by(source,
             timeperiod,
             isd_timeperiod) %>%
    dplyr::summarize(n_seeds = length(unique(sampling_seed)),
              seeds = toString(unique(sampling_seed)))

  expect_true(all(multi_draws_summ$n_seeds == 5))
  expect_true(length(unique(multi_draws$sampling_seed) == 20))
  multi_nc_draws2 <- make_nochange_sims(granby,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)
  multi_nc_draws3 <- make_nochange_sims(granby, n_isd_draws = 2, ndraws = 5, initial_draw_seed = 1989)

  expect_false(all(multi_nc_draws$total_biomass == multi_nc_draws2$total_biomass)) #xpect no
  expect_true(dplyr::all_equal(multi_nc_draws, multi_nc_draws3)) # expect y

  actual_ssims <- ssims_wrapper(granby, "actual", n_isd_draws = 2, ndraws = 5)
  nc_ssims <- ssims_wrapper(granby, "nc", n_isd_draws = 2, ndraws = 5)
  nsc_ssims <- ssims_wrapper(granby, "nsc", n_isd_draws = 2, ndraws = 5)


  # expect all true
  actual_summ <- summarize_sims(multi_actual_draws)
  nc_summ <- summarize_sims(multi_nc_draws)
  nsc_summ <- summarize_sims(multi_nsc_draws)

  expect_true(dplyr::all_equal(actual_ssims, actual_summ))
  expect_true(dplyr::all_equal(nc_ssims, nc_summ))
  expect_true(dplyr::all_equal(nsc_ssims, nsc_summ))



})


