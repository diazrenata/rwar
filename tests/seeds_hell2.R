library(rwar)
library(BBSsize)
library(dplyr)
dat = granby

# This should always reproduce with the same seed and not with a different seed. The default seed is 1989.
dat_gmms <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
dat_gmms2 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
dat_gmms3 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 22)
dat_gmms4 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2)

all_equal(dat_gmms$begin, dat_gmms2$begin) # expect yes
all_equal(dat_gmms$end, dat_gmms2$end) # expect yes

all(dat_gmms$begin$density == dat_gmms3$begin$density) # expect no
all(dat_gmms$end$density == dat_gmms3$end$density) # expect no
all(dat_gmms$end$mass == dat_gmms3$end$mass) # expect yes


all(dat_gmms4$begin$density == dat_gmms3$begin$density) # expect no
all(dat_gmms4$end$density == dat_gmms3$end$density) # expect no
all(dat_gmms4$end$mass == dat_gmms3$end$mass) # expect yes

draws <- draw_communities(granby, sampling_gmms  =dat_gmms)

multi_draws <- draw_communities_wrapper(granby, ndraws = 5, sampling_gmms = dat_gmms)

multi_draws_summ <- multi_draws %>%
  group_by(source,
           timeperiod,
           isd_timeperiod) %>%
  summarize(n_seeds = length(unique(sampling_seed)),
            seeds = toString(unique(sampling_seed)))

all(multi_draws_summ$n_seeds == 5)
length(unique(multi_draws$sampling_seed) == 20)

multi_actual_draws <- make_actual_sims(granby,n_isd_draws = 2, ndraws= 5)
all(multi_actual_draws[,1:22] == multi_draws)

multi_nc_draws <- make_nochange_sims(granby,n_isd_draws = 2, ndraws= 5)
multi_nc_draws2 <- make_nochange_sims(granby,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)

multi_draws_summ <- multi_nc_draws %>%
  group_by(source,
           timeperiod,
           isd_timeperiod) %>%
  summarize(n_seeds = length(unique(sampling_seed)),
            seeds = toString(unique(sampling_seed)))

all(multi_draws_summ$n_seeds == 5)
length(unique(multi_draws$sampling_seed) == 20)


multi_nsc_draws <- make_nosizechange_sims(granby, n_isd_draws = 2, ndraws = 5)
multi_nsc_draws2 <- make_nosizechange_sims(granby,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)
multi_nsc_draws3 <- make_nosizechange_sims(granby, n_isd_draws = 2, ndraws = 5, initial_draw_seed = 1989)

multi_nsc_draws$total_biomass == multi_nsc_draws2$total_biomass #xpect no
all_equal(multi_nsc_draws, multi_nsc_draws3) # expect y

multi_draws_summ <- multi_nsc_draws %>%
  group_by(source,
           timeperiod,
           isd_timeperiod) %>%
  summarize(n_seeds = length(unique(sampling_seed)),
            seeds = toString(unique(sampling_seed)))

all(multi_draws_summ$n_seeds == 5)
length(unique(multi_draws$sampling_seed) == 20)
multi_nc_draws2 <- make_nochange_sims(granby,n_isd_draws = 2, ndraws= 5, initial_draw_seed = 2000)
multi_nc_draws3 <- make_nochange_sims(granby, n_isd_draws = 2, ndraws = 5, initial_draw_seed = 1989)

multi_nc_draws$total_biomass == multi_nc_draws2$total_biomass #xpect no
all_equal(multi_nc_draws, multi_nc_draws3) # expect y


library(ggplot2)


ggplot(multi_actual_draws, aes(year, total_biomass, color  = source)) + geom_point()
ggplot(multi_nc_draws, aes(year, total_biomass, color  = source)) + geom_point()
ggplot(multi_nsc_draws, aes(year, total_biomass, color  = source)) + geom_point()


actual_ssims <- ssims_wrapper(granby, "actual", n_isd_draws = 2, ndraws = 5)
nc_ssims <- ssims_wrapper(granby, "nc", n_isd_draws = 2, ndraws = 5)
nsc_ssims <- ssims_wrapper(granby, "nsc", n_isd_draws = 2, ndraws = 5)


# expect all true
actual_summ <- summarize_sims(multi_actual_draws)
nc_summ <- summarize_sims(multi_nc_draws)
nsc_summ <- summarize_sims(multi_nsc_draws)

all_equal(actual_ssims, actual_summ)
all_equal(nc_ssims, nc_summ)
all_equal(nsc_ssims, nsc_summ)
