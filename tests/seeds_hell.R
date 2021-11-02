library(rwar)
library(BBSsize)
library(dplyr)
dat = granby

# This should always reproduce with the same seed and not with a different seed or no seed
dat_gmms <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
# dat_gmms2 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 1989)
# dat_gmms3 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2, initial_isd_seed = 22)
# dat_gmms4 <- rwar::construct_sampling_gmm(dat, n_isd_draws = 2)
#
# all_equal(dat_gmms$begin, dat_gmms2$begin) # expect yes
# all_equal(dat_gmms$end, dat_gmms2$end) # expect yes
#
# all(dat_gmms$begin$density == dat_gmms3$begin$density) # expect no
# all(dat_gmms$end$density == dat_gmms3$end$density) # expect no
# all(dat_gmms$end$mass == dat_gmms3$end$mass) # expect yes
#
#
# all(dat_gmms4$begin$density == dat_gmms3$begin$density) # expect no
# all(dat_gmms4$end$density == dat_gmms3$end$density) # expect no
# all(dat_gmms4$end$mass == dat_gmms3$end$mass) # expect yes

as <- make_actual_sims(granby, initial_isd_seed_gmm = 1989, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 22, raw_isd_seed = 13, initial_draw_seed = 27)

as2 <- make_actual_sims(granby, initial_isd_seed_gmm = 1989, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 22, raw_isd_seed = 13, initial_draw_seed = 27)

all_equal(as, as2) #expect yes

# different sims should have different values

as_sim1 <- filter(as, sim_iteration == "1")
as_sim2 <- filter(as, sim_iteration == "2")

as_sim1$total_biomass == as_sim2$total_biomass # expect no for sims and yes for raws


# change initial_isd_seed_sim - should not matter
as3 <- make_actual_sims(granby, initial_isd_seed_gmm = 1989, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 15, raw_isd_seed = 13, initial_draw_seed = 27)

all_equal(as, as3) #expect yes


# change raw_isd_seed - should change "raw" values but nothing else

as4 <- make_actual_sims(granby, initial_isd_seed_gmm = 1989, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 15, raw_isd_seed = 19, initial_draw_seed = 27)

as_sims <- filter(as, source != "raw")
as4_sims <- filter(as4, source != "raw")

all_equal(as_sims, as4_sims)


as_raw <- filter(as, source == "raw")
as4_raw <- filter(as4, source == "raw")

all(as_raw$total_biomass == as4_raw$total_biomass) # expect false


# change initial_draw_seed - should change "sim" but not raw values

as5 <- make_actual_sims(granby, initial_isd_seed_gmm = 1989, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 15, raw_isd_seed = 19, initial_draw_seed = 25)

all(as_sims$total_biomass == as5_sims$total_biomass) # expect false
all(as_sims$total_abundance == as5_sims$total_abundance) # expect true

as5_raw <- filter(as5, source == "raw")

all_equal(as5_raw, as4_raw) # expect true

as5_sim1 <- filter(as5, sim_iteration == "1")
as5_sim2 <- filter(as5, sim_iteration == "2")

as5_sim1$total_biomass == as5_sim2$total_biomass # expect no for sims and yes for raws

# change gmm seed - should change "sim" but not "raw"


as6 <- make_actual_sims(granby, initial_isd_seed_gmm = 2021, n_isd_draws = 2, ndraws = 2, initial_isd_seed_sim = 22, raw_isd_seed = 13, initial_draw_seed = 27)

as6$total_biomass == as$total_biomass # expect true for raw and false for sims



