#' Get an ISD
#'
#' Tiny wrapper for `BBSsize::simulate_isd_ts` to get just the ISD component of that output.
#'
#' @param ts_comp MATSS-shaped (list with $abundance, $covariates, $metadata)
#' @param isd_seed can provide, to standardize results (e.g. for testing). if not provided, one will be drawn randomly and kept track of so you can replicate later.
#'
#' @return isd
#' @export
#'
#' @importFrom BBSsize simulate_isd_ts
just_isd <- function(ts_comp, isd_seed = NULL) {
  library(BBSsize)
  BBSsize::simulate_isd_ts(ts_comp, isd_seed = isd_seed)$isd
}
#' Construct GMMs for sampling
#'
#' Constructs well-resolved smooths for the begin and end ISDs to use for resampling. Repeatedly draws ISDs by sampling from species' intraspecific body size distributions and combines these repeated samples in to one large sample. Then fits a GMM to this sample to get a very smooth version of the GMM-ified ISD. This allows us to resample ISDs with noninteger ratios of species abundances in a consistent way.
#'
#' @param ts_comp matss-shaped dataset
#' @param n_isd_draws how many versions of the ISD to use to construct the density function. There isn't too much variation between them, defaults to 5. Adding more increases time nonlinear due to the GMM, do not recommend much more than 10.
#' @param initial_isd_seed don't recommend providing unless testing
#' @param begin_years defaults 1988-1992
#' @param end_years defaults 2014-2018
#'
#' @return list of two gmm densities as dataframes with columns mass, density, timeperiod
#' @importFrom dplyr filter bind_rows mutate
#' @export
construct_sampling_gmm <- function(ts_comp, n_isd_draws = 5, initial_isd_seed = NULL, begin_years = 1988:1992, end_years = 2014:2018) {

  # For debugging
  if(is.null(initial_isd_seed)) {
    set.seed(NULL)
    initial_isd_seed <- sample(1:1000000, 1)
  }

  set.seed(initial_isd_seed)

  isd_seeds <- sample.int(10000000, size = n_isd_draws)

  # Draw the ISD for the whole timeseries (draws masses for all individuals ever "seen") repeatedly.
  ts_isds_many <- lapply(isd_seeds, just_isd, ts_comp = ts_comp)
  names(ts_isds_many) <- 1:n_isd_draws

  # Combine over repeated draws to smooth out sampling variability
  ts_isds_many <- dplyr::bind_rows(ts_isds_many, .id = "sim")

  # Pull "begin" and "end" time periods
  ts_isds_begin <- dplyr::filter(ts_isds_many, year %in% begin_years)
  ts_isds_end <- dplyr::filter(ts_isds_many, year %in% end_years)

  # Fitting GMMs can be slightly sensitive to starting seed. There's no reason to be fancy about this, but for testing you want this to be fixed.
  set.seed(1977)

  # Fit the GMMS to the repeatedly-sampled ISDs. This will give dataframes with mass and density columns.
  ts_gmm_begin <- add_gmm(ts_isds_begin)
  ts_gmm_end <- add_gmm(ts_isds_end)

  set.seed(NULL)

  # Add timeperiod columns because you want to be sure about that!
  ts_gmm_begin <- dplyr::mutate(ts_gmm_begin, timeperiod = "begin")
  ts_gmm_end <- dplyr::mutate(ts_gmm_end, timeperiod = "end")

  return(list(begin = ts_gmm_begin, end = ts_gmm_end))

}


#' Draw individuals from a GMM
#'
#' Draw individual body masses from a GMM density function for the ISD
#'
#' @param nind how many to draw
#' @param sampling_gmm a dataframe with mass (on log scale) and density columns.
#' @param draw_seed don't recommmend providing
#'
#' @return drawn (log) body masses, and metainfo for whether the sampling gmm is from begin or end and the seed used when drawing
#' @export
#'
draw_individuals <- function(nind, sampling_gmm, draw_seed = NULL) {

  # The whole point of having draw_seed is for reproducibility at a kind of post-mortem level. I don't ever want to constrain this sampling process to come out the same way every time when running sims (e.g. to have the first 10 draws always be the same), but I do want to be able to set this for testing & debugging purposes.
  if(is.null(draw_seed)) {
    set.seed(NULL)
    draw_seed <- sample(1:1000000, 1)
  }

  set.seed(draw_seed)
  # The actual sampling
  drawn_individuals <- data.frame(
    logmass = sample(sampling_gmm$mass, nind, T, sampling_gmm$density),
    isd_timeperiod = sampling_gmm$timeperiod[1],
    sampling_seed = draw_seed
  )

  set.seed(NULL)
  drawn_individuals
}

#' Add drawn individuals
#'
#' Adds individuals drawn from a sampling GMM to an ISD-shaped dataframe.
#'
#' @param timeperiod_isd an ISD-shaped dataframe. what matters is the number of rows and that it have mass and id columns.
#' @param sampling_gmm begin or end, just one
#' @param draw_seed Don't recommend providing one.
#'
#' @return dataframe with body mass (and energy use) for those individuals drawn from the sampling_gmm provided
#' @export
#'
#' @importFrom dplyr select bind_cols mutate
#' @importFrom BBSsize estimate_b
add_drawn_individuals <- function(timeperiod_isd, sampling_gmm, draw_seed = NULL) {

  # Draw new body masses for the individuals in timeperiod_isd using the sampling gmm provided

  ntodraw = nrow(timeperiod_isd)

  drawn <- draw_individuals(ntodraw, sampling_gmm, draw_seed)

  # Swap out the masses in timeperiod_isd. Because the gmm is on log masses, before we get energy use and biomass we have to convert back to the raw values.
  timeperiod_isd <- timeperiod_isd %>%
    dplyr::select(-mass, -id) %>%
    dplyr::bind_cols(drawn) %>%
    dplyr::mutate(mass = exp(logmass)) %>%
    dplyr::mutate(energy = BBSsize::estimate_b(mass))

  timeperiod_isd
}

#' Draw communities
#'
#' High-level fxn to get yearly estimates for the state variables under the scenarios of 1) the ISD is conserved begin-end and 2) the ISD changes begin-end. A way to partition change in biomass or energy use expected given changes in abundance and sampling variation, and changes that must be partially driven by changes in the size shift.
#'
#' @param ts_comp matss-style dataset
#' @param begin_years default 1988:1992.
#' @param end_years default 2014:2018.
#' @param draw_seed don't provide unless you are doing something specific
#' @param sampling_gmms do provide for speed! result of construct_sampling_gmms
#' @param initial_isd_seed don't provide unless you have a specific reason to
#' @param raw_isd_seed don't provide unless you have a reason to
#'
#' @return df of sim inds
#' @export
#'
#' @importFrom BBSsize simulate_isd_ts
#' @importFrom dplyr filter mutate bind_rows group_by summarize ungroup n bind_cols
draw_communities <- function(ts_comp, begin_years = 1988:1992, end_years = 2014:2018, draw_seed = NULL, sampling_gmms = NULL, initial_isd_seed = NULL, raw_isd_seed = NULL) {

  # If GMM density smooths to sample from are not provided, get them
  # When running at scale, you want to provide them: 1) to ensure that all sims are being drawn from the same sampling GMMS, 2) for speed, it takes a long time to fit a GMM to a lot of ISD draws.

  if(is.null(sampling_gmms)) {

    sampling_gmms <- construct_sampling_gmm(ts_comp, begin_years = begin_years, end_years = end_years, initial_isd_seed = NULL)

  }

  # This doesn't actually matter, I am using the raw isds to get the correct shape data frame. I am keeping it this way, even though it's inefficient, because in the back of my mind I think I might want to pull out raw estimates for biomass and energy use (i.e. drawing individuals directly from normal distributions for each year, rather than drawing from the ISD combining over all years). It'll be a tiny change but I'm not 100% sure what I want to do with that info yet so holding off.

  if(is.null(raw_isd_seed)) {
    set.seed(NULL)
    raw_isd_seed <- sample(1:1000000, 1)
  }

  # Here I am sampling ISDs to get dfs of the correct shape to then sample new body masses from different density fxns.

  raw_isd <- BBSsize::simulate_isd_ts(ts_comp, isd_seed = raw_isd_seed)$isd

  begin_isd <- dplyr::filter(raw_isd, year %in% begin_years) %>%
    dplyr::mutate(timeperiod = "begin")
  end_isd <- dplyr::filter(raw_isd, year %in% end_years) %>%
    dplyr::mutate(timeperiod = "end")


  # Draw individuals for each time period from the MATCHING density functions
  # This will destroy interannual, intratimeperiod variation in the size structure, which we're OK with (the point of using 5-year intervals is to smooth out species accumulation)
  # I do not think it is a good idea to provide draw_seed, that will constrain things to come out the same in weird ways at this scale.
  begin_individuals <- add_drawn_individuals(begin_isd, sampling_gmms$begin, draw_seed = draw_seed)

  end_individuals <- add_drawn_individuals(end_isd, sampling_gmms$end, draw_seed = draw_seed)

  actual_individuals <- dplyr::bind_rows(begin_individuals, end_individuals) %>%
    dplyr::mutate(source = "actual")


  # Now draw individuals for each time period with scrambled ISDs. Specifically, draw for the beginning from the beginning ISD. But then also draw the end from the beginning ISD. This gives an "end" ISD pretending that the ISD didn't change from the beginning.
  begin_individuals_sim <- add_drawn_individuals(begin_isd, sampling_gmms$begin, draw_seed = draw_seed)

  end_individuals_sim <- add_drawn_individuals(end_isd, sampling_gmms$begin, draw_seed = draw_seed)


  sim_individuals <- dplyr::bind_rows(begin_individuals_sim, end_individuals_sim) %>%
    dplyr::mutate(source = "sim")


  # Go ahead and pull the raw state variable estimates too..

  raw_individuals <- dplyr::bind_rows(begin_isd, end_isd) %>%
    dplyr::mutate(energy = BBSsize::estimate_b(mass),
                  source = "raw",
                  isd_timeperiod = "raw",
                  sampling_seed = NA)


  all_individuals <- dplyr::bind_rows(actual_individuals, sim_individuals, raw_individuals)


  # Summarize individuals to get toal abundance, biomass, and energy use per year for each sim scenario.
  # And add route-level identifying info.
  all_svs <- all_individuals %>%
    dplyr::group_by(year, timeperiod, isd_timeperiod, sampling_seed, isd_seed, source) %>%
    dplyr::summarize(total_abundance = dplyr::n(),
                     total_biomass = sum(mass),
                     total_energy = sum(energy)) %>%
    dplyr::ungroup() %>%
    dplyr::bind_cols(as.data.frame(ts_comp$metadata$location)) %>%
    dplyr::mutate(matssname = paste0("bbs_rtrg_", route, "_", statenum))

  all_svs
}

#' Wrapper for draw_communities
#'
#' To facilitate pipelines. Just runs draw_communities `ndraws` times and tacks on a column keeping track of which sim.
#'
#' @param ts_comp matss-style dataset
#' @param begin_years default 1988:1992.
#' @param end_years default 2014:2018.
#' @param ndraws default 100
#' @param draw_seed don't provide unless you are doing something specific
#' @param sampling_gmms do provide for speed! result of construct_sampling_gmms
#' @param initial_isd_seed don't provide unless you have a specific reason to
#' @param raw_isd_seed don't provide unless you have a reason to
#'
#' @return sims
#' @export
#' @importFrom dplyr bind_rows
draw_communities_wrapper <- function(ts_comp, begin_years = 1988:1992, end_years = 2014:2018, ndraws = 100, draw_seed = NULL, sampling_gmms = NULL, initial_isd_seed = NULL, raw_isd_seed = NULL) {

  # Run draw communities ndraws times.
  drawn_communities <- replicate(ndraws, draw_communities(ts_comp, begin_years = begin_years, end_years = end_years, draw_seed = draw_seed, sampling_gmms = sampling_gmms, initial_isd_seed = initial_isd_seed, raw_isd_seed = raw_isd_seed), simplify = F)

  names(drawn_communities) <- 1:ndraws

  dplyr::bind_rows(drawn_communities, .id = "sim_iteration")

}
#' Fit brms to currencies
#'
#' @param some_sims dataframe with columns matssname, timeperiod, source, year, total_energy, total_biomass
#' @param cores how many cores to use. if on hpg, use ONE. if local, do what you want.
#' @param iter how many iterations. at scale, I've been using 4000 to be generous.
#'
#' @return list of brm fit on total_energy, on total_biomass, and dataset name
#' @export
#'
#' @importFrom brms brm
#' @importFrom dplyr filter
fit_brms <- function(some_sims, cores = 1, iter = 8000, thin =2) {

  # something in rwar as I currently have it is locking the namespace and interfering with drake, at least locally. this is not my best work but it gets rwar out of the namespace if it's attached.
  is_rwar_attached = any(grepl("rwar", names(sessionInfo()[7]$otherPkgs)))
  if(is_rwar_attached) {
    detach("package:rwar", unload = T)
  }

  # sims returns estimates of the raw values, which we don't want for the model fit (we jsut want the ones that come from drawing from the densityGMMS)
  justsims <- dplyr::filter(some_sims, source %in% c("actual", "sim")) # remove raw


  # Fit a brm on total_energy
  te_brm <- brms::brm(total_energy ~ (timeperiod * source) + (1 | year), data = justsims, cores = cores, iter = iter, thin = thin)

  # Fit the brm on total_biomass
  tb_brm <- brms::brm(total_biomass ~ (timeperiod * source) + (1 | year), data = justsims, cores = cores, iter = iter, thin = thin)

  # keep track of what dataset this is
  md <- some_sims$matssname[1]

  return(list(
    te_brm = te_brm,
    tb_brm = tb_brm,
    matssname =md
  ))

}



#' Extract estimates from brms posteriors
#'
#' Extract estimates for a bunch of quantities of interest from the posterior of fitted brms.
#'
#' @param some_brms list with brms te_brm, tb_brm, and matssname. Result of running `fit_brms()`
#'
#' @return dataframe with estimates. each row is a draw from the posterior. draws from the 2 models are stacked, group by `currency` to get each model.
#' @export
#'
extract_brm_ests <- function(some_brms){

  e_ests <- extract_ests(some_brms$te_brm, "energy", matssname = some_brms$matssname)
  b_ests <- extract_ests(some_brms$tb_brm, "biomass", matssname = some_brms$matssname)

  return(
    rbind(
      e_ests,
      b_ests
    )
  )

}


#' Extract estimates from one brm
#'
#' This is the workhorse function for extracting estimates from the posterior and using them to calculate various quantities.
#'
#' @param a_brm one of the currency brms
#' @param brm_currency which currency "energy" or "biomass"
#' @param matssname dataset name
#'
#' @return dataframe with estimates, each row is a draw from the posterior
#' @export
#'
#' @importFrom tidybayes tidy_draws
#' @importFrom dplyr mutate row_number filter rename select group_by_all mutate ungroup
extract_ests <- function(a_brm, brm_currency = NULL, matssname = NULL){


  # Get all draws from the posterior and get just the terms we want
  td <- tidybayes::tidy_draws(a_brm) %>%
    #  select_at(vars(starts_with("b"))) %>%
    dplyr::mutate(rowindex = dplyr::row_number()) #%>% # and get a row index to keep draws together, I'm not sure if this matters but I'll do it
  # dplyr::filter(rowindex > max(rowindex) / 2) # remove warmup

  td_ests <- td %>%
    dplyr::rename(timeperiodend_sourcesim = `b_timeperiodend:sourcesim`) %>%
    dplyr::select(rowindex, b_Intercept, b_timeperiodend, b_sourcesim, timeperiodend_sourcesim) %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(
      estimated_actual_begin = sum(b_Intercept), # estimated beginning value
      estimated_actual_end = sum(b_Intercept, b_timeperiodend), # estimated end value
      estimated_sim_begin = sum(b_Intercept, b_sourcesim), # estimated beginning value from sims. we expect this to be equal to the estimated beginning value, any change is just sampling error.
      estimated_sim_end = sum(b_Intercept, b_timeperiodend, timeperiodend_sourcesim, b_sourcesim),
      estimated_actual_change_ratio = (estimated_actual_end - estimated_actual_begin) / estimated_actual_begin, # this is a measure of the magnitude of the change from beginning to end. the sign is going to be increase (positive) or decrease. the magnitude is the % increase. so .1 = added 10% of starting (biomass or energy) to get to the end. -.2 = lost 20% of starting (biomass or energy) between begin and end.
      estimated_sim_change_ratio = (estimated_sim_end - estimated_sim_begin) / estimated_sim_begin, # same measure but having drawn the end values using the beginning isd. this is the amount of change expected due only to changes in the numbers of individuals observed in each time period. by comparing estimated_actual_change_ratio to estimated_sim_change_ratio, I believe we get an estimate of both the significance and magnitude of decoupling of (biomass or energy) and numerical abundance due to changes in the size spectrum.
      estimated_actual_change = estimated_actual_end - estimated_actual_begin, # the "slope" assuming x = 0 or 1 for begin or end. aka the absolute change from end to begin.
      estimated_sim_change = estimated_sim_end - estimated_sim_begin, # absolute change from end to begin due to abundance change
      estimated_change_ratio_deviation = estimated_actual_change_ratio - estimated_sim_change_ratio, # deviation of change ratios from 1:1
      estimated_change_deviation = estimated_actual_change - estimated_sim_change # deviation of actual change from 1:1

    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(currency = brm_currency,
                  matssname = matssname)


  return(td_ests)
}

#' Lower cutoff for a vector
#'
#' @param vector vector of values
#' @param lower_cutoff defaults 0.025
#'
#' @return cutoff value
#' @export
#'
lower_quantile <- function(vector, lower_cutoff = 0.025) {
  as.numeric(quantile(vector, probs = lower_cutoff))
}


#' Upper cutoff for a vector
#'
#' @param vector vector of values
#' @param upper_cutoff defaults 0.975
#'
#' @return cutoff value
#' @export
#'
upper_quantile <- function(vector, upper_cutoff = .975) {
  as.numeric(quantile(vector, probs = upper_cutoff))
}

#' Summarize brm ests
#'
#' Summarize draws from the posterior to give mean and upper/lower quantile estimates for quantities of interest (parameters, derived estimates)
#'
#' @param some_ests dataframe with columns matssname, currency, rowindex, and parameters/values of interest. result of extract_brm_ests
#'
#' @return dataframe summarized by matssname and currency to get mean and upper/lower qs.
#' @export
#'
#' @importFrom dplyr select group_by summarize_all ungroup
summarize_brm_ests <- function(some_ests) {

  td_route_ests_summary <- some_ests %>%
    dplyr::select(-rowindex) %>%
    dplyr::group_by(matssname, currency) %>%
    dplyr::summarize_all(.funs = list(mean = mean,
                                      lower = lower_quantile,
                                      upper = upper_quantile,
                                      median = median)) %>%
    dplyr::ungroup()

  td_route_ests_summary

}

