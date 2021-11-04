#' Perform species shuffling for local null model
#'
#' @param ts_dat matss-style
#' @param null_mod_seed or null
#'
#' @return matss-style, but shuffled
#' @export
#'
shuffle_local <- function(ts_dat, null_mod_seed = NULL) {


  orig_species <- colnames(ts_dat$abundance)

  set.seed(null_mod_seed)

  shuffled_species <- sample(orig_species, size = length(orig_species), replace = F)

  shuffled_dat <- ts_dat

  colnames(shuffled_dat$abundance) <- shuffled_species

  return(shuffled_dat)
}

#' Local null model
#'
#' Reassign species IDs using all species ever observed on a single route.
#'
#' @param ts_dat matssdat
#' @param null_mod_seed seed to control reshuffling of species
#' @param sim_index for wrapping
#' @param begin_years default 1988:1992
#' @param end_years default 2014:2018
#' @param isd_seed default 1989, passed to ssims_wrapper
#' @param n_isd_draws default 10, passed to ssims_wrapper
#' @param ndraws default 100, passed to ssims_wraper
#' @param return_actual return actual (unmodified) results?
#'
#' @return df of sims
#' @export
#' @importFrom dplyr left_join
local_null_model <- function(ts_dat, null_mod_seed = 1989, sim_index = 1, begin_years = 1988:1992, end_years = 2014:2018, isd_seed = 1989, n_isd_draws = 10, ndraws = 100, return_actual = FALSE) {

  if(is.null(null_mod_seed)) {
    set.seed(NULL)
    null_mod_seed <- sample.int(1000000000, 1)
  }


  shuffled_dat <- shuffle_local(ts_dat = ts_dat, null_mod_seed = null_mod_seed)


  if(return_actual) {

    shuffled_dat <- ts_dat
    sim_index <- -99
  }

  null_sims <- ssims_wrapper(shuffled_dat, initial_isd_seed_gmm =isd_seed, initial_draw_seed = null_mod_seed, n_isd_draws = n_isd_draws, ndraws = ndraws, simtype = "actual")

  total_abund <- data.frame(
    year = ts_dat$covariates$year,
    total_abundance = rowSums(ts_dat$abundance)
  ) %>%
    dplyr::filter(year %in% c(1988:1992, 2014:2018))

  null_sims <- dplyr::left_join(null_sims, total_abund)

  results <- null_sims %>%
    dplyr::mutate(
      null_mod_type = "local",
      null_mod_seed = null_mod_seed,
      sim_index = sim_index
    )


  results
}

#' Wrapper for local null model
#'
#' @param ts_dat ts dat
#' @param n_null_model_sims number of sims
#' @param begin_years defaults good
#' @param end_years defaults good
#' @param isd_seed default good
#' @param n_isd_draws default 10
#' @param ndraws default 100
#' @param initial_null_model_seed default 1989
#' @param return_actual return actual?
#'
#' @return df of many nms
#' @export
#' @importFrom dplyr bind_rows
local_null_model_wrapper <- function(ts_dat, n_null_model_sims = 100, begin_years = 1988:1992, end_years = 2014:2018, isd_seed = 1989, n_isd_draws = 10, ndraws = 100, initial_null_model_seed = 1989, return_actual = F) {

  null_model_seeds <- initial_null_model_seed:(initial_null_model_seed + n_null_model_sims - 1)

  repeated_nulls <- list()

  for(i in 1:n_null_model_sims) {

    repeated_nulls[[i]] <- local_null_model(ts_dat, null_mod_seed = null_model_seeds[i], sim_index = i, begin_years = begin_years, end_years = end_years, isd_seed = isd_seed, ndraws = ndraws, n_isd_draws = n_isd_draws, return_actual = return_actual)
  }

  results <- dplyr::bind_rows(repeated_nulls)

  results

}

#' Return ratio of mean e/b begin:end
#'
#' @param nm_results result of null_model_wrapper
#'
#' @return df of summary results
#' @export
#' @importFrom dplyr mutate group_by summarize ungroup
#' @importFrom tidyr pivot_longer pivot_wider
summarize_null_models <- function(nm_results) {
  nm_summary <- nm_results %>%
    dplyr::mutate(mean_energy = total_energy / total_abundance,
                  mean_biomass = total_biomass / total_abundance) %>%
    dplyr::group_by(matssname, simtype,null_mod_type,  null_mod_seed, timeperiod, source) %>%
    dplyr::summarize(mean_energy = mean(mean_energy),
              mean_biomass = mean(mean_biomass)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(-c(matssname, simtype, null_mod_type, null_mod_seed, timeperiod, source), names_to = "e_or_b", values_to = "avg") %>%
    tidyr::pivot_wider(id_cols = c(matssname, simtype,null_mod_type,  null_mod_seed, e_or_b, timeperiod), names_from = source, values_from = avg) %>%
    dplyr::mutate(ratio = currency / abundance) %>%
    tidyr::pivot_wider(id_cols = c(matssname, simtype, null_mod_type, null_mod_seed, e_or_b), names_from = timeperiod, values_from = ratio) %>%
    dplyr::mutate(time_ratio = end / begin)

  nm_summary
}


#' Compare null models to reality
#'
#' @param actual_summary actual
#' @param null_summary null
#'
#' @return df
#' @export
#'
#' @importFrom dplyr select rename left_join group_by summarize distinct
compare_actual_null_models <- function(actual_summary, null_summary) {

  actual_vals <- actual_summary %>%
    dplyr::select(e_or_b, time_ratio, matssname) %>%
    dplyr::rename(time_ratio_actual = time_ratio)

  null_scores <- null_summary %>%
    dplyr::left_join(actual_vals) %>%
    dplyr::group_by(matssname, e_or_b) %>%
    dplyr::summarize(ses = ses(time_ratio_actual, time_ratio),
                     percentile = percentile_score(time_ratio_actual, time_ratio)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  null_scores

}

