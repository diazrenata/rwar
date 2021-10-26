#' LOO compare on brms
#'
#' @param brms_fits for one currency
#'
#' @return a dataframe of results from loo_compare
#' @export
#'
#' @importFrom brms add_criterion loo_compare
#' @importFrom dplyr mutate row_number
compare_brms <- function(brms_fits) {

  brms_fits <- lapply(brms_fits, brms::add_criterion, criterion  = "loo")

  brms_comparison <- brms::loo_compare(brms_fits[[1]], brms_fits[[2]],brms_fits[[3]], model_names = names(brms_fits)) %>%
    as.data.frame() %>%
    dplyr::mutate(model = row.names(.),
                  rank = dplyr::row_number())

  return(brms_comparison)

}

#' Pull comparisons for two sets of brms
#'
#' @param some_brms_fits a list
#'
#' @return a dataframe
#' @export
#'
#' @importFrom dplyr bind_rows mutate
compare_both_brms <- function(some_brms_fits) {

  biomass <- compare_brms(some_brms_fits$tb_brms)
  energy <- compare_brms(some_brms_fits$te_brms)

  both_comparisons <- dplyr::bind_rows(biomass = biomass, energy = energy, .id = "currency") %>%
    dplyr::mutate(matssname = some_brms_fits$matssname,
                  simtype = some_brms_fits$simtype[1])


  return(both_comparisons)

}

#' Fit multiple candidate brms
#'
#' @param some_sims dataframe with columns matssname, timeperiod, source, total_energy, total_biomass
#' @param cores how many cores to use. if on hpg, use ONE. if local, do what you want.
#' @param iter how many iterations. at scale, I've been using 4000 to be generous.
#'
#' @return list of brm fit on total_energy, on total_biomass, and dataset name
#' @export
#'
#' @importFrom brms brm prior
#' @importFrom dplyr filter
fit_brms <- function(some_sims, cores = 1, iter = 8000, thin =2) {


  # Fit a brm on total_energy
  te_brm_full <- brms::brm(total_energy ~ (timeperiod * source) , data = some_sims, cores = cores, iter = iter, thin = thin)
  te_brm_nosource <- brms::brm(total_energy ~ (timeperiod), data = some_sims, cores = cores, iter = iter, thin = thin)
  te_brm_notime <- brms::brm(total_energy ~ 1, data = some_sims, cores = cores, iter = iter, thin = thin)

  te_brms = list(
    te_brm_full = te_brm_full,
    te_brm_nosource = te_brm_nosource,
    te_brm_notime = te_brm_notime
  )


  # Fit the brm on total_biomass
  tb_brm_full <- brms::brm(total_biomass ~ (timeperiod * source) , data = some_sims, cores = cores, iter = iter, thin = thin)
  tb_brm_nosource <- brms::brm(total_biomass ~ (timeperiod) , data = some_sims, cores = cores, iter = iter, thin = thin)
  tb_brm_notime <- brms::brm(total_biomass ~ 1 , data = some_sims, cores = cores, iter = iter, thin = thin)


  tb_brms = list(
    tb_brm_full = tb_brm_full,
    tb_brm_nosource = tb_brm_nosource,
    tb_brm_notime = tb_brm_notime
  )

  return(list(
    te_brms = te_brms,
    tb_brms = tb_brms,
    matssname =some_sims$matssname[1],
    simtype = some_sims$simtype[1]
  ))

}


#' Select the simplest model within 1 se of the best model
#'
#' @param some_compares from compare_both_brms
#'
#' @return df with model names for best model for each currency
#' @export
#' @importFrom dplyr mutate arrange filter group_by ungroup row_number
loo_select <- function(some_compares) {


  winners <- some_compares %>%
    dplyr::mutate(model_complexity = ifelse(grepl("full", model), 3,
                                     ifelse(grepl("source", model), 2, 1))) %>%
    dplyr::arrange(matssname, currency, simtype, rank) %>%
    dplyr::mutate(in_one_se = (elpd_diff + se_diff ) >= 0) %>%
    dplyr::filter(in_one_se) %>%
    dplyr::group_by(currency) %>%
    dplyr::arrange(model_complexity) %>%
    dplyr::mutate(model_rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(model_rank == 1)

  return(winners)
}

#' Extract draws from winning models
#'
#' @param some_winners df with rows for energy and biomass, column "model" name of best model
#' @param some_models list of models from fit_brms
#'
#' @return df of draws from both models
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom tidybayes tidy_draws
winner_draws <- function(some_winners, some_models) {

  winner_energy_mod <- some_winners %>%
    dplyr::filter(currency == "energy") %>%
    dplyr::select(model) %>%
    as.character()

  winner_biomass_mod <- some_winners %>%
    dplyr::filter(currency == "biomass") %>%
    dplyr::select(model) %>%
    as.character()

  te_draws <- tidybayes::tidy_draws(some_models$te_brms[[winner_energy_mod]]) %>%
    dplyr::mutate(currency = "energy", modtype = winner_energy_mod)

  tb_draws <- tidybayes::tidy_draws(some_models$tb_brms[[winner_biomass_mod]]) %>%
    dplyr::mutate(currency = "biomass", modtype = winner_biomass_mod)

  all_draws <- dplyr::bind_rows(te_draws, tb_draws) %>%
    dplyr::mutate(matssname = some_winners$matssname[1],
                  simtype = some_winners$simtype[1])

  return(all_draws)
}


#' Extract qis from draws
#'
#' @param some_draws df of draws
#'
#' @return df of qis
#' @export
#'
#' @importFrom tidybayes median_qi
#' @importFrom dplyr group_by ungroup
winner_qis <- function(some_draws) {

  some_qis <- some_draws %>%
    dplyr::group_by(currency, modtype, matssname, simtype) %>%
    tidybayes::median_qi(.width = c(.95, .99)) %>%
    dplyr::ungroup()


  return(some_qis)
}
