#' Select the simplest model within 1 se of the best model
#'
#' @param some_compares from compare_both_stanarms
#'
#' @return df with model names for best model for each currency
#' @export
#' @importFrom dplyr mutate arrange filter group_by ungroup row_number
loo_select <- function(some_compares) {


  winners <- some_compares %>%
    dplyr::mutate(model_complexity = ifelse(grepl("full", model), 3,
                                     ifelse(grepl("source", model), 2, 1))) %>%
    dplyr::arrange(matssname, currency, simtype, rank) %>%
    dplyr::group_by(matssname, currency, simtype) %>%
    # dplyr::mutate(best_elpd = elpd_loo[1],
    #               best_se = se_elpd_loo[1]) %>%
    # dplyr::ungroup() %>%
    # dplyr::mutate(one_se_cutoff = best_elpd - best_se) %>%
    dplyr::mutate(in_one_se = (elpd_diff + se_diff) >= 0) %>%
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
#' @param some_models list of models from fit_stanlm
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

  te_draws <- tidybayes::tidy_draws(some_models$te_stanlms[[winner_energy_mod]]) %>%
    dplyr::mutate(currency = "energy", modtype = winner_energy_mod)

  tb_draws <- tidybayes::tidy_draws(some_models$tb_stanlms[[winner_biomass_mod]]) %>%
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


#' Extract model diagnostics
#'
#' for many unsupervised stanlms
#'
#' @param some_fits some_fits
#'
#' @return data frame of rhats, neffs, divergent
#' @export
#'
#' @importFrom brms nuts_params rhat neff_ratio
#' @importFrom dplyr mutate filter group_by ungroup bind_rows summarize bind_cols
extract_diagnostics <- function(some_fits) {

  te_diagnostics <- list()

  for(mod_name in names(some_fits$te_stanlms)) {
    nuts <- brms::nuts_params(some_fits$te_stanlms[[mod_name]])
    rhats <- brms::rhat(some_fits$te_stanlms[[mod_name]]) %>%
      t() %>%
      as.data.frame()
    colnames(rhats) <- paste0("rhat_", colnames(rhats))

    neffs <- brms::neff_ratio(some_fits$te_stanlms[[mod_name]]) %>%
      t() %>%
      as.data.frame()
    colnames(neffs) <- paste0("neff_", colnames(neffs))

    divergents <- nuts %>%
      dplyr::filter(Parameter == "divergent__") %>%
      dplyr::summarize(divergent_sum = sum(Value)) %>%
      dplyr::ungroup()

    te_diagnostics[[mod_name]] <- dplyr::bind_cols(neffs, rhats) %>%
      dplyr::mutate(divergent_sum = divergents$divergent_sum,
                    model = mod_name,
                    currency = "energy")
  }
  tb_diagnostics <- list()

  for(mod_name in names(some_fits$tb_stanlms)) {
    nuts <- brms::nuts_params(some_fits$tb_stanlms[[mod_name]])
    rhats <- brms::rhat(some_fits$tb_stanlms[[mod_name]]) %>%
      t() %>%
      as.data.frame()
    colnames(rhats) <- paste0("rhat_", colnames(rhats))

    neffs <- brms::neff_ratio(some_fits$tb_stanlms[[mod_name]]) %>%
      t() %>%
      as.data.frame()
    colnames(neffs) <- paste0("neff_", colnames(neffs))

    divergents <- nuts %>%
      dplyr::filter(Parameter == "divergent__") %>%
      dplyr::summarize(divergent_sum = sum(Value)) %>%
      dplyr::ungroup()

    tb_diagnostics[[mod_name]] <- dplyr::bind_cols(neffs, rhats) %>%
      dplyr::mutate(divergent_sum = divergents$divergent_sum,
                    model = mod_name,
                    currency = "biomass")
  }

  all_diagnostics <- dplyr::bind_rows(
    dplyr::bind_rows(te_diagnostics),
    dplyr::bind_rows(tb_diagnostics)) %>%
    dplyr::mutate(matssname = some_fits$matssname,
                  simtype = some_fits$simtype)

  return(all_diagnostics)
}


#' Fit using rstanarm
#'
#' @param some_sims sims
#'
#' @return list of fits
#' @export
#'
#' @importFrom rstanarm stan_glm
fit_stanlm <- function(some_sims) {

  # Fit a stanlm on total_energy
  te_stanlm_full <- rstanarm::stan_glm(total_energy ~ (timeperiod * source) , data = some_sims, iter =8000, thin = 4)
  te_stanlm_nosource <- rstanarm::stan_glm(total_energy ~ (timeperiod), data = some_sims, iter =8000, thin = 4)
  te_stanlm_notime <- rstanarm::stan_glm(total_energy ~ 1, data = some_sims, iter =8000, thin = 4)

  te_stanlms = list(
    te_stanlm_full = te_stanlm_full,
    te_stanlm_nosource = te_stanlm_nosource,
    te_stanlm_notime = te_stanlm_notime
  )


  # Fit the stanlm on total_biomass
  tb_stanlm_full <- rstanarm::stan_glm(total_biomass ~ (timeperiod * source) , data = some_sims, iter = 8000, thin = 4)
  tb_stanlm_nosource <- rstanarm::stan_glm(total_biomass ~ (timeperiod) , data = some_sims, iter = 8000, thin = 4)
  tb_stanlm_notime <- rstanarm::stan_glm(total_biomass ~ 1 , data = some_sims, iter = 8000, thin = 4)


  tb_stanlms = list(
    tb_stanlm_full = tb_stanlm_full,
    tb_stanlm_nosource = tb_stanlm_nosource,
    tb_stanlm_notime = tb_stanlm_notime
  )

  return(list(
    te_stanlms = te_stanlms,
    tb_stanlms = tb_stanlms,
    matssname =some_sims$matssname[1],
    simtype = some_sims$simtype[1]
  ))

}

#' Compare fits
#'
#' @param stanlms_fits fits
#'
#' @return df
#' @export
#'
#' @importFrom rstanarm loo loo_compare
#' @importFrom dplyr mutate row_number
compare_stanarms <- function(stanlms_fits) {

  stanlms_loos<- lapply(stanlms_fits, rstanarm::loo, k_threshold= .7)

  stanlms_comparison <- rstanarm::loo_compare(stanlms_loos) %>%
    as.data.frame() %>%
    dplyr::mutate(model = row.names(.),
                  rank = dplyr::row_number())

  return(stanlms_comparison)

}

#' COmpare all fits
#'
#' @param some_stanlms_fits fits
#'
#' @return df
#' @export
#'
#' @importFrom dplyr bind_rows mutate
compare_both_stanarms<- function(some_stanlms_fits) {

  biomass <- compare_stanarms(stanlms_fits = some_stanlms_fits$tb_stanlms)
  energy <- compare_stanarms(stanlms_fits = some_stanlms_fits$te_stanlms)

  both_comparisons <- dplyr::bind_rows(biomass = biomass, energy = energy, .id = "currency") %>%
    dplyr::mutate(matssname = some_stanlms_fits$matssname,
                  simtype = some_stanlms_fits$simtype[1])


  return(both_comparisons)

}

