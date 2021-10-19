#' Calculate standardized effect size
#'
#' @param actual_value value
#' @param null_distribution distribution to calculate relative to
#'
#' @return ses
#' @export
#'
ses <- function(actual_value, null_distribution) {

  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)

  (actual_value - null_mean) / null_sd

}
#' Calculate percentile score
#'
#' @param actual_value value
#' @param null_distribution distribution to calculate relative to
#'
#' @return percentile
#' @export
#'
percentile_score <- function(actual_value, null_distribution) {

  sum(null_distribution < actual_value) / length(null_distribution)

}

#' Summarize actual v. null results
#'
#' @param null_results df
#' @param actual_results df
#'
#' @return df
#' @export
#'
#' @importFrom dplyr select rename left_join group_by summarize ungroup n
summarize_null_results <- function(null_results, actual_results) {

  null_use <- null_results %>%
    dplyr::select(statenum, route, routename, isd_turnover,mean_energy_raw_ratio,mean_biomass_raw_ratio) %>%
    dplyr::rename(null_isd_turnover = isd_turnover,
                  null_mean_energy_raw_ratio = mean_energy_raw_ratio,
                  null_mean_biomass_raw_ratio= mean_biomass_raw_ratio) %>%
    dplyr::left_join(
      dplyr::select(actual_results, statenum, route, routename, isd_turnover, mean_biomass_raw_ratio, mean_energy_raw_ratio))

  null_summarized <- null_use %>%
    dplyr::group_by(statenum, route, routename) %>%
    dplyr::summarize(
      ses = ses(unique(isd_turnover), null_isd_turnover),
      perc = percentile_score(unique(isd_turnover), null_isd_turnover),
      ses_e = ses(unique(mean_energy_raw_ratio), null_mean_energy_raw_ratio),
      ses_b = ses(unique(mean_biomass_raw_ratio), null_mean_biomass_raw_ratio),
      perc_e = percentile_score(unique(mean_energy_raw_ratio), null_mean_energy_raw_ratio),
      perc_b = percentile_score(unique(mean_biomass_raw_ratio), null_mean_biomass_raw_ratio),
      nsims = dplyr::n()
    ) %>%
    dplyr::ungroup()

  null_summarized
}
