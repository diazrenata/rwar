#' Calculate state variables for every year given ISDs
#'
#' @param ts_isds result of BBSsize::simulate_isd_ts()
#'
#' @return data frame of state variables for each year, route region info, and seed used when simulating the ISD (if one was)
#' @export
#'
#' @importFrom dplyr filter mutate
get_annual_svs <- function(ts_isds) {

  annual_isds <- lapply(ts_isds$covariates$year, FUN = function(a_year, all_isds) return(list(isd = dplyr::filter(all_isds, year == a_year), covariates = ts_isds$covariates, metadata = ts_isds$metadata)), all_isds = ts_isds$isd)

  isd_summaries <- lapply(annual_isds, rwar::summarize_isd)

  dplyr::bind_rows(isd_summaries) %>%
    dplyr::mutate(year = as.numeric(years))

}

#' Fit a lm to the full timeseries
#'
#' @param ts_svs result of get_annual_svs, must include time as `year`
#' @param response_name the name of the variable to compare to year, must be the name of a column in ts_svs
#'
#' @return data frame of response_name, p-value, r2, slope, and fitted ratio (last value / first fitted value) for `lm(response_name ~ year)`
#' @export
#'
fit_timeseries_lm <- function(ts_svs, response_name) {

  if(!(response_name %in% colnames(ts_svs))) {
    stop("Response variable missing")
  }

  if(!("year" %in% colnames(ts_svs))) {
    stop("year missing")
  }

  ts_d <- ts_svs[, c("year", response_name)]

  colnames(ts_d)[2] <- "response"

  ts_lm <- lm(response ~ year, data = ts_d)

  lm_p <- anova(ts_lm)[1,5]

  lm_r2 <- summary(ts_lm)$r.squared

  lm_slope <- coef(ts_lm)[2]

  lm_fitted <- ts_lm$fitted.values

  first_year <- 1
  last_year <- length(lm_fitted)

  lm_fitted_ratio <- lm_fitted[last_year] / lm_fitted[first_year]

  return(data.frame(
    response = response_name,
    p_ts = lm_p,
    r2_ts = lm_r2,
    slope_ts = lm_slope,
    fitted_ratio_ts = lm_fitted_ratio
  ))

}

#' Fit linear models to the 5 state variables for the full timeseries
#'
#' @param ts_svs dataframe, must contain columns `year`, `abundance`, `biomass`, `energy`, `mean_biomass`, `mean_energy`
#'
#' @return dataframe of lm results
#' @export
#'
#' @importFrom dplyr bind_rows
fit_all_timeseries_lms <- function(ts_svs) {

  n_lm <- fit_timeseries_lm(ts_svs, "abundance")
  e_lm <- fit_timeseries_lm(ts_svs, "energy")
  b_lm <- fit_timeseries_lm(ts_svs, "biomass")
  me_lm <- fit_timeseries_lm(ts_svs, "mean_energy")
  mb_lm <- fit_timeseries_lm(ts_svs, "mean_biomass")


  all_lm <- dplyr::bind_rows(
    n_lm,
    e_lm,
    b_lm,
    me_lm,
    mb_lm
  )

  return(all_lm)
}

#' Pull first and last 5 years
#'
#' @param ts_svs dataframe, must have col `year`
#' @param begin_years vector of begin years, or will use first 5 years of ts
#' @param end_years vector of end years, or will use last 5 years of ts
#'
#' @return ts filtered to years in begin or end and with column `timeperiod` for "begin", "end"
#' @export
#'
#' @importFrom dplyr filter mutate
pull_caps <- function(ts_svs, begin_years = NULL, end_years = NULL){

  if(is.null(begin_years)) {
    begin_years <- ts_svs$year[ (min(ts_svs$year)):(max(ts_svs$year) + 4)]
  }

  if(is.null(end_years)) {
    end_years <- ts_svs$year[ (max(ts_svs$year) - 4):(max(ts_svs$year))]
  }

  if(min(end_years) < max(begin_years)) {
    stop("End overlaps beginning")
  }

  dplyr::filter(ts_svs, year %in% c(begin_years, end_years)) %>%
    dplyr::mutate(timeperiod = ifelse(year > max(begin_years), "end", "begin"))

}


#' Fit a lm to the first and last 5 years ("caps")
#'
#' @param caps_svs result of pull_caps(get_annual_svs), must include timeperiod as `timeperiod`
#' @param response_name the name of the variable to compare to year, must be the name of a column in caps_svs
#'
#' @return data frame of response_name, p-value, r2, slope, and fitted ratio (last value / first fitted value) for `lm(response_name ~ timeperiod)`
#' @export
#'
fit_caps_lm <- function(caps_svs, response_name) {

  if(!(response_name %in% colnames(caps_svs))) {
    stop("Response variable missing")
  }

  if(!("timeperiod" %in% colnames(caps_svs))) {
    stop("timeperiod missing")
  }

  caps_d <- caps_svs[, c("timeperiod", response_name)]

  colnames(caps_d)[2] <- "response"

  caps_lm <- lm(response ~ timeperiod, data = caps_d)

  lm_p <- anova(caps_lm)[1,5]

  lm_r2 <- summary(caps_lm)$r.squared

  lm_slope <- coef(caps_lm)[2]

  lm_fitted <- caps_lm$fitted.values

  first_year <- 1
  last_year <- length(lm_fitted)

  lm_fitted_ratio <- lm_fitted[last_year] / lm_fitted[first_year]

  return(data.frame(
    response = response_name,
    p_caps = lm_p,
    r2_caps = lm_r2,
    slope_caps = lm_slope,
    fitted_ratio_caps = lm_fitted_ratio
  ))

}



#' Fit linear models to the 5 state variables for the first and last 5 years ("caps")
#'
#' @param caps_svs result of pull_caps(get_annual_svs), must include timeperiod as `timeperiod` and cols `abundance`, `biomass`, `energy`, `mean_biomass`, `mean_energy`
#'
#' @return dataframe of lm results
#' @export
#'
#' @importFrom dplyr bind_rows
fit_all_caps_lms <- function(caps_svs) {

  n_lm <- fit_caps_lm(caps_svs, "abundance")
  e_lm <- fit_caps_lm(caps_svs, "energy")
  b_lm <- fit_caps_lm(caps_svs, "biomass")
  me_lm <- fit_caps_lm(caps_svs, "mean_energy")
  mb_lm <- fit_caps_lm(caps_svs, "mean_biomass")


  all_lm <- dplyr::bind_rows(
    n_lm,
    e_lm,
    b_lm,
    me_lm,
    mb_lm
  )

  return(all_lm)
}

#' Compute raw change from the first five to the last five years
#'
#' @param caps_svs result of pull_caps(get_annual_svs), must include timeperiod as `timeperiod` and cols `abundance`, `biomass`, `energy`, `mean_biomass`, `mean_energy`
#'
#' @return dataframe of `response` and `mean value`
#' @export
#'
#' @importFrom dplyr group_by summarize mutate ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
compute_raw_sv_change <- function(caps_svs) {


  raw_change <- caps_svs %>%
    dplyr::group_by(timeperiod) %>%
    dplyr::summarize(energy = sum(energy),
              abundance = sum(abundance),
              biomass = sum(biomass)) %>%
    dplyr::mutate(mean_energy = energy / abundance,
           mean_biomass = biomass / abundance) %>%
    dplyr::ungroup()

  raw_results <- raw_change[2, 2:6] / raw_change[1, 2:6]

  raw_results <- raw_results %>%
    tidyr::pivot_longer(tidyselect::everything(), names_to = "response", values_to = "raw_ratio")


  return(raw_results)
}
