#' Calculate state variables for every year given ISDs
#'
#' @param ts_isds result of BBSsize::simulate_isd_ts()
#'
#' @return data frame of state variables for each year, and seed used when simulating the ISD (if one was)
#' @export
#'
#' @importFrom dplyr group_by summarize n ungroup mutate
get_annual_svs <- function(ts_isds) {

  ts_isds %>%
    dplyr::mutate(ind_energy = estimate_b(mass)) %>%
    dplyr::group_by(year, isd_seed) %>%
    dplyr::summarize(abundance = dplyr::n(),
                     energy = sum(ind_energy),
                     biomass = sum(mass),
                     # median_energy = median(ind_energy),
                     #  median_biomass = median(mass),
                     mean_energy = mean(ind_energy),
                     mean_biomass = mean(mass)) %>%
    dplyr::ungroup()
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
#' @importFrom tidyr pivot_wider
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
  )%>%
    tidyr::pivot_wider(names_from = response, values_from = c(2:5))

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
    begin_years <- sort(unique(ts_svs$year))[1:5]
  }

  if(is.null(end_years)) {
    end_years <-  sort(unique(ts_svs$year))[ (length(unique(ts_svs$year)) - 4):length(unique(ts_svs$year))]
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
#' @importFrom tidyr pivot_wider
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
  ) %>%
    tidyr::pivot_wider(names_from = response, values_from = c(2:5))

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

  colnames(raw_results) <- paste0(colnames(raw_results), "_raw_ratio")


  return(raw_results)
}


#' Compute ISD turnover
#'
#' @param ts_isds result of simulate_isd_ts
#' @param begin_years optional
#' @param end_years optional
#'
#' @return df of turnover
#' @export
#' @importFrom dplyr filter mutate group_by summarize ungroup select bind_rows
compare_isds <- function(ts_isds, begin_years = NULL, end_years = NULL) {


  if(is.null(begin_years)) {
    begin_years <- sort(unique(ts_isds$year))[1:5]
  }

  if(is.null(end_years)) {
    end_years <- sort(unique(ts_isds$year))[ (length(unique(ts_isds$year)) - 4):length(unique(ts_isds$year))]
  }

  if(min(end_years) < max(begin_years)) {
    stop("End overlaps beginning")
  }

  begin_isd <- dplyr::filter(ts_isds, year %in% begin_years)
  end_isd <- dplyr::filter(ts_isds, year %in% end_years)
  #
  #   compare_isds <- dplyr::bind_rows(begin_isd, end_isd) %>%
  #     dplyr::mutate(timeperiod = ifelse(year %in% end_years, "end", "begin")) %>%
  #     dplyr::group_by(timeperiod) %>%
  #     dplyr::summarize(mean_size = mean(mass),
  #                      median_size = median(mass)) %>%
  #     dplyr::ungroup() %>%
  #     tidyr::pivot_wider(names_from = timeperiod, values_from = c(mean_size, median_size)) %>%
  #     dplyr::mutate(
  #       mean_size_ratio = mean_size_end / mean_size_begin,
  #       median_size_ratio = median_size_end / median_size_begin
  #     )

  begin_gmm <- add_gmm(begin_isd) %>%
    dplyr::mutate(timeperiod = "begin")
  end_gmm <- add_gmm(end_isd) %>%
    dplyr::mutate(timeperiod = "end")


  compare_gmms <- dplyr::bind_rows(begin_gmm, end_gmm)


  isd_overlap <- compare_gmms %>%
    dplyr::group_by(mass) %>%
    dplyr::summarize(mindensity = min(density)) %>%
    dplyr::ungroup() %>%
    dplyr::select(mindensity) %>%
    dplyr::summarize(isd_turnover =1- sum(mindensity))

  isd_overlap

}


#' Compare species composition
#'
#' @param ts_comp MATSS dataset
#' @param begin_years optional
#' @param end_years optional
#'
#' @return dataframe
#' @export
#'
#' @importFrom dplyr group_by summarize ungroup mutate select bind_rows
#' @importFrom vegan vegdist
compare_species_composition <- function(ts_comp, begin_years = NULL, end_years = NULL) {

  if(is.null(begin_years)) {
    begin_years <- ts_comp$covariates$year[1:5]
  }

  if(is.null(end_years)) {
    end_years <- ts_comp$covariates$year[ (length(ts_comp$covariates$year) - 4):length(ts_comp$covariates$year)]
  }

  if(min(end_years) < max(begin_years)) {
    stop("End overlaps beginning")
  }

  begin_rows <- which(ts_comp$covariates$year %in% begin_years)


  end_rows <- which(ts_comp$covariates$year %in% end_years)
  begin_composition <- colSums(ts_comp$abundance[begin_rows, ])

  end_composition <- colSums(ts_comp$abundance[end_rows, ])

  begin_relabund <- begin_composition / sum(begin_composition)

  end_relabund <- end_composition / sum(end_composition)

  relabund <- data.frame(
    begin = begin_relabund,
    end = end_relabund,
    beginsp = names(begin_relabund),
    endsp = names(end_relabund)
  )

  relabund_change <- relabund %>%
    dplyr::group_by(beginsp) %>%
    dplyr::summarize(minRel = min(begin, end)) %>%
    dplyr::ungroup() %>%
    dplyr::select(minRel) %>%
    dplyr::summarize(sp_turnover = 1-sum(minRel))


  be_matrix <- dplyr::bind_rows(begin_composition, end_composition)

  be_diss <- vegan::vegdist(be_matrix)

  relabund_change <- relabund_change %>%
    dplyr::mutate(bcd = be_diss)

  relabund_change


}

#' Interaction of sv change
#'
#' @param caps_svs caps_svs
#'
#' @return df
#' @export
#'
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_wider
interaction_lms <- function(caps_svs) {

  rangescale <- function(vect) {

    vectrange <- max(vect) - min(vect)

    vect = (vect - min(vect)) / vectrange

    vect

  }

  caps_svs <- caps_svs %>%
    dplyr::mutate(energy = rangescale((energy)),
                  abundance = rangescale((abundance)),
                  biomass = rangescale(biomass))

  caps_long <- caps_svs %>%
    tidyr::pivot_longer(c(-year, -timeperiod), names_to = "currency", values_to = "val")%>%
    dplyr::filter(currency %in% c("energy", "abundance", "biomass"))

  cap_lm <- lm(val ~ timeperiod * currency, data = caps_long)

  cap_lm_results <- summary(cap_lm)

  cap_lm_ps <- cap_lm_results$coefficients


  cap_p <- pf(cap_lm_results$fstatistic[1], cap_lm_results$fstatistic[2], cap_lm_results$fstatistic[3], lower.tail = F)
  cap_lm_results_wide <- cap_lm_ps %>%
    as.data.frame() %>%
    dplyr::mutate(coef_name = row.names(.)) %>%
    tidyr::pivot_wider(names_from = coef_name, values_from = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) %>%
    dplyr::mutate(overall_p = cap_p,
                  overall_r2 = cap_lm_results$r.squared)

  cap_lm_results_wide


}

#' Run and collect all core analyses
#'
#' @param ts_comp a MATSS_style dataset
#' @param begin_years optional
#' @param end_years optional
#' @param isd_seed optional
#'
#' @return results
#' @export
#'
#' @importFrom dplyr bind_cols mutate
#' @importFrom BBSsize simulate_isd_ts
all_core_analyses <- function(ts_comp, begin_years = NULL, end_years = NULL, isd_seed = NULL) {

  ts_isd <- BBSsize::simulate_isd_ts(ts_comp, isd_seed = isd_seed)
  ts_svs <- get_annual_svs(ts_isd$isd)
  ts_lms <- fit_all_timeseries_lms(ts_svs)
  caps_svs <- pull_caps(ts_svs, begin_years, end_years)
  caps_lms <- fit_all_caps_lms(caps_svs)
  i_lms <- interaction_lms(caps_svs)
  raw_ratios <-  compute_raw_sv_change(caps_svs)
  set.seed(1977)
  isd_turn <- compare_isds(ts_isd$isd, begin_years, end_years)
  comp_turn <- compare_species_composition(ts_comp, begin_years, end_years)


  all_results <- dplyr::bind_cols(ts_lms, caps_lms, raw_ratios, i_lms, isd_turn, comp_turn, as.data.frame(ts_comp$metadata$location)) %>%
    dplyr::mutate(beginyears = toString(begin_years),
                  endyears = toString(end_years))
  return(all_results)
}
