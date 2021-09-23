#' Shuffle species for regional null model
#'
#' @param ts_dat matss-style
#' @param ranges_dat ranges
#' @param null_mod_seed or null
#'
#' @return matss-style with additional metadata
#' @export
#'
#' @importFrom dplyr left_join
shuffle_regional <- function(ts_dat, ranges_dat, null_mod_seed) {

  ts_locationdat <- ts_dat$metadata$location

  ts_speciesoverlap <- dplyr::left_join(ts_locationdat, ranges_dat)

  overlapping_spp <- ts_speciesoverlap$id
  encountered_spp <- ts_dat$metadata$species_table$id
  all_ts_spp = c(overlapping_spp, encountered_spp) %>% unique()

  overlap_rich = length(overlapping_spp)
  encountered_rich = length(encountered_spp)
  overlap_not_in_encountered = sum(!(overlapping_spp %in% encountered_spp))

  set.seed(null_mod_seed)

  shuffled_species <- sample(all_ts_spp, size = (encountered_rich), replace = F)

  shuffled_dat <- ts_dat

  colnames(shuffled_dat$abundance) <- shuffled_species

  shuffled_dat$metadata$regional_pool <- list(overlap_rich = overlap_rich,
                                              encountered_rich = encountered_rich,
                                              overlap_not_in_encountered = overlap_not_in_encountered)

  return(shuffled_dat)
}

#' Regional null model
#'
#' Reassign species IDs using all species whose ranges overlap a route, or were ever observed on that route.
#'
#' @param ts_dat a MATSS-shaped dataset
#' @param ranges_dat a dataframe of route & species encountered each route
#' @param null_mod_seed seed to use when reshuffling species. if none given will be drawn.
#' @param sim_index sim index
#' @param begin_years to pass to all_core_analyses
#' @param end_years to pass to all_core_analyses
#' @param isd_seed to pass to all_core_analyses
#'
#' @return results of all_core_analyses
#' @export
#'
#' @importFrom dplyr mutate
regional_null_model <- function(ts_dat, ranges_dat = NULL, null_mod_seed = NULL, sim_index = NULL, begin_years = NULL, end_years = NULL, isd_seed = NULL) {

  if(is.null(ranges_dat)) {
    stop("Need ranges data")
  }

  if(is.null(null_mod_seed)) {
    set.seed(NULL)
    null_mod_seed <- sample.int(1000000000, 1)
  }


  shuffled_dat <- shuffle_regional(ts_dat = ts_dat, ranges_dat = ranges_dat, null_mod_seed = null_mod_seed)

  results <- all_core_analyses(shuffled_dat, begin_years, end_years, isd_seed)

  results <- results %>%
    dplyr::mutate(
      null_mod_type = "regional",
      null_mod_seed = null_mod_seed,
      sim_index = sim_index,
      overlap_richness = shuffled_dat$metadata$regional_pool$overlap_rich,
      local_richness = shuffled_dat$metadata$regional_pool$encountered_rich,
      regionally_added = shuffled_dat$metadata$regional_pool$overlap_not_in_encountered
    )


  results
}

#' Wrapper for regional null model
#'
#' @param ts_dat matss dataset
#' @param ranges_dat ranges dat
#' @param nsims default 100
#' @param begin_years pass
#' @param end_years pass
#'
#' @return df
#' @export
#'
#' @importFrom dplyr bind_rows
regional_null_model_wrapper <- function(ts_dat, ranges_dat = NULL, nsims = 100, begin_years = NULL, end_years = NULL) {

  repeated_nulls <- lapply(1:nsims, regional_null_model, ts_dat = ts_dat, ranges_dat = ranges_dat, null_mod_seed = NULL, begin_years = begin_years, end_years = end_years, isd_seed = NULL)

  results <- dplyr::bind_rows(repeated_nulls)

  results

}
