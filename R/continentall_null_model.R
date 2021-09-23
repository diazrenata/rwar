#' Shuffle species for continental null model
#'
#' @param ts_dat matss sytle
#' @param null_mod_seed seed
#'
#' @return shuffled
#' @export
shuffle_continental <- function(ts_dat, null_mod_seed) {

  orig_species <- colnames(ts_dat$abundance)

  all_species <- BBSsize::sd_table$id

  set.seed(null_mod_seed)

  shuffled_species <- sample(all_species, size = length(orig_species), replace = F)

  shuffled_dat <- ts_dat

  colnames(shuffled_dat$abundance) <- shuffled_species

  return(shuffled_dat)
}

#' Continental null model
#'
#' Reassign species IDs using all species ever observed on a single route.
#'
#' @param ts_dat a MATSS-shaped dataset
#' @param null_mod_seed seed to use when reshuffling species. if none given will be drawn.
#' @param sim_index sim index to track
#' @param begin_years to pass to all_core_analyses
#' @param end_years to pass to all_core_analyses
#' @param isd_seed to pass to all_core_analyses
#'
#' @return results of all_core_analyses
#' @export
#'
#' @importFrom dplyr mutate
continental_null_model <- function(ts_dat, null_mod_seed = NULL, sim_index = NULL, begin_years = NULL, end_years = NULL, isd_seed = NULL) {

  if(is.null(null_mod_seed)) {
    set.seed(NULL)
    null_mod_seed <- sample.int(1000000000, 1)
  }

  shuffled_dat <- shuffle_continental(ts_dat = ts_dat, null_mod_seed = null_mod_seed)

  results <- all_core_analyses(shuffled_dat, begin_years, end_years, isd_seed)

  results <- results %>%
    dplyr::mutate(
      null_mod_type = "continental",
      null_mod_seed = null_mod_seed,
      sim_index = sim_index
    )


  results
}

#' Wrapper for continental null model
#'
#' @param ts_dat matss dataset
#' @param nsims default 100
#' @param begin_years pass
#' @param end_years pass
#'
#' @return df
#' @export
#'
#' @importFrom dplyr bind_rows
continental_null_model_wrapper <- function(ts_dat, nsims = 100, begin_years = NULL, end_years = NULL) {

  repeated_nulls <- lapply(1:nsims, continental_null_model, ts_dat = ts_dat, null_mod_seed = NULL, begin_years = begin_years, end_years = end_years, isd_seed = NULL)

  results <- dplyr::bind_rows(repeated_nulls)

  results

}
