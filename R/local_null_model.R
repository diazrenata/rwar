#' Local null model
#'
#' Reassign species IDs using all species ever observed on a single route.
#'
#' @param ts_dat a MATSS-shaped dataset
#' @param null_mod_seed seed to use when reshuffling species. if none given will be drawn.
#' @param begin_years to pass to all_core_analyses
#' @param end_years to pass to all_core_analyses
#' @param isd_seed to pass to all_core_analyses
#'
#' @return results of all_core_analyses
#' @export
#'
#' @importFrom dplyr mutate
local_null_model <- function(ts_dat, null_mod_seed = NULL, sim_index = NULL, begin_years = NULL, end_years = NULL, isd_seed = NULL) {

  if(is.null(null_mod_seed)) {
    set.seed(NULL)
    null_mod_seed <- sample.int(1000000000, 1)
  }


  orig_species <- colnames(ts_dat$abundance)

  set.seed(null_mod_seed)

  shuffled_species <- sample(orig_species, size = length(orig_species), replace = F)

  shuffled_dat <- ts_dat

  colnames(shuffled_dat$abundance) <- shuffled_species

  results <- all_core_analyses(shuffled_dat, begin_years, end_years, isd_seed)

  results <- results %>%
    dplyr::mutate(
      null_mod_type = "local",
      null_mod_seed = null_mod_seed,
      sim_index = sim_index
    )


  results
}
