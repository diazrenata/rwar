#' SVs from ISD
#'
#' @param isd list
#'
#' @return df of svs
#' @export
#' @importFrom dplyr mutate summarize n
summarize_isd <- function(isd) {

  isd_sv <- isd$isd %>%
    dplyr::mutate(energy = BBSsize::estimate_b(mass)) %>%
    dplyr::summarize(
      richness = length(unique(id)),
      abundance = dplyr::n(),
      energy = sum(energy),
      biomass = sum(mass),
      mean_energy = sum(energy) / dplyr::n(),
      mean_biomass = sum(mass) / dplyr::n(),
      years = toString(unique(year)))

  isd_sv_with_id <- isd_sv %>%
    dplyr::mutate(route = isd$metadata$route[1],
                  region = isd$metadata$region[1],
                  sim_seed = ifelse("seed" %in% names(isd$metadata), isd$metadata$seed, NA))

  return(isd_sv_with_id)

}
