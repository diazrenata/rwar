#' Get portal rats
#'
#' @return list
#' @export
#'
#' @importFrom cvlt get_rodents_annual
get_portal_rats <- function() {

  individual_rats <- cvlt::get_rodents_annual()

  ratyears <- individual_rats$covariates$year

  keepyears <- c(1980:2010)

  keeprows <- which(ratyears %in% keepyears)

  individual_rats$abundance <- individual_rats$abundance[keeprows, ]
  individual_rats$covariates <- individual_rats$covariates[keeprows, ]

  individual_rats$metadata$route <- 1977
  individual_rats$metadata$region <- 1977
  individual_rats$metadata$location.routename <- "PORTAL RATS"
  individual_rats$metadata$location.bcr <- 1977
  individual_rats$metadata$location.statenum <- 1977


  return(individual_rats)
}

#' Get rat sds
#'
#' @return df
#' @export
#'
#' @importFrom portalr summarize_individual_rodents
#' @importFrom dplyr filter group_by summarize ungroup rename
get_portal_sd_dat <- function() {

  ind_rats <- portalr::summarise_individual_rodents()

  ind_rats <- ind_rats %>%
    dplyr::filter(!is.na(wgt), !is.na(species)) %>%
    dplyr::group_by(species) %>%
    dplyr::summarize(mean_mass = mean(wgt),
                     mean_sd = sd(wgt)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(id = species)

  return(ind_rats)

}
