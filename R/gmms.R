#' Fit a GMM
#'
#' @param size_vect vector of sizes (logged or not)
#' @param max_G number of Gaussians
#'
#' @return GMM fit with up to max_G Gaussians
#' @export
#' @importFrom mclust densityMclust
fit_gmm <- function(size_vect, max_G = 15) {
  library(mclust)
  this_gmm <- mclust::densityMclust(size_vect, G = c(1:max_G) )

  return(this_gmm)

}

#' Add GMM density to an ISD
#'
#' Fits to log of mass
#'
#' @param isd from BBSsize::simulate_isd_ts
#' @param max_size default 15000 but changeable for mammals
#' @param max_G max n gaussians
#'
#' @return dataframe with columns mass, density
#' @export
#'
#' @importFrom dplyr mutate
add_gmm <- function(isd, max_size = 15000, max_G = 15) {

  isd <- isd %>%
    dplyr::mutate(logmass = log(mass))

  gmm <- fit_gmm(isd$logmass, max_G)


  gmm_isd <- data.frame(logmass = seq(0, log(max_size), length.out = 1000))
  gmm_isd$dens <- predict(gmm, newdata = gmm_isd$logmass)


  isd_gmm <- data.frame(
    mass = gmm_isd$logmass,
    density = (gmm_isd$dens)/ sum(gmm_isd$dens)
  )

  return(isd_gmm)

}
