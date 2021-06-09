#' Fit a GMM
#'
#' @param size_vect vector of sizes (logged or not)
#'
#' @return GMM fit with up to 15 Gaussians
#' @export
#' @importFrom mclust densityMclust
fit_gmm <- function(size_vect) {
  library(mclust)
  this_gmm <- mclust::densityMclust(size_vect, G = c(1:15) )

  return(this_gmm)

}

#' Add GMM density to an ISD
#'
#' Fits to log of mass
#'
#' @param isd from BBSsize::simulate_isd_ts
#'
#' @return dataframe with columns mass, density
#' @export
#'
#' @importFrom dplyr mutate
add_gmm <- function(isd) {

  isd <- isd %>%
    dplyr::mutate(logmass = log(mass))

  gmm <- fit_gmm(isd$logmass)


  gmm_isd <- data.frame(logmass = seq(0, log(15000), length.out = 1000))
  gmm_isd$dens <- predict(gmm, newdata = gmm_isd$logmass)


  isd_gmm <- data.frame(
    mass = gmm_isd$logmass,
    density = (gmm_isd$dens)/ sum(gmm_isd$dens)
  )

  return(isd_gmm)

}
