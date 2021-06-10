#' Reassign species randomly
#'
#' Will preserve turnover, SAD, but mess with ISD agnostic to size.
#'
#' @param dataset a dataset
#' @param seed optional
#'
#' @return modified dataset
#' @export
shuffle_species <- function(dataset, seed = NULL) {

  if(is.null(seed)) {
    seed <- sample.int(1000000000, size = 1)
  }

  abundance_reshuffled <- dataset$abundance

  newspecies <- sample(colnames(abundance_reshuffled), size = ncol(abundance_reshuffled), replace = F)

  colnames(abundance_reshuffled) <- newspecies

  shuffled_dat <- dataset
  shuffled_dat$abundance <- abundance_reshuffled

  shuffled_dat$metadata$seed = seed

  return(shuffled_dat)

}
