#' Slice a dataset to given start and end
#'
#' @param dataset from MATSS
#' @param startyear startyear
#' @param endyear endyear
#'
#' @return dataset subset to start:end years
#' @export
#'
slice_dataset <- function(dataset, startyear, endyear) {

  keep_years = startyear:endyear

  keep_rows <- which(dataset$covariates$year %in% keep_years)

  dataset$abundance <- dataset$abundance[keep_rows,]
  dataset$covariates <- dataset$covariates[keep_rows,]

  return(dataset)

}
