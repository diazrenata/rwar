
#' Get GMMs or KDEs fit to first and last 5 years of data
#'
#' @param dataset from BBSsize/MATSS
#' @param smooth_method "gmm" or "kde", currently only implemented for "gmm"
#'
#' @return dataframe of mass, density start, density end, and density change, plus route identifying info
#' @export
#'
#' @importFrom dplyr mutate select distinct bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom BBSsize simulate_isd_ts

get_begin_end_smooths <- function(begin_end_isds, smooth_method = "gmm") {

  isd1 <- begin_end_isds$isd1
  isd2 <- begin_end_isds$isd2

  if(smooth_method == "kde") {
    isd1_kde <- add_kde(isd1$isd) %>%
      dplyr::mutate(chunk = "start")
    isd2_kde <- add_kde(isd2$isd) %>%
      dplyr::mutate(chunk = "end")

    both_isds <- dplyr::bind_rows(isd1_kde, isd2_kde)

  } else if(smooth_method == "gmm") {
    isd1_gmm <- add_gmm(isd1$isd) %>%
      dplyr::mutate(chunk = "start")
    isd2_gmm <- add_gmm(isd2$isd) %>%
      dplyr::mutate(chunk = "end")

    both_isds <- dplyr::bind_rows(isd1_gmm, isd2_gmm)
  }

  both_isds <- both_isds %>%
    tidyr::pivot_wider(id_cols = mass, names_from = chunk, values_from = density) %>%
    dplyr::mutate(density_diff = end - start,
                  density_rat = end / start)

  both_isds_with_id <- as.data.frame(begin_end_isds$metadata) %>%
    dplyr::select(route, region, location.bcr, location.statenum, location.routename) %>%
    dplyr::distinct() %>%
    cbind(both_isds) %>%
    dplyr::mutate(startyears = begin_end_isds$metadata$startyears,
                  endyears =  begin_end_isds$metadata$endyears)
  return(both_isds_with_id)

}

#' Get sv summaries for start/end
#'
#' @param begin_end_isds list
#'
#' @return df
#' @export
#' @importFrom dplyr bind_rows
get_begin_end_svs <- function(begin_end_isds) {

  start = summarize_isd(begin_end_isds$isd1)
  end = summarize_isd(begin_end_isds$isd2)

  long= (dplyr::bind_rows(start = start, end = end, .id = "timechunk"))

  return(long)
}

#' Get ISDs for first and last 5 years
#'
#' @param dataset dataset
#'
#' @return list
#' @export
#'
#' @importFrom BBSsize simulate_isd_ts
#'
get_begin_end_isds <- function(dataset) {

  library(BBSsize)

  startyears <- dataset$covariates$year[1:5]

  endyears <- dataset$covariates$year[(nrow(dataset$covariates)-4) : nrow(dataset$covariates)]


  isd1 <- BBSsize::simulate_isd_ts(dataset, censusyears = startyears)
  isd2 <- BBSsize::simulate_isd_ts(dataset, censusyears = endyears)

  new_metadata <- dataset$metadata
  new_metadata$startyears = toString(startyears)
  new_metadata$endyears = toString(endyears)

  return(list(isd1 = isd1, isd2 = isd2, metadata = new_metadata))

}

#' Get smooths and svs for a pair of isds
#'
#' @param begin_end_isds list
#'
#' @return list of dfs for sv, smooths, metadata
#' @export
#'
process_begin_end_isds <- function(begin_end_isds) {

  gmm_smooths <- get_begin_end_smooths(begin_end_isds)

  svs <- get_begin_end_svs(begin_end_isds)

  return(list(
    smooths = gmm_smooths,
    svs = svs,
    metadata = begin_end_isds$metadata
  ))

}

#' Estimate b (energy)
#'
#' @param body_mass mass
#' @return estimated b (metabolic rate) using pars from Fristoe 2015
#' @export

estimate_b <- function(body_mass) {
  return(10.5 * (body_mass ^ .713))
}

#' SVs from ISD
#'
#' @param isd list
#'
#' @return df of svs
#' @export
#' @importFrom dplyr mutate summarize n
summarize_isd <- function(isd) {

  isd_sv <- isd$isd %>%
    dplyr::mutate(energy = estimate_b(mass)) %>%
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
                  region = isd$metadata$region[1])

  return(isd_sv_with_id)

}

#' Calculate overlap
#'
#' After NEON paper
#'
#' @param begin_end_smooths df
#'
#' @return df
#' @export
#'
#' @importFrom dplyr select group_by mutate ungroup distinct
overlap = function(begin_end_smooths) {

  smooths <- begin_end_smooths %>%
    dplyr::select(mass, start, end) %>%
    dplyr::group_by(mass) %>%
    dplyr::mutate(mindensity = min(start, end)) %>%
    dplyr::ungroup()

  overlap = sum(smooths$mindensity)

  overlap_df <- begin_end_smooths %>%
    dplyr::select(route, region, location.bcr, location.statenum, location.routename, startyears, endyears) %>%
    dplyr::mutate(overlap = overlap) %>%
    dplyr::distinct()

  return(overlap_df)
}


#' Get compositional turnover
#'
#' @param dataset dataset
#'
#' @return df
#' @export
#'
#' @importFrom dplyr mutate group_by summarize ungroup distinct
#'
get_begin_end_composition <- function(dataset) {

  startyears <- dataset$covariates$year[1:5]

  endyears <- dataset$covariates$year[(nrow(dataset$covariates)-4) : nrow(dataset$covariates)]


  start <- dataset$abundance[ which(dataset$covariates$year %in% startyears), ]
  end <- dataset$abundance[ which(dataset$covariates$year %in% endyears), ]

  start <- data.frame(
    species = colnames(start),
    abundance = colSums(start),
    timechunk = "start"
  )

  start <- start %>%
    dplyr::mutate(total = sum(abundance)) %>%
    dplyr::mutate(relative = abundance / total)


  end <- data.frame(
    species = colnames(end),
    abundance = colSums(end),
    timechunk = "end"
  )

  end <- end %>%
    dplyr::mutate(total = sum(abundance)) %>%
    dplyr::mutate(relative = abundance / total)


  composition <- dplyr::bind_rows(start,end)

  composition_turnover <- composition %>%
    dplyr::group_by(species) %>%
    dplyr::summarize(min_relative = min(relative)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(composition_overlap = sum(min_relative))%>%
    dplyr::mutate(route = dataset$metadata$route[1],
                  region = dataset$metadata$region[1],
                  location.bcr = dataset$metadata$location$bcr)


  return(composition_turnover)

}
