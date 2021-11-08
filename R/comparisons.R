
#' Compare ISD, sp composition
#'
#' @param dat matss
#' @param initial_isd_seed 1989
#' @param shuffle_seed 1989
#' @param begin_years usual
#' @param end_years usual
#' @param nshuffles n for ks test
#' @param n_isd_draws n for gmms
#'
#' @return dataframe
#' @export
#' @importFrom dplyr bind_cols
be_comparison <- function(dat, initial_isd_seed = 1989, shuffle_seed = 1989, begin_years = 1988:1992, end_years = 2014:2018, nshuffles = 500, n_isd_draws = 5) {

  ks_comp <- ks_comparison(dat, initial_isd_seed, begin_years, end_years, shuffle_seed, nshuffles)
  overlap_comp <- overlap_comparison(dat, initial_isd_seed, n_isd_draws, begin_years, end_years)

  sp_comp <- sp_comparison(dat, begin_years, end_years)


  out <- data.frame(
    matssname = paste0("bbs_rtrg_", dat$metadata$route, "_", dat$metadata$region)
  ) %>%
    dplyr::bind_cols(ks_comp, overlap_comp, sp_comp)

  out
}

#' Compare species composition begin-end
#'
#' @param dat matss
#' @param begin_years usual
#' @param end_years usual
#'
#' @return df
#' @export
#' @importFrom vegan vegdist
#' @importFrom dplyr mutate group_by ungroup group_by_all summarize
#' @importFrom tidyr pivot_longer pivot_wider
sp_comparison <- function(dat, begin_years = 1988:1992, end_years = 2014:2018) {

  begin_comp <- dat$abundance[ which(dat$covariates$year %in% begin_years), ]
  end_comp <- dat$abundance[ which(dat$covariates$year %in% end_years), ]

  begin_totals <- colSums(begin_comp)
  end_totals <- colSums(end_comp)

  commatrix <- dplyr::bind_rows(begin_totals, end_totals) %>% as.matrix()

  bcd <- vegan::vegdist(commatrix)

  sp_overlap <- dplyr::bind_rows(begin_totals, end_totals) %>%
    dplyr::mutate(timeperiod = c("begin", "end")) %>%
    tidyr::pivot_longer(-timeperiod, names_to = "sp", values_to= "abund") %>%
    dplyr::group_by(timeperiod) %>%
    dplyr::mutate(total_abund = sum(abund)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(relAbund = abund / total_abund) %>%
    tidyr::pivot_wider(id_cols = sp,
                       names_from = timeperiod,
                       values_from = relAbund) %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(minRelAbund = min(begin, end)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(sp_overlap = sum(minRelAbund))

  sp_overlap

}

#' Compare isd overlap
#'
#' @param dat mats
#' @param initial_isd_seed passed to sampling gmm
#' @param n_isd_draws passed to sampling gmm
#' @param begin_years ditto
#' @param end_years ditto
#'
#' @return df
#' @export
#'
#' @importFrom dplyr left_join group_by_all mutate ungroup summarize
overlap_comparison <- function(dat, initial_isd_seed = 1989, n_isd_draws = 5, begin_years = 1988:1992, end_years = 2014:2018) {

  sampling_gmms <- construct_sampling_gmm(dat, n_isd_draws = n_isd_draws, initial_isd_seed = initial_isd_seed, begin_years = begin_years, end_years = end_years)

  isd_overlap <- sampling_gmms$begin %>%
    dplyr::left_join(sampling_gmms$end, by = "mass") %>%
    dplyr::group_by_all() %>%
    dplyr::mutate(mindensity = min(density.x, density.y)) %>%
    dplyr::ungroup() %>%
    dplyr::summarize(overlap = sum(mindensity))

  isd_overlap
}

#' KS comparison bootstrapped
#'
#' @param dat matss
#' @param initial_isd_seed to draw isd
#' @param begin_years usual
#' @param end_years usual
#' @param shuffle_seed for shuffles
#' @param nshuffles how many bootstraps
#'
#' @return df
#' @export
#'
#' @importFrom dplyr bind_rows filter
ks_comparison <- function(dat, initial_isd_seed = 1989, begin_years = 1988:1992, end_years = 2014:2018, shuffle_seed= 1989, nshuffles = 500) {

  raw_isds <- just_isd(dat, initial_isd_seed)

  begin_raw_isd <- raw_isds %>%
    dplyr::filter(year %in% begin_years)
  end_raw_isd <- raw_isds %>%
    dplyr::filter(year %in% end_years)

  begin_mean_size <- mean(sqrt(begin_raw_isd$mass))
  end_mean_size <- mean(sqrt(end_raw_isd$mass))

  begin_mean_e <- mean(sqrt(estimate_b(begin_raw_isd$mass)))
  end_mean_e <- mean(sqrt(estimate_b(end_raw_isd$mass)))

  begin_median_size <- median(begin_raw_isd$mass)
  end_median_size <- median(end_raw_isd$mass)

  begin_median_e <- median(estimate_b(begin_raw_isd$mass))
  end_median_e <- median(estimate_b(end_raw_isd$mass))

  actual_ks_test = ks.test(begin_raw_isd$mass, end_raw_isd$mass)

  begin_nind = nrow(begin_raw_isd)
  end_nind = nrow(end_raw_isd)

  all_ind = c(begin_raw_isd$mass, end_raw_isd$mass)

  total_nind = sum(begin_nind, end_nind)

  shuffles <- list()

  set.seed(shuffle_seed)

  for(i in 1:nshuffles){

    shuffle_begin <- sample.int(total_nind, size = begin_nind, replace = F)

    begin_masses <- all_ind[shuffle_begin]
    end_masses <- all_ind[-shuffle_begin]

    shuffle_ks <- ks.test(begin_masses, end_masses)

    shuffles[[i]] <- data.frame(
      it = i,
      d = shuffle_ks$statistic,
      p = shuffle_ks$p.value
    )
  }

  shuffles <- dplyr::bind_rows(shuffles)

  ks_results <- data.frame(
    actual_d = actual_ks_test$statistic,
    actual_p = actual_ks_test$p.value,
    d_ses = ses(actual_ks_test$statistic, shuffles$d),
    d_percentile = percentile_score(actual_ks_test$statistic, shuffles$d),
    begin_mean_size = begin_mean_size,
    end_mean_size  = end_mean_size,
    begin_mean_e = begin_mean_e,
    end_mean_e = end_mean_e,
    begin_median_size = begin_median_size,
    end_median_size  = end_median_size,
    begin_median_e = begin_median_e,
    end_median_e = end_median_e
  )

  set.seed(NULL)


  ks_results
}
