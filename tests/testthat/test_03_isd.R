library(BBSsize)

h <- BBSsize::hartland
h_isd <- BBSsize::simulate_isd_ts(h, isd_seed = 1977)

test_that("isd_overlap works without years specified", {

  set.seed(1977)
  overlap <- compare_isds(h_isd$isd)

  begin_isd <- h_isd$isd %>%
    dplyr::filter(year %in% c(1994:1998))

  set.seed(1977)
  begin_gmm <- add_gmm(begin_isd)

  end_isd <- h_isd$isd %>%
    dplyr::filter(year %in% c(2014:2018))

  end_gmm <- add_gmm(end_isd)

  compare_densities <- dplyr::left_join(begin_gmm, end_gmm, by = c("mass")) %>%
    dplyr::group_by(mass) %>%
    dplyr::mutate(minDensity = min(density.x, density.y)) %>%
    dplyr::ungroup()

  expect_true(overlap$isd_turnover[1] == 1 - sum(compare_densities$minDensity))

}
)

test_that("isd_overlap works with years specified", {

  set.seed(1977)
  overlap <- compare_isds(h_isd$isd, begin_years = c(2000:2004), end_years = c(2010:2014))

  begin_isd <- h_isd$isd %>%
    dplyr::filter(year %in% c(2000:2004))

  set.seed(1977)
  begin_gmm <- add_gmm(begin_isd)

  end_isd <- h_isd$isd %>%
    dplyr::filter(year %in% c(2010:2014))

  end_gmm <- add_gmm(end_isd)

  compare_densities <- dplyr::left_join(begin_gmm, end_gmm, by = c("mass")) %>%
    dplyr::group_by(mass) %>%
    dplyr::mutate(minDensity = min(density.x, density.y)) %>%
    dplyr::ungroup()

  expect_true(overlap$isd_turnover[1] == 1 - sum(compare_densities$minDensity))

  expect_error(compare_isds(h_isd$isd, begin_years = c(2000:2004), end_years = c(1994:1998)))

}
)
