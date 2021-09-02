library(BBSsize)

h <- BBSsize::hartland

test_that("composition works no years given", {

  spcomp <- compare_species_composition(h)

  start_comp <- h$abundance[1:5, ]
  end_comp <- h$abundance[21:25, ]

  start_totals <- colSums(start_comp)
  start_relatives <- start_totals / sum(start_totals)

  end_totals <- colSums(end_comp)
  end_relatives <- end_totals / sum(end_totals)

  rel_comparison <- data.frame(
    species1 = names(start_relatives),
    startRel = start_relatives,
    species2 = names(end_relatives),
    endRel = end_relatives)

  rel_comparison <- rel_comparison %>%
    dplyr::group_by(species1) %>%
    dplyr::mutate(minRel = min(startRel, endRel)) %>%
    dplyr::ungroup()

  expect_true(all(rel_comparison$species1 == rel_comparison$species2))

  expect_true(spcomp$sp_turnover == 1 - sum(rel_comparison$minRel))


  compmatrix <- matrix(nrow = 2, data = c(start_totals, end_totals), byrow = T)

  expect_true(all(compmatrix[2, ] == end_totals))
  expect_true(all(compmatrix[1, ] == start_totals))

  bcd <- vegan::vegdist(compmatrix)

  expect_true(bcd == spcomp$bcd)

})

test_that("composition works with years given", {

  spcomp <- compare_species_composition(h, begin_years = c(2000:2004), end_years = c(2010:2014))

  start_comp <- h$abundance[7:11, ]
  end_comp <- h$abundance[17:21, ]

  start_totals <- colSums(start_comp)
  start_relatives <- start_totals / sum(start_totals)

  end_totals <- colSums(end_comp)
  end_relatives <- end_totals / sum(end_totals)

  rel_comparison <- data.frame(
    species1 = names(start_relatives),
    startRel = start_relatives,
    species2 = names(end_relatives),
    endRel = end_relatives)

  rel_comparison <- rel_comparison %>%
    dplyr::group_by(species1) %>%
    dplyr::mutate(minRel = min(startRel, endRel)) %>%
    dplyr::ungroup()

  expect_true(all(rel_comparison$species1 == rel_comparison$species2))

  expect_true(spcomp$sp_turnover == 1 - sum(rel_comparison$minRel))


  compmatrix <- matrix(nrow = 2, data = c(start_totals, end_totals), byrow = T)

  expect_true(all(compmatrix[2, ] == end_totals))
  expect_true(all(compmatrix[1, ] == start_totals))

  bcd <- vegan::vegdist(compmatrix)

  expect_true(bcd == spcomp$bcd)

})
