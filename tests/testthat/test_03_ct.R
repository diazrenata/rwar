test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

g <- BBSsize::granby
h <- BBSsize::hartland

test_that("ct works", {

  not_g <- core_transient(g, F)
  core_g <- core_transient(g, T)

  g_manual <- data.frame(
    sp = colnames(g$abundance),
    occur = colSums(g$abundance > 0),
    nobs  = nrow(g$abundance)
  ) %>%
    dplyr::mutate(prop_obs = occur/nobs) %>%
    dplyr::mutate(is_transient = prop_obs < 1/3,
                  is_core = prop_obs > 2/3)

  expect_true(ncol(not_g$abundance) == sum(!(g_manual$is_transient)))
  expect_true(all(colnames(not_g$abundance) == dplyr::filter(g_manual, !is_transient)$sp))
  expect_true(all(sort(colnames(not_g$abundance)) == sort(not_g$metadata$species_table$id)))

  expect_true(ncol(core_g$abundance) == sum((g_manual$is_core)))
  expect_true(all(colnames(core_g$abundance) == dplyr::filter(g_manual, is_core)$sp))
  expect_true(all(sort(colnames(core_g$abundance)) == sort(core_g$metadata$species_table$id)))

  not_h <- core_transient(h, F)
  core_h <- core_transient(h, T)

  h_manual <- data.frame(
    sp = colnames(h$abundance),
    occur = colSums(h$abundance > 0),
    nobs  = nrow(h$abundance)
  ) %>%
    dplyr::mutate(prop_obs = occur/nobs) %>%
    dplyr::mutate(is_transient = prop_obs < 1/3,
                  is_core = prop_obs > 2/3)

  expect_true(ncol(not_h$abundance) == sum(!(h_manual$is_transient)))
  expect_true(all(colnames(not_h$abundance) == dplyr::filter(h_manual, !is_transient)$sp))
  expect_true(all(sort(colnames(not_h$abundance)) == sort(not_h$metadata$species_table$id)))

  expect_true(ncol(core_h$abundance) == sum((h_manual$is_core)))
  expect_true(all(colnames(core_h$abundance) == dplyr::filter(h_manual, is_core)$sp))
  expect_true(all(sort(colnames(core_h$abundance)) == sort(core_h$metadata$species_table$id)))

  expect_true(ncol(core_h$abundance) < ncol(not_h$abundance))
  expect_true(ncol(not_h$abundance) < ncol(h$abundance))


})


test_that("jt works", {

  not_g <- core_transient(g, F)
  core_g <- core_transient(g, T)

  g_manual <- data.frame(
    sp = colnames(g$abundance),
    occur = colSums(g$abundance > 0),
    nobs  = nrow(g$abundance)
  ) %>%
    dplyr::mutate(prop_obs = occur/nobs) %>%
    dplyr::mutate(is_transient = prop_obs < 1/3,
                  is_core = prop_obs > 2/3)

  expect_true(ncol(not_g$abundance) == sum(!(g_manual$is_transient)))
  expect_true(all(colnames(not_g$abundance) == dplyr::filter(g_manual, !is_transient)$sp))
  expect_true(all(sort(colnames(not_g$abundance)) == sort(not_g$metadata$species_table$id)))

  expect_true(ncol(core_g$abundance) == sum((g_manual$is_core)))
  expect_true(all(colnames(core_g$abundance) == dplyr::filter(g_manual, is_core)$sp))
  expect_true(all(sort(colnames(core_g$abundance)) == sort(core_g$metadata$species_table$id)))

  not_h <- core_transient(h, F)
  core_h <- core_transient(h, T)

  h_manual <- data.frame(
    sp = colnames(h$abundance),
    occur = colSums(h$abundance > 0),
    nobs  = nrow(h$abundance)
  ) %>%
    dplyr::mutate(prop_obs = occur/nobs) %>%
    dplyr::mutate(is_transient = prop_obs < 1/3,
                  is_core = prop_obs > 2/3)

  expect_true(ncol(not_h$abundance) == sum(!(h_manual$is_transient)))
  expect_true(all(colnames(not_h$abundance) == dplyr::filter(h_manual, !is_transient)$sp))
  expect_true(all(sort(colnames(not_h$abundance)) == sort(not_h$metadata$species_table$id)))

  expect_true(ncol(core_h$abundance) == sum((h_manual$is_core)))
  expect_true(all(colnames(core_h$abundance) == dplyr::filter(h_manual, is_core)$sp))
  expect_true(all(sort(colnames(core_h$abundance)) == sort(core_h$metadata$species_table$id)))

  expect_true(ncol(core_h$abundance) < ncol(not_h$abundance))
  expect_true(ncol(not_h$abundance) < ncol(h$abundance))


})

