test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

g <- BBSsize::granby

test_that("just_isd", {

  granby_isd <- just_isd(g, isd_seed = 1989)

  expect_true(nrow(granby_isd) == sum(g$abundance))
  expect_true(unique(granby_isd$isd_seed) == 1989)
  expect_true(all(granby_isd$mass > 0))

  granby_isd2 <- just_isd(g, isd_seed = 1977)
  granby_isd3 <- just_isd(g, isd_seed = 1977)

  expect_false(all(granby_isd$mass == granby_isd2$mass))

  expect_true(all(granby_isd2 == granby_isd3))

})

test_that("sampling_gmms" {

  sgmms <- construct_sampling_gmm(g, initial_isd_seed = 1989)
  sgmms2 <- construct_sampling_gmm(g, initial_isd_seed = 1977)

  sgmms3 <- construct_sampling_gmm(g, initial_isd_seed = 1977)

  sgmms4 <-  construct_sampling_gmm(g)

  expect_true(sum(sgmms$begin$density) == 1)
  expect_true(sum(sgmms$end$density) == 1)

  expect_true(all(sgmms$begin$mass == sgmms2$begin$mass))
  expect_false(all(sgmms$begin$density == sgmms2$begin$density))

  expect_true(all(sgmms2$begin == sgmms3$begin))

  expect_false(all(sgmms2$begin$density == sgmms4$begin$density))

})
