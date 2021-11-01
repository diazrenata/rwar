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

})

