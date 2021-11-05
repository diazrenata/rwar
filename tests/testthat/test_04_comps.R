test_that("trivial", {

  expect_true(TRUE)


}
)

library(BBSsize)

g <- BBSsize::granby

test_that("comps works", {


  comps <- be_comparison(g)


  expect_false(anyNA(comps))
  expect_true(nrow(comps) == 1)

})



