library(BBSsize)

h <- BBSsize::hartland

test_that("all core RUNS no years given", {

  all_res <- all_core_analyses(h)

}
)

test_that("all core RUNS years given", {

  all_res <- all_core_analyses(h, begin_years= c(2000:2004), end_years = c(2010:2014))


  expect_error(all_core_analyses(h, begin_years = c(2000:2004), end_years = c(1994:1998)))
}
)
