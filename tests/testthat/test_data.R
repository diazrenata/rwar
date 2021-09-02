test_that("trivial", {

  expect_true(TRUE)


}
)
test_that("portal annual works", {

  portal_annual <- get_rodents_annual()

  expect_true(length(portal_annual) == 3)

  expect_true(names(portal_annual)[1] == "abundance")
  expect_true(names(portal_annual)[2] == "covariates")
  expect_true(names(portal_annual)[3] == "metadata")

  expect_true(is.data.frame(portal_annual$abundance))
  expect_true(is.data.frame(portal_annual$covariates))

  expect_equal(nrow(portal_annual$abundance), nrow(portal_annual$covariates))

  expect_true("year" %in% colnames(portal_annual$covariates))

  expect_equal(ncol(portal_annual$abundance), 21)
  expect_equal(ncol(portal_annual$covariates), 2)

  expect_false(anyNA(portal_annual$abundance))

  expect_true(all(portal_annual$covariates$year == c(1979:2018)))
}
)

test_that("subset_data_one works as expected", {

  portal_annual <- get_rodents_annual()

  subset_data_firstyear = subset_data_one(portal_annual, test_timestep = 1, buffer_size = 2)

  expect_true(length(subset_data_firstyear) == 5)
  expect_true(subset_data_firstyear$test$covariates$year == 1979)
  expect_true(subset_data_firstyear$train$covariates$year[1] == 1982)
  expect_true(nrow(subset_data_firstyear$train$abundance) == 37)

  subset_data_middleyear = subset_data_one(portal_annual, test_timestep = 10, buffer_size = 2)

  expect_true(subset_data_middleyear$test$covariates$year == 1988)
  expect_true(subset_data_middleyear$train$covariates$year[1] == 1979)
  expect_true(nrow(subset_data_middleyear$train$abundance) == 35)

  subset_data_lastyear = subset_data_one(portal_annual, test_timestep = 40, buffer_size = 2)

  expect_true(subset_data_lastyear$test$covariates$year == 2018)
  expect_true(subset_data_lastyear$train$covariates$year[1] == 1979)
  expect_true(nrow(subset_data_lastyear$train$abundance) == 37)


}
)

test_that("subset_data_all works as expected", {
  portal_annual <- get_rodents_annual()

  all_subsets <- subset_data_all(portal_annual, buffer_size = 2)

  expect_true(length(all_subsets) == 40)


  subset_data_firstyear = subset_data_one(portal_annual, test_timestep = 1, buffer_size = 2)

  for(i in 1:5) {
    expect_equal(all_subsets[[1]][[i]], subset_data_firstyear[[i]])

    }




  subset_data_middleyear = subset_data_one(portal_annual, test_timestep = 10, buffer_size = 2)
  for(i in 1:5) {
    expect_equal(all_subsets[[10]][[i]], subset_data_middleyear[[i]])

  }
}
)
