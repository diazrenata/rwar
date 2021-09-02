library(testthat)
library(cvlt)

abund_dat <- get_rodents_annual()

# Create a data subset

abund_dat_subset <- subset_data_one(abund_dat, 10, 2)

# Fit an LDA

lda_on_full <- LDA_set_user_seeds(
  document_term_table = abund_dat_subset$full$abundance,
  topics = 2,
  seed = 2)[[1]]

subsetted_lda <- subset_lda(lda_on_full, abund_dat_subset)


test_that("get abund probs works", {


  fitted_ts <- LDATS::TS_on_LDA(subsetted_lda,
                         document_covariate_table = as.data.frame(abund_dat_subset$train$covariates),
                         timename = "year",
                         formulas = ~1,
                         nchangepoints = 0,
                         control = LDATS::TS_control(nit = 100))[[1]]


  abund_probabilities <- get_abund_probabilities(abund_dat_subset, subsetted_lda, fitted_ts)

  betas <- exp(subsetted_lda@beta)
  thetas <- multinom_theta(abund_dat_subset, fitted_ts, sim = 10)

  one_probs <- thetas %*% betas

  expect_true(all(one_probs == abund_probabilities[[1]]))

  expect_true(all(dim(one_probs) == dim(abund_dat$abundance)))

})


test_that("one loglik works", {



  fitted_ts <- LDATS::TS_on_LDA(subsetted_lda,
                                document_covariate_table = as.data.frame(abund_dat_subset$train$covariates),
                                timename = "year",
                                formulas = ~1,
                                nchangepoints = 0,
                                control = LDATS::TS_control(nit = 100))[[1]]
  abund_probabilities <- get_abund_probabilities(abund_dat_subset, subsetted_lda, fitted_ts)

  one_probs <- abund_probabilities[[15]]

  test_row <- abund_dat_subset$test$abundance

  test_index <- abund_dat_subset$test_timestep

  prob_row <- one_probs[test_index, ]

  ll <- dmultinom(test_row, prob = prob_row, log = T)

  one_test_ll <- get_one_test_loglik(abund_dat_subset, one_probs)

  expect_true(ll == one_test_ll)

  })

test_that("all loglik work", {

  fitted_ts <- LDATS::TS_on_LDA(subsetted_lda,
                                document_covariate_table = as.data.frame(abund_dat_subset$train$covariates),
                                timename = "year",
                                formulas = ~1,
                                nchangepoints = 0,
                                control = LDATS::TS_control(nit = 100))[[1]]

  abund_probabilities <- get_abund_probabilities(abund_dat_subset, subsetted_lda, fitted_ts)
  one_probs <- abund_probabilities[[15]]
  one_test_ll <- get_one_test_loglik(abund_dat_subset, one_probs)

  all_ll <- get_test_loglik(abund_dat_subset, abund_probabilities)

  expect_true(all_ll[15] == one_test_ll)

}
)
