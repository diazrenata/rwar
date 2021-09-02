library(testthat)
library(rwar)
test_dir("testthat", reporter = c("check", "progress"))

test_check("rwar")
