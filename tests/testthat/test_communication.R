test_that("communication Rcpp structures can be copied to C++ structures",{
  control <- medfate::defaultControl()
  expect_type(medfate:::.testControlListToStructure(control), "integer")
})