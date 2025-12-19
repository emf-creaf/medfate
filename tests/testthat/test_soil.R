library(medfate)



test_that("default soils can be created",{
  expect_s3_class(defaultSoilParams(), "data.frame")
  df <- defaultSoilParams()
  expect_s3_class(soil(df), "soil")
  df_minimal <- df
  df_minimal$om <- NULL
  df_minimal$ph <- NULL
  df_minimal$nitrogen <- NULL
  expect_s3_class(soil(df_minimal), "soil")
})

test_that("soils can be redefined",{
  out_widths <- c(1000,2000)
  df <- defaultSoilParams()
  df_minimal <- df
  df_minimal$om <- NULL
  df_minimal$ph <- NULL
  df_minimal$nitrogen <- NULL
  expect_s3_class(soil_redefineLayers(df, widths = out_widths), "data.frame")
  expect_s3_class(soil_redefineLayers(soil(df), widths = out_widths), "soil")
  expect_s3_class(soil_redefineLayers(df_minimal, widths = out_widths), "data.frame")
  expect_s3_class(soil_redefineLayers(soil(df_minimal), widths = out_widths), "soil")
})

test_that("soils can be summarized",{
  out_widths <- c(1000,2000)
  df <- defaultSoilParams()
  df_minimal <- df
  df_minimal$om <- NULL
  df_minimal$ph <- NULL
  df_minimal$nitrogen <- NULL
  s_minimal <- soil(df_minimal)
  # s_minimal$om <- NULL
  s_minimal$ph <- NULL
  s_minimal$nitrogen <- NULL
  expect_type(capture.output(summary(soil(df), widths = out_widths)), "character")
  expect_type(capture.output(summary(soil(df_minimal), widths = out_widths)), "character")
  expect_type(capture.output(summary(s_minimal, widths = out_widths)), "character")
})