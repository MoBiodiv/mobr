### NOTE: the line values for the Fire_lat_longs datafile are NOT including the header row
### row 1 is considered to be the first line with lat/loing values

context("Testing of the sphere_dist function")

test_that("sphere distance function values are equivalent to known true values", {

  # Matrix of first 5 long/lat values for the Fire_lat_long dataset
  firstFive_fire_lat_longs = matrix(c(-92.75974,
                                      -92.76202,
                                      -92.76155,
                                      -92.76958,
                                      -92.76854,
                                      37.97213,
                                      37.97941,
                                      37.97928,
                                      37.96129,
                                      37.96147), 
                                      nrow = 5, ncol = 2)
  
  # Distance matrix for first five Fire_lat_long correct values
  firstFive_dist = matrix(c(0.0, 0.0001327012, 0.0001284380, 0.0002514478, 0.0002378107,
                            0.0001327012, 0.0, 0.0000082031, 0.0003408448, 0.0003317490,
                            0.0001284380, 0.0000082031, 0.0, 0.0003417858, 0.0003323211,
                            0.0002514478, 0.0003408448, 0.0003417858, 0.0, 0.0000178100,
                            0.0002378107, 0.0003317490, 0.0003323211, 0.0000178100, 0.0))

  
  # EXPECTATION 1
  # Compute sphere_dist using first five Fire_lat_longs data values
  # Compare it to the correct dist matrix above - firstFive_dist
  expect_equal(as.list(round(sphere_dist(firstFive_fire_lat_longs), 10)), as.list(firstFive_dist))
  
  
  # Matrix of Fire_lat_long data rows 20 to 24
  # in Long/Lat format for function to work properly
  SecondSet_fire_lat_longs = matrix(c(-92.5877,
                                      -92.58775,
                                      -92.58838,
                                      -92.58806,
                                      -92.58658,
                                      38.05347,
                                      38.06267,
                                      38.06311,
                                      38.06386,
                                      38.06178), 
                                      nrow = 5, ncol = 2)
  
  SecondSet_dist = matrix(c(0.0, 0.0001605724, 0.0001686162, 0.0001814351, 0.0001461874,
                            0.0001605724, 0.0, 0.0000128285, 0.0000213759, 0.0000246104,
                            0.0001686162, 0.0000128285, 0.0, 0.0000140918, 0.0000374318,
                            0.0001814351, 0.0000213759, 0.0000140918, 0.0, 0.0000435969,
                            0.0001461874, 0.0000246104, 0.0000374318, 0.0000435969, 0.0))

  # EXPECTATION 2
  # Compute sphere_dist using Fire_lat_longs data values 20 - 24
  # Compare it to the correct dist matrix above - SecondSet_dist
  expect_equal(as.list(round(sphere_dist(SecondSet_fire_lat_longs), 10)), as.list(SecondSet_dist))
  
  # Matrix of Fire_lat_long data rows 40 to 44
  # in Long/Lat format for function to work properly
  ThirdSet_fire_lat_longs = matrix(c(-92.59498,
                                      -92.61025,
                                      -92.61071,
                                      -92.61199,
                                      -92.61252,
                                      38.11853,
                                      38.06127,
                                      38.06185,
                                      38.0618,
                                      38.06209), 
                                      nrow = 5, ncol = 2)
  
  ThirdSet_dist = matrix(c(0.0, 0.0010292792, 0.0010212643, 0.0010274289, 0.0010248742,
                           0.0010292792, 0.0, 0.0000126027, 0.0000298647, 0.0000397119,
                           0.0010212643, 0.0000126027, 0.0, 0.0000209049, 0.0000298291,
                           0.0010274289, 0.0000298647, 0.0000209049, 0.0, 0.0000100203,
                           0.0010248742, 0.0000397119, 0.0000298291, 0.0000100203, 0.0))
  
  # EXPECTATION 3
  # Compute sphere_dist using Fire_lat_longs data values 40 - 44
  # Compare it to the correct dist matrix above - ThirdSet_dist
  expect_equal(as.list(round(sphere_dist(ThirdSet_fire_lat_longs), 10)), as.list(ThirdSet_dist))
})


