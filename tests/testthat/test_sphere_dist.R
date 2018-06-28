### NOTE: the Fire_lat_longs data file from the mobr data folder must be uploaded and in ones current working directory for
### the following code to run since the file is used in testing.

context("Testing of the sphere_dist function")

test_that("matrix values are accurate provided 2 columns of lat and longs data", {
  
  # Remove first column of 'plotunique' data from the Fire_lat_longs data file
  test_fire_lat_longs = Fire_lat_longs[ ,2:3]
  
  # Compute sphere_dist matrix to be tested against- assume correct values are calcuated
  test_dist1 = sphere_dist(test_fire_lat_longs)
  
  # EXPECTATION 1
  # Compute sphere_dist using Fire_lat_longs data and compare it to the test_dist1 matrix created above
  # This test should always pass since both distance matrices are calculated from the same data
  expect_equal(sphere_dist(test_fire_lat_longs), test_dist1)
  
})
