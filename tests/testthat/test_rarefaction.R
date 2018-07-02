context("Testing of the rarefaction function")

test_that("Calcualted SAD values using rarefaction function are equivalent to Huberts (1971) paper values", {
  
  # Caluculate SAD values that should be equivalent to the SAD values provided in Table 2 of Stuart H. Huberts 
  # paper The Nonconcept on Species Diversity:A Critique and Alternative Parameters (1971)
  # SAD object numbers (i.e. sad1, sad2,...) correspond with the SAD values provided in Table 2
  
  # NOTE: sad3 values can be found in Table 1 of Howard L. Sanders paper Marine Benthic Diversity: A Comparative Study (1968)
  # Link: http://terascan.smast.umassd.edu/nasdata/archive/Bisagni_Data/students/abrunner/My%20Documents/Classes/ECOS630/papers/Sanders1968.pdf 

  # First SAD value
  sad1 = rep(10, 100)
  # Calculate rarefaction value and convert to dataframe for double value reference
  sad1DF = as.data.frame(rarefaction(sad1, 'indiv', 100))
  
  # Second SAD value
  sad2 = c(76, rep(50, 5), rep(20, 27-7), rep(5, 77-27), rep(1, 101-77))
  # Calculate rarefaction value and convert to dataframe for double value reference
  sad2DF = as.data.frame(rarefaction(sad2, 'indiv', 100))
  
  # Third SAD value
  sad3 = c(365, 112, 81, 61, 55, 46, 40, 38, 29, 23, 21, 15, 13, 12, 10, 8, 7, 7, 6, 6, 5, 5, 5, 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1)
  # Calculate rarefaction value and convert to dataframe for double value reference
  sad3DF = as.data.frame(rarefaction(sad3, 'indiv', 100))
  
  # Fourth SAD value
  sad4 = c(505, rep(5, 99))
  # Calculate rarefaction value and convert to dataframe for double value reference
  sad4DF = as.data.frame(rarefaction(sad4, 'indiv', 100))
  
  # Fifth SAD value
  sad5 = c(901, rep(1, 99))
  # Calculate rarefaction value and convert to dataframe for double value reference
  sad5DF = as.data.frame(rarefaction(sad5, 'indiv', 100))
  
  # EXPECTATION 1
  # Check equivalence of first SAD value with value from table - 65.3
  # Calculated rarefaction value has been rounded to nearest tenth to match table value
  expect_equal(round(sad1DF$`rarefaction(sad1, "indiv", 100)` , digits = 1), 65.3)
  
  # EXPECTATION 2
  # Check equivalence of second SAD value with value from table - 46.5
  # Calculated rarefaction value has been rounded to nearest tenth to match table value
  expect_equal(round(sad2DF$`rarefaction(sad2, "indiv", 100)`, digits = 1), 46.5)
  
  # EXPECTATION 3
  # Check equivalence of second SAD value with value from table - 20.4
  # Calculated rarefaction value has been rounded to nearest tenth to match table value
  expect_equal(round(sad3DF$`rarefaction(sad3, "indiv", 100)`, digits = 1), 20.4)
  
  # EXPECTATION 4
  # Check equivalence of second SAD value with value from table - 41.6
  # Calculated rarefaction value has been rounded to nearest tenth to match table value
  expect_equal(round(sad4DF$`rarefaction(sad4, "indiv", 100)`, digits = 1), 41.6)
  
  # EXPECTATION 5
  # Check equivalence of second SAD value with value from table - 10.9
  # Calculated rarefaction value has been rounded to nearest tenth to match table value
  expect_equal(round(sad5DF$`rarefaction(sad5, "indiv", 100)`, digits = 1), 10.9)
})
