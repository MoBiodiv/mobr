context("Testing of the rarefaction function")

test_that("mobr rarefaction function values are equivalent to true/accepted values", {
  
  # INDIVIDUAL BASED RAREFACTION TESTING (EXPECTATIONS 1 THROUGH 5)
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
  
  # SPACIAL-BASED RAREFACTION TESTING (EXPECTATIONS 6 THROUGH 9)
  # Compare mobr spacial-based rarefaction function against hand calculated values
  
  # Community matrix used for calculation
  comms = cbind(
    c(1, 1, 1, 0, 0, 0),
    c(0, 1, 1, 1, 0, 0),
    c(0, 0, 1, 1, 1, 0),
    c(0, 0, 0, 0, 1, 1),
    c(0, 0, 1, 0, 1, 0),
    c(0, 0, 0, 0, 1, 1),
    c(0, 0, 0, 0, 0, 1))
  
  # Calculate rarefaction values using mobr rarefaction function
  # Convert values to dataframe for easy referencing in test functions
  rarefaction_values = as.data.frame(rarefaction(comms, 'spat', coords = 1:6, latlong=F))
  
  # NOTE:
  # Rarefaction values to compare are at indeces: 1, 3, 5, and 6
  # Calculations for these values should remain constant with no deviation from the specified double/integer value
  
  # EXPECTATION 6
  # Check equivalence of rarefaction function generated value at index 1 to the decimal: 2.667 
  # Rounding was done to 3 decimal places
  expect_equal(round(rarefaction_values[1,1], digits = 3), 2.667)
  
  # EXPECTATION 7
  # Check equivalence of rarefaction function generated value at index 3 to the integer: 5 
  expect_equal(rarefaction_values[3,1], 5)
  
  # EXPECTATION 8
  # Check equivalence of rarefaction function generated value at index 5 to the double: 6.500 
  # Rounding was done to 3 decimal places
  expect_equal(round(rarefaction_values[5,1], digits = 3), 6.500)
  
  # EXPECTATION 9
  # Check equivalence of rarefaction function generated value at index 6 to the integer: 7 
  expect_equal(rarefaction_values[6,1], 7)

})
