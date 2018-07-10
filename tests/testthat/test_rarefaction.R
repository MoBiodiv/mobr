context("Testing of the rarefaction function")

test_that("Test mobr rarefaction values against true/accepted values", {
  
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
  
  # SAMPLE-BASED RAREFACTION TESTING (EXPECTATIONS 6 THROUGH 9)
  # Compare mobr sample-based rarefaction function against the vegan:specaccum rarefaction function
  
  # Read in files for the next 4 expectations
  # Files are located in the testthat folder
  load("~/Mobr/tests/testthat/tank_comm.rda")
  load("~/Mobr/tests/testthat/inv_comm.rda")
  load("~/Mobr/tests/testthat/fire_comm.rda")
  coffee_comm = read.csv("~/Mobr/tests/testthat/coffee_comm.csv")
  
  # mobr sample-based rarefaction value calculation & set up
  # tank_comm
  mobrRareTank = as.data.frame(rarefaction(tank_comm, "samp"))
  mobrRareTankVal = mobrRareTank$`rarefaction(tank_comm, "samp")`[30]       
  
  # inv_comm
  mobrRareInv = as.data.frame(rarefaction(inv_comm, "samp"))
  mobrRareInvVal = mobrRareInv$`rarefaction(inv_comm, "samp")`[100] 
  
  # fire_comm
  mobrRareFire = as.data.frame(rarefaction(fire_comm, "samp"))
  mobrRareFireVal = mobrRareFire$`rarefaction(fire_comm, "samp")`[52]       
  
  # coffee_comm
  mobrRareCoffee = as.data.frame(rarefaction(coffee_comm[2:length(coffee_comm)], "samp"))
  mobrRareCoffeeVal = mobrRareCoffee$`rarefaction(coffee_comm[2:length(coffee_comm)], "samp")`[6]
  
  # Vegan::specaccum value calculation & set up
  # tank_comm
  specRareTank = vegan::specaccum(tank_comm, method = "rarefaction")
  specRareTankVal = specRareTank[[4]][30]
  
  # inv_comm
  specRareInv = vegan::specaccum(inv_comm, method = "rarefaction")
  specRareInvVal = specRareInv[[4]][100]
  
  # fire_comm
  specRareFire = vegan::specaccum(fire_comm, method = "rarefaction")
  specRareFireVal = specRareFire[[4]][52]
  
  # coffee_comm
  specRareCoffee = vegan::specaccum(coffee_comm[2:length(coffee_comm)], method = "rarefaction")
  specRareCoffeeVal = specRareCoffee[[4]][6]
  
  # EXPECTATION 6
  # Check equivalence of LAST richness value given for mobr and vegan::specaccum rarefaction function
  # tank_comm richness value for last (30th) value should equal 47.000
  # value has been rounded to 3 decimal points
  expect_equal(round(mobrRareTankVal, digits = 3), 
               round(specRareTankVal, digits = 3))
  
  # EXPECTATION 7
  # Check equivalence of LAST richness value given for mobr and vegan::specaccum rarefaction function
  # inv_comm richness value for last (100th) value should equal 111.000
  # value has been rounded to 3 decimal points
  expect_equal(round(mobrRareInvVal, digits = 3), 
               round(specRareInvVal, digits = 3))
  
  # EXPECTATION 8
  # Check equivalence of LAST richness value given for mobr and vegan::specaccum rarefaction function
  # tank_comm richness value for last (52nd) value should equal 21.000
  # value has been rounded to 3 decimal points
  expect_equal(round(mobrRareFireVal, digits = 3), 
               round(specRareFireVal, digits = 3))
  
  # EXPECTATION 9
  # Check equivalence of LAST richness value given for mobr and vegan::specaccum rarefaction function
  # inv_comm richness value for last (6th) value should equal 21.000
  # value has been rounded to 3 decimal points
  expect_equal(round(mobrRareCoffeeVal, digits = 3), 
               round(specRareCoffeeVal, digits = 3))


})
