# Unit Tests

To run all tests, run `RUN_unit_tests.m` in Matlab. This function runs the five tests described below. Test figures will be saved in a `Results` folder. At the end of tests, the number of tests that passed will be displayed.

### testGenerateDensity
The purpose of this test is to ensure that the traffic density tool 1) computes densities correctly and 2) properly accounts for airspace classes. This test reads in simple fabricated data where cells in altitude bins 1, 2, 3 are assigned total time/maximum occupancy 1, 2, 3, respectively. Data is NaN elsewhere. The test checks that density (both yearly averages and maximum observed occupancy) for both cooperative and noncooperative intruders has been properly computed. The output is separated into altitude layers and airspace classes.

### testTrajectory
The purpose of this test is to ensure that the traffic density tool correctly computes collision risk for an input test trajectory. The input trajectory is straight-and-level flight. The test checks that the collision rate for both cooperative and noncooperative intruders has been properly computed.

### testArea
The purpose of this test is to verify the traffic density tool's capability to filter based on position. This test uses a data set that only contains positive data for latitudes between 48.5° and 49.5°, and longitudes between -126.5° and -125.5°. (The data are 0s elsewhere). The observed occupancy is densest in the center becoming less dense toward the edges of the region. The test passes if the resulting plots show a concentric pattern with the greatest density in the center.

### testTemporal
The purpose of this test is to verify the traffic density tool's capability to filter based on time. This test uses data sets that are specifically designed to test the tool's capability to filter by month, day, and hour. See the TestData [README](./TestData/README.md) for more information about these data sets. In each test, the data is queried for a specific month, day, or hour. The test passes if the output only contains data for the specific times queried.

### testMissingData
The purpose of this test is to verify that the traffic density tool correctly outputs an estimate of the radar code coverage. This test is performed with two data sets: one where there is 0% coverage and one where there is 50% coverage. The test passes if the tool outputs the correct estimate of radar coverage corresponding to each data set.

**Note:** When the tests are run, warnings from `LoadData.m` will appear. This is because `LoadData.m` checks to ensure the traffic density tables containing 2018/2019 data have been properly loaded, whereas the unit tests use fabricated test data. These warnings are expected and can be ignored.

## <a name="diststatement"></a> Distribution Statement
DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

This material is based upon work supported by the National Aeronautics and Space Administration under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Aeronautics and Space Administration.

© 2019-2023 Massachusetts Institute of Technology.

Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
