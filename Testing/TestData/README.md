# TestData

This directory contains data used by the Traffic Density Tool unit tests:

* **cellOrigTest.mat** - basic cell data structure containing one-year of data produced by `generateTestCellData`. This structure is modified by `generateTestData` to generate the various test data sets described below.
* **emptyCellCoverage.mat**, **emptyTest.mat** - used by the `testMissingData` unit test. In this test, every cell has 0% coverage/NaN data.
* **halfEmptyCellCoverage.mat**, **halfEmptyTest.mat** - used by the `testMissingData` unit test. In this test, every other cell has 0% coverage/NaN data.
* **job_test.dat** - contains meta-data associated with the input traffic density information used for testing: in particular, the size and extent of the input traffic density data (spatial and temporal bin definition).
* **testAirspaceA.mat**, **testAirspaceB.mat**, **testAirspaceC.mat**, **testAirspaceD.mat**, **testAirspaceO.mat**, **testAirspaceAll.mat** - airspace data used for testing. In testAirspaceAll.mat, all airspaces (A, B, C, D, O) have 20% coverage. In testAirspaceX.mat, Airspace X has 100% coverage; all other airspaces have 0% coverage.
* **testCellCoverage.mat** - cell coverage data used for testing: every cell has 100% coverage.
* **testData.mat** - simple fabricated data set containing one year's worth of data for 3 altitudes bins: 0-10000 ft, 10000-20000 ft, 20000-30000 ft. Cells in altitude bins 1, 2, 3 were assigned total time/maximum occupancy equal to 1, 2, 3, respectively. Data was NaN elsewhere. Used in `testGenerateDensity` unit test
* **testDataArea.mat** - variation of testData.mat that only contains positive data for latitudes between 48.5° and 49.5°, longitudes between -126.5° and -125.5° (0s elsewhere). Data is densest in the center becoming less dense toward the edges of the region. Used to test position-filtering capability in the `testArea` unit test
* **testDataDay.mat** - variation of testData.mat, where the total time/max occupancy of each cell is equal to the cell's day (e.g., 1 for Sunday, 2 for Monday, etc.). Used to test time-filtering capability in the `testTemporal` unit test
* **estDataHour.mat** - variation of testData.mat, where the total time/max occupancy of each cell is equal to the cell's hour (e.g., 1 for Hour Bin 1 (Midnight - 3AM), 2 for Hour Bin 2 (3AM - 6AM), etc.). Used to test time-filtering capability in the `testTemporal` unit test
* **testDataMonth.mat** - variation of testData.mat, where the total time/max occupancy of each cell is equal to the cell's month (e.g., 1 for January, 2 for February, etc.). Used to test time-filtering capability in the `testTemporal` unit test
* **testTrack.csv** - a straight-and-level test trajectory. Used to verify correct calculation of collision risk in the `testTrajectory` unit test
* **volumeMatrix.mat** - matrix containing the volume of each cell. Used to verify correct calculation of metrics in the `testGenerateDensity` and `testTemporal` unit tests

## <a name="diststatement"></a> Distribution Statement
DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

This material is based upon work supported by the National Aeronautics and Space Administration under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Aeronautics and Space Administration.

© 2019-2020 Massachusetts Institute of Technology.

Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
