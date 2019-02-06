# Work Log for A-CURE

## Current state as of 2019/02/04

Work has been done to familiarise self with the UM by following the NCAS tutorial and the UKCA tutorial. This has covered Stash additions, rose suites, cylc and branch creation. A section on testing with rose-stem had to be skipped since it will mainly only work on met office systems.
A test branch of UM vn11.2 has been created with the variable `REAL :: ukca_aeros_volc_so4` added to the UKCA system. The procedure for adding a variable, and getting it to display in the namelist, is shown in the [UKCA tutorial](http://www.ukca.ac.uk/wiki/index.php/UKCA_Chemistry_and_Aerosol_Tutorials_at_vn10.9) part 11.
The test branch has been given the name "vn11.2_ACURE-Training" and can be accessed from `fcm:um.x_br/dev/christophersymonds/vn11.2_ACURE-Training`. This branch uses the new _aeros_ variable and gives it an initial value.
An attempt at compilation is needed, as is investigation of the ensemble rose suite. Possible use of UKCA release jobs for running, once those have been examined.
