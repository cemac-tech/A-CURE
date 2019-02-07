# Work Log for A-CURE

## Current state as of 2019/02/04

Work has been done to familiarise self with the UM by following the NCAS tutorial and the UKCA tutorial. This has covered Stash additions, rose suites, cylc and branch creation. A section on testing with rose-stem had to be skipped since it will mainly only work on met office systems.

A test branch of UM vn11.2 has been created with the variable `REAL :: ukca_aeros_volc_so4` added to the UKCA system. The procedure for adding a variable, and getting it to display in the namelist, is shown in the [UKCA tutorial](http://www.ukca.ac.uk/wiki/index.php/UKCA_Chemistry_and_Aerosol_Tutorials_at_vn10.9) part 11. Details of the process were written up by HB and a copy of this write-up is included in the ACURE Task 1 folder

The test branch has been given the name "vn11.2_ACURE-Training" and can be accessed from `fcm:um.x_br/dev/christophersymonds/vn11.2_ACURE-Training`. This branch uses the new _aeros_ variable and gives it an initial value.

An attempt at compilation is needed, as is investigation of the ensemble rose suite. Possible use of UKCA release jobs for running, once those have been examined.

## 2019/02/07

### Rose suites of particular use:


|Rose Suite|Description|
|-|-|
|u-bf163|A copy of the UKESM 1.0 beta suite with ensembles implemented|
|u-bf737|A copy of the 1-day UKCA prototype release job for UKESM (also works with vn11.2). Set up to run the UM with ACURE variable `ukca_aeros_volc_so2`

The rose suite from the UKESM 1day proposed release job was copied, and should be run with the command

`cd ~/roses/u-bf737
rose edit -M ~/um/branches/vn11.2_ACURE-Training/rose-meta/`

and the meta field in um should be set to `um-atmos/HEAD`

This enables the modified meta data to be used. The current states of the rose suite and the code branch have both been committed to the repository. The branch and rose suite have been tested, and run correctly, showing the existance of the ACURE variable in the namelist printout.

### The Ensemble Rose Suite

Work will continue this afternoon on understanding the ensemble rose suite made by Dave Sexton (u-bf163). One quesation which needs answering is in what way do the required ensemble perturbations differ from those used by Dave Sexton?

The Ensembles seem to be controlled through modification of the suite.rc file. This sets up cylc cycles with different ensemble members in each. This is an idea to investigate further.
