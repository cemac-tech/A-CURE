# Converting a Rose Suite for Perturbed Parameter Ensemble (PPE) Experiments

To convert a rose suite for use in a PPE experiment, or for use in any general ensembling system, there are two main changes. The first is the creation of a python script to change the initial parameters for each run, and the inclusion of this script in the cylc task list. The second change is the creation of a parameterised task set within cylc to create the ensemble and keep the data output streams separate. The example I will be referring to in this document is the ACURE PPE suite [u-bi318], within which these changes have been made to convert the atmosphere-only UKESM AMIP suite to an ensembling suite.

## 1. Creation of the perturbation script

The perturbation script is designed to create a set of optional configuration files for the UM which can be used to overwrite the values of certain parameters. One file is created for each ensemble member with its own parameter set, with the parameters themselves read from a dataframe held in a csv file. Throughout the perturbation script there are various points at which the *acure* flag is expected as a part of the variable names. If adapting the script for use with a different experiment then these areas of code will need to be modified.

The perturbation script is called by cylc as a task prior to the fork of atmos_main jobs into the ensembled runs, and takes as its input a set of four arguments:
* The location of a template configuration file, exported with the name TMPL_LOC by Cylc. The template file is to contain every parameter to be perturbed as a logical set to false with the exception of the regional parameters, as each parameter has an associated logical in the um namelists.
* The location of the optional configuration files folder, typically `app/um/opt` in the Cylc run directory, exported with the name CONF_LOC by Cylc. This directory is where the output files from the script will be placed, for later reading by the atmos_main task in the UM run.
* The ensemble number list, referred to as ENS_LST by Cylc. This takes the form of a space separated list of ensemble members which can be itterated through to create the required specific optional configuration files.
* The location of the PPE dataframe. This dataframe is used by the script to provide the
PPE variables their values and logical switches. A discussion of the required layout of the dataframe file is included later.

In its function, the perturbation script iterates through the elements of the ensemble list. Each ensemble member corresponds to a row in the dataframe. For each ensemble member, the script reads the relevant row, and for variables that have a value other than "x" (used to indicate no perturbation) the replacement value is written to the output configuration file on the line following the corresponding logical switch, ie the section of the template file which reads:
```
l_acure_dry_dep_ait=.false.
l_acure_dry_dep_acc=.false.
l_acure_dry_dep_so2=.false.
l_acure_kappa_oc=.false.
l_acure_sig_w=.false.
```
would be combined with a line in the dataframe reading:
```
|acure_dry_dep_ait|acure_dry_dep_acc|acure_dry_dep_so2|acure_kappa_oc|acure_sig_w|
----------------------------------------------------------------------------------
|          x      |         0.26    |      17.54      |        x     |   0.001   |
```
to give an optional configuration file which contains the lines:
```
l_acure_dry_dep_ait=.false.
l_acure_dry_dep_acc=.true.
acure_dry_dep_acc=0.26
l_acure_dry_dep_so2=.true.
acure_dry_dep_so2=17.54
l_acure_kappa_oc=.false.
l_acure_sig_w=.true.
acure_sig_w=0.001
```

The parameter values written to this output file overwrite the corresponding parameters in the `app/um/rose-app.conf` file from which the UM obtains values input in the GUI - as such if ensembling is enabled the **values from the dataframe completely supercede those from the GUI**.

### Input File Formats ###

There are two input files for the perturbation script - the *template.conf* file and the *ppe_dataframe* file. Both of these files have format rules described below.

#### Template File ####

* The template file must have in it all the logical switches used for the perturbations.
* Typically this will be one per perturbed parameter, however for regional perturbations there is a one-to-many relationship.
* The logical switches must be divided and placed under a heading relating to the namelist in which the logical switch and the associated real variable reside within the code.
* All the logical switches must be set to `.false.`.

The listing below gives an example:
```
[namelist: run_ukca]
l_acure_bl_nuc=.false.
l_acure_ait_width=.false.
l_acure_cloud_ph=.false.
l_acure_carb_ff_ems=.false.
l_acure_carb_bb_ems=.false.
l_acure_carb_res_ems=.false.
l_acure_anth_so2=.false.
l_acure_carb_ff_diam=.false.
l_acure_carb_bb_diam=.false.


[namelist: run_radiation]
l_acure_bparam=.false.
l_acure_two_d_fsd_factor=.false.

[namelist: run_cloud]
l_acure_dbsdtbs_turb_0=.false.
```

#### PPE dataframe file

* The PPE dataframe file should be a csv file
* There should be one row per ensemble member
* If a parameter is not required as part of an ensemble members parameter set, the value in the dataframe should be "x", otherwise a real number is expected.
* There should be the same number of columns in the dataframe as there are logical switches in the template file, plus an additional 22 rows for the four regional parameters.
* Each column should have a header field with the name of the *real* parameter being perturbed. These parameters should match the name of the associated logical, omitting the `l_` characters, i.e. `l_acure_bl_nuc` in the template file would have a corresponding column in the dataframe with the header `acure_bl_nuc`. The columns need not be in the same order as in the template file however.


## 2. Modification of the cylc suite to allow ensembles

The cylc suite.rc file controls the way in which tasks are combined and chained to allow a full UM run to be completed. Generally this is done through constructing a cylc "graph" which, once compiled from the suite.rc file, may take the form
```
    [[dependencies]]
        [[[ R1 ]]]
            graph = fcm_make_um => fcm_make2_um => recon => atmos_main
        [[[ R1 ]]]
            graph = install_cold => install_ancil => recon => atmos_main
        [[[ R1 ]]]
            graph = fcm_make_pp => fcm_make2_pp => postproc
        [[[ R1 ]]]
            graph = fcm_make_pptransfer => fcm_make2_pptransfer => pptransfer
        [[[ P1D ]]]
            graph = atmos_main => postproc => pptransfer => housekeeping
```
Each of the graph entries is equivalent to a single branch of the cylc graph, shown in the picture below
![Cylc Graph](ppe_no_ens.png)
