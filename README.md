# [ACURE](https://gtr.ukri.org/projects?ref=NE%2FP013406%2F1) #

The Aerosol-Cloud Uncertainty REduction project

Ken Carslaw and Leighton Regayre UoL

<hr>

## Background information ##

Looking at running UKCA Perturbed Parameter Ensemble. Updating version of [UM](http://cms.ncas.ac.uk/wiki/UM) used. Previously they were working the [UMUI](http://cms.ncas.ac.uk/wiki/UM/RunningUMOnArcher) (vn8.4) and we're going
to move to running the [UM Rose Suites](http://cms.ncas.ac.uk/wiki/RoseCylc) v11.1 .

A PPE is a perturbed parameter ensemble, so rather than perturbing a
start file, each member will have different values for some (maybe
only one) parameters (same start file).

Detailed UM documentation can be found [here](https://code.metoffice.gov.uk/doc/um/latest/umdp.html) (login required)

[CEMAC's Role is outlined in the pdf](CEMAC_priorities_for_ACURE_v2_1-1.pdf)

<hr>

## Deliverable

A set of rose suites and a branch of the UKESM or UM which allows ensembles to be run

## Tasks

Tasks to produce the deliverable will include:

### Preliminary test
* Create a branch of the UKESM which contains a single new variable in the UKCA section
* Implement some science behind this new variable
* Create an ensemble rose suite which can perturb the variable

### Full system
* As above but with a full set of 33 variables, some of which are new, others of which are already in the UKCA module, and some which would replace existing variables and give new, slightly modified, behaviour

## Resources

Currently, various resources have been used. Among these are:

* The UKCA tutorials which instruct on the addition of variables to the UKCA section of the UM
* A pair of ensemble rose suites which could be repurposed, one by david case with local copy u-bf163 and one by grenville with local copy on PUMA u-al652-ens-l (see email)

## Progress

The work is divided into sections. Each of these sections has a separate readme, linked to here. Brief descriptions of the different sections are listed below.

1. [Namelist Additions](ACURE-Namelist-Test/Notes_by_HB.md)
2. [Simple Ensemble Creation](Ensemble_Generation_1/Part1-Create_Ensemble.md)
3. [Simple Perturbation Creation](Ensemble_Generation_1/Part2-Add_A_Perturbation.md)
4. [Perturbations from a Dataframe](Ensemble_Generation_2/Part3-Extended_to_dataframe.md)
5. [Incorporation in the AMIP suite](Ensemble_Generation_3/Part4-Implement_In_AMIP.md)
6. [Addition of Regional Perturbations](Regional_Perturbation/Regional_Perturbations.md)
7. [Expansion of PPE to special cases](Regional_Perturbation/Special_Case_Parameters.md)


[The first section](ACURE-Namelist-Test/Notes_by_HB.md) deals with the initial test of adding a namelist entry to the UM, and is described in notes made by Helen Burns.<br>
[The second section](Ensemble_Generation_1/Part1-Create_Ensemble.md) deals with how to create an ensemble of identical atmosphere tasks in a copy of a UKCA release job `u-bf737`.<br>
[The third section](Ensemble_Generation_1/Part2-Add_A_Perturbation.md) deals with how to introduce simple perturbations, hard coded on a single variable.<br>
[The fourth section](Ensemble_Generation_2/Part3-Extended_to_dataframe.md) expands this primitive dataframe system to a proper PPE style perturbation system, drawing on a csv file which provides values for the various PPE parameters. This was carried out based on a fully coupled UKESM suite ( either `u-bg935` or `u-bh079`). Due to conflicts with the reassigning of the {{DATAM}} parameter in the suite.rc file, this suite could only be run on newer versions of cylc such as is found on the new *pumatest* machine. Accompanying files are included in the [folder](Ensemble_Generation_2) that gathers this aspect of the work<br>
[The fifth section](Ensemble_Generation_3/Part4-Implement_In_AMIP.md) implements these features into a copy of the AMIP suite `u-bi318`. Unexpectedly, the AMIP suite, on which the ACURE work will be based, uses vn11.1 of the UM. This required redoing the code changes to introduce a set of dummy variables and switches into the UKCA namelist. To make the work easier, and to limit the amount of work needed when the AMIP suite upgrades to vn11.3, a python script was written which prepares the necessary code changes so that they can be pasted into the relevant files. Extra changes were made to the suite.rc file to allow ensembling to be properly disabled. Further changes were made to ensure that the post processing of the output files was done separately by making these parameterised tasks just like the atmos_main task<br>
[The sixth section](Regional_Perturbation/Regional_Perturbations.md) deals with the expansion of the PPE system to include regional perturbations. These perturbations will read a set of five or six perturbation parameters for the four parameters `anth_so2`, `carb_bb_ems`, `carb_ff_ems`, `carb_res_ems`, each of these five or six parameters relating to a single region (ie `acure_anth_so2_region_1`, `acure_anth_so2_region_2`, `acure_anth_so2_region_3`, `acure_anth_so2_region_4`,`acure_anth_so2_region_5`). A "regional mask", modelled on the land-sea mask, is created for each of the regionally perturbed parameters which is read in as an ancillary file and then referenced in the code to apply the regional perturbations to the correct grid cells.<br>
[The seventh section](Regional_Perturbation/Special_Case_Parameters.md) is a short section dealing with the inclusion of special case parameters. These parameters work differently in the perturbation scheme to the standard UKCA parameters, which come in logical-real pairs. The regional perturbations for example have 5 or 6 real parameters controlled by a single logical. Parameters that are perturbed directly such as bparam or c_r_correl have a "l_acure_*" format logical in their own namelist, but no corresponding extra parameter. One parameter, the pcalc_index parameter, changes a filename pointer rather than adding a real parameter.<br>

## Summary of Conversion of a Rose Suite for Ensembling

A [summary](rose_suite_doc.md) of the work needed to convert a rose suite to an ensembling PPE suite has been produced which combines all the aspects of the modification. This was needed as the details given in the earlier documents all depend on the previous document, and so none of them are "stand alone" documents. As such, the summary document covers all aspects of converting a general suite into a general ensembling suite. To fully understand all the work done to produce the final ACURE rose suite `u-bi318`, this document should be read in combination with the documents on [regional perturbations](Regional_Perturbation/Regional_Perturbations.md) and on [special parameters](Regional_Perturbation/Special_Case_Parameters.md).
