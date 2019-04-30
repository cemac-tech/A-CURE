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

The work is divided into sections.

1. The first section deals with the initial test of adding a namelist entry to the UM, and is described in notes made by HB [here](ACURE-Namelist-Test/Notes_by_HB.md)
2. The second section deals with how to create an ensemble of identical atmosphere tasks in a copy of a UKCA release job `u-bf737`, details of which are [here](Ensemble_Generation_1/Part1-Create_Ensemble.md)
3. The third section deals with how to introduce simple perturbations, hard coded on a single variable, details of which are [here](Ensemble_Generation_1/Part2-Add_A_Perturbation.md)
4. The fourth section expands this primitive dataframe system to a proper PPE style perturbation system, drawing on a csv file which provides values for the various PPE parameters. This was carried out based on a fully coupled UKESM suite ( either `u-bg935` or `u-bh079`). Due to conflicts with the reassigning of the {{DATAM}} parameter in the suite.rc file, this suite could only be run on newer versions of cylc such as is found on the new *pumatest* machine. Deatils of the changes made to the code to allow for this improved perturbed ensembling is found [here](Ensemble_Generation_2/Part3-Extended_to_dataframe.md). Accompanying files are included in the [folder](Ensemble_Generation_2) that gathers this aspect of the work
5. The fifth section implements these features into a copy of the AMIP suite `u-bi318`. Unexpectedly, the AMIP suite, on which the ACURE work will be based, uses vn11.1 of the UM. This required redoing the code changes to introduce a set of dummy variables and switches into the UKCA namelist. To make the work easier, and to limit the amount of work needed when the AMIP suite upgrades to vn11.3, a python script was written which prepares the necessary code changes so that they can be pasted into the relevant files. Extra changes were made to the suite.rc file to allow ensembling to be properly disabled. Those changes are currently being tested.
