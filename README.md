# [ACURE](https://gtr.ukri.org/projects?ref=NE%2FP013406%2F1) #

The Aerosol-Cloud Uncertainty REduction project

Ken Carslaw and Leighton Regayre UoL

<hr>

## Background information ##

Looking at running UKCA Perturbed Parameter Ensemble. Updating version of [UM](http://cms.ncas.ac.uk/wiki/UM) used. Previously they were working the [UMUI](http://cms.ncas.ac.uk/wiki/UM/RunningUMOnArcher) (vn8.4) and we're going
to move to running the [UM Rose Suites](http://cms.ncas.ac.uk/wiki/RoseCylc) v10.8 .

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
* A pair of ensemble rose suites which could be repurposed, one by david case with local copy u-bf163 and one by grenville with local copy u-al652-ens-l (see email)
