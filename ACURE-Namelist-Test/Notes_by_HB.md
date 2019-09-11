# ACURE Task 1 - Name list change:

<hr>

Introduce namelist from leightons selection i.e. add to option mod. Initial ukcas_aeros_step consider the output channel and switches.

PRINT "IM NAMELIST X and I am switched on"

Extra print setting.

Then PRINT "IM NAMELIST X and I am switched on, my setting is Y"

*Consider: enhancing lookup values, look up tables, cataloguing member ID. Generating sub tasks*

Follow [coding standards](https://code.metoffice.gov.uk/doc/um/vn11.2/papers/umdp_003.pdf)

<hr>

Following the [namelist tutorial](https://code.metoffice.gov.uk/doc/um/latest/um-training/working-practices-tutorial.html)

From Variables.list I have chosen AEROS_VOLC_SO2.

## Open Ticket ##

Following [working practices](https://code.metoffice.gov.uk/doc/um/latest/um-training/working-practices.html)

* Summary: adding namelist items for ACURE
* type: task
* milestone: not for builds
* severity: minor (normally set by reviewer)
* keywords: macro, SC0138 (UKCA License agreement number)
* owner: self

ticket [#4653](https://code.metoffice.gov.uk/trac/um/ticket/4653#ticket)
created

## Checkout Branch ##

`fcm branch-create -k 4653 ACURE_Training fcm:um.x_tr@vn11.2`

include #ticketno in the log message which will automatically link in trac

`fcm checkout fcm:um.x_br/dev/christophersymonds/vn11.2_ACURE_Training`

** note multiple branch checkouts is diskspace intensive, all version controlled so no need to hold on to local copies!**

## namelist change ##

The namelist entry for aerosols will be in the UKCA

`vn11.2_ACURE_Training/src/atmosphere/UKCA/ukca_option_mod.F90`

* Added in a switch  l_ukca_aeros_volc_so2 default to FALSE
* And the value ukca_aeros_volc_so2 valid from 0 to 50
* A print statement prints l_ukca_aeros_volc_so2 = True/FALSE and
  ukca_aeros_volc_so2 = X alongside the rest of the namelist entries

## metadata change ##

`vn11.2_ACURE_Training/rose-meta/um-atmos/HEAD/rose-meta.conf`

add in

[namelist:run_ukca=l_ukca_aeros_volc_so2]
[namelist:run_ukca=ukca_aeros_volc_so2]

following the example of bio ems scaling parameter

* check metadata is valid: `rose metadata-check -C ./rose-meta/um-atmos/HEAD`
* reformat the metadata file so it is in the correct order:
  `rose config-dump -C ./rose-meta/um-atmos/HEAD`
* can edit
`rose config-edit -C ./rose-stem/app/um_ukca_eg -M ./rose-meta`
set meta-data version to be um-atmos/HEAD in the **um_ukca_eg**
and fix the red warning flags by adding to config.
* `fcm commit` with ticket number and mosrs username

## macros ##

As we have added a new variable, we need to write an upgrade macro to add it to all UM apps, so that when a user upgrades the new variable is present.

`rose-meta/um-atmos/versionXX_YY.py` (`version112_113.py `)

* `fcm commit` with ticket number

## Testing ##

*Changes which affect the Rose meta-data (either directly in the rose-meta.conf file or through an upgrade macro) need a test branch. This is because your code and meta-data changes will be committed to the trunk, but the resulting changes to rose stem should not be committed (the code reviewer will generate those based on the head of the trunk while merging your branch). To make sure that the effects of your tests are separate from your development branch, but are still under version control, you create a branch-of-branch from your development to contain the results of your meta-data changes.*


**NB: tutorial suggests rose stem test - they don't work on archer so skip that, might as well just test from dev branch and revert if necessary**

0. create test branch
  `fcm branch-create --branch-of-branch -k (ticket number) (new test branch name) -t test fcm:um.x_br/dev/(user)/(UM version)_(old dev branch name)`

  i.e.

  `fcm branch-create --branch-of-branch -k 4653 ACURE_TRAINING_TEST -t test fcm:um.x_br/dev/christophersymonds/vn11.2_ACURE-Training`

1. `./admin/rose-stem/update_all.py --path=$PWD --um=vn11.2_t4653`
  vn11.2_t4653 should match the AFTER_TAG in versionXX_YY.py
  * takes a while
2. If the macro is valid it will run with out errors
  * `fcm diff -g` will show changes in rose stem apps
  * errors will need to be fixed in both test and dev branches
  * can test in dev and use `fcm revert`
3. `fcm commit`

Note - due to a few mistakes in code, I needed to create 4 different test branches before the macro would run properly.

## Documenting ##

On trac add:

'''Ticket Summary: '''[wiki:ticket/<ticket number>/TicketSummary]

and create a ticket summary following ticket summary template

## Testing on in ukca release job ##

**u-bf766** (my rose suite of ukca release job 1 day)

```bash
cd ~/roses/u-bf766
rose edit -M ~/vn11.2_ACURE-Training/rose-meta
```
edit
fcm_make - env- sources + um source
+ /home/c.c.symonds/vn11.2_ACURE-Training
