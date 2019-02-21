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
|u-bf737|A copy of the 1-day UKCA prototype release job for UKESM (also works with vn11.2). Set up to run the UM with ACURE variable `ukca_aeros_volc_so2`. Used for further development of including ensembles.|
|u-bf766|A copy of u-bf737 made after ACURE namelist changes were confirmed to be working but no ensemble changes had been made. Good for rolling back.|

The rose suite from the UKESM 1day proposed release job was copied, and should be run with the command

`cd ~/roses/u-bf737
rose edit -M ~/um/branches/vn11.2_ACURE-Training/rose-meta/`

and the meta field in um should be set to `um-atmos/HEAD`

This enables the modified meta data to be used. The current states of the rose suite and the code branch have both been committed to the repository. The branch and rose suite have been tested, and run correctly, showing the existance of the ACURE variable in the namelist printout.

**I need to figure out exactly how to make the metadata from a branch the default, so that the rose suite can be used on it's own
drawing from the remote branch on the repository**

### The Ensemble Rose Suite

* Work continues on understanding the ensemble rose suite made by Dave Sexton (u-bf163). One quesation which needs answering is in what way do the required ensemble perturbations differ from those used by Dave Sexton?

* The Ensembles seem to be controlled through modification of the suite.rc file. This sets up cylc cycles with different ensemble members in each. This is an idea to investigate further.

* A problem with the original suite u-bd149 is that for some reason the archer username has been hard coded to be dcase (that would be David Case). In addition to this the variable $ARCHER_USERNAME had been set to be the same as the login name for puma

* A further problem, this one being more general, is that the username I have on puma, `c.c.symonds` causes a Jinja2 error due to the dots. The error caused for this has the form 

```[FAIL] cylc validate -v --strict u-bf163 # return-code=1, stderr=  
[FAIL] Jinja2Error:
[FAIL]   File "/usr/local/python/lib/python2.6/site-packages/Jinja2-2.7.3-py2.6.egg/jinja2/environment.py", 
                line 397, in getattr
[FAIL]     return getattr(obj, attribute)
[FAIL] UndefinedError: 'c' is undefined
```
* The suite as it would run is held in `cylc-run/u-bf163/suite.rc.processed`. Comparison against the original suite.rc file should be very illuminating.

## 2019/02/11

### Setting up an Ensemble

To use ensembles, you must use parameterisation in the suite configuration file `suite.rc`

Changes were made to the `u-bf737` suite, with the aim of implementing ensembles. Currently the ensembles are generated but cannot run to completion. The reason for this will be investigated shortly. The output of fcm diff is contained below

```ini
Index: app/fcm_make_um/rose-app.conf
===================================================================
--- app/fcm_make_um/rose-app.conf       (revision 105592)
+++ app/fcm_make_um/rose-app.conf       (working copy)
@@ -43,4 +43,4 @@
 stash_version=1A
 timer_version=3A
 um_rev=vn11.2
-um_sources=fcm:um.x_br/dev/christophersymonds/vn11.2_ACURE-Training
+um_sources=branches/dev/christophersymonds/vn11.2_ACURE-Training
Index: app/um/rose-app.conf
===================================================================
--- app/um/rose-app.conf        (revision 105592)
+++ app/um/rose-app.conf        (working copy)
@@ -15,7 +15,7 @@
 !!DR_HOOK_PROFILE=$CYLC_TASK_WORK_DIR/drhook.prof.%d
 !!DR_HOOK_PROFILE_LIMIT=-10.0
 !!DR_HOOK_PROFILE_PROC=-1
-ENS_MEMBER=0
+ENS_MEMBER=${ENS_MEMBER}          ######### This assignment uses a valuie exported from cylc
 HISTORY=$DATAM/${RUNID}.xhist     ######### The ENS_MEMBER variable in um/rose-app.conf is always present
 PRINT_STATUS=PrStatus_Normal
 RCF_PRINTSTATUS=PrStatus_Normal
Index: meta/rose-meta.conf
===================================================================
--- meta/rose-meta.conf (revision 105592)
+++ meta/rose-meta.conf (working copy)
@@ -194,6 +194,20 @@
 title=Model output data and restart directory
 type=character
 
 ######################################################## Addition of A-CURE Ensembles to the GUI
 
+[jinja2:suite.rc=ENSEMBLE_SIZE]
+compulsory=true
+description=Size of the A-CURE ensemble
+help=This is the ACURE ensemble variable
+    =
+    =This variable is currently only present for testing the rose suite,
+    =and does not have any impact on the values of parameters
+ns=cycle
+range=1:50
+sort-key=999
+title=A-CURE Ensemble Size
+type=integer
+
+
 [jinja2:suite.rc=EXTRACT_HOST]
 compulsory=true
 description=Host to use for code extraction
Index: rose-suite.conf
===================================================================
--- rose-suite.conf     (revision 105592)
+++ rose-suite.conf     (working copy)
@@ -19,6 +19,7 @@
 CLOCK='PT10M'
 !!CORE='broadwell'
 DATAM='History_Data'
+ENSEMBLE_SIZE=3                     ############### Default ensemble size
 !!EXTRACT_HOST='linux'
 !!FUNDING='hccp'
 !!FUNDING_OTHER=''
@@ -53,7 +54,7 @@
 RUNID=true
 !!RUNID_USR=''
 RUNLEN='other'
-RUNLEN_OTHER='P1D'
+RUNLEN_OTHER='P3D'
 SITE='archer'
 !!SUBPROJECT='esms'
 !!SUBPROJECT_OTHER=''
Index: suite.rc
===================================================================
--- suite.rc    (revision 105592)
+++ suite.rc    (working copy)
@@ -27,6 +27,9 @@
     [[events]]
         mail events = shutdown

########################################## This is where the cylc parameterisation occurs

+    [[parameters]]
+        ensemble = {{ range(ENSEMBLE_SIZE) | join(', ') }}
+
 [scheduling]
 
################################### In scheduling, all the atmos_main are replaced by atmos_main<ensemble> 
 
     cycling mode        = {{CALENDAR}}
@@ -62,7 +65,7 @@
 
         {% if RUN %}
         [[[ R1 ]]]
-            graph = recon => atmos_main
+            graph = recon => atmos_main<ensemble>
         {% endif %}
 {% endif %}
 
@@ -69,20 +72,20 @@
 {% if RUN %}
         {% if BUILD_UM %}
         [[[ R1 ]]]
-            graph = {{FCMUM_LAST}} => atmos_main
+            graph = {{FCMUM_LAST}} => atmos_main<ensemble>
         {% endif %}
 
         [[[ R1 ]]]
-            graph = install_ancil => atmos_main
+            graph = install_ancil => atmos_main<ensemble>
 
         [[[ {{RESUB}} ]]]
-            graph = atmos_main => housekeeping
+            graph = atmos_main<ensemble> => housekeeping
 
         {% if POSTPROC %}
         [[[R1]]]
             graph = fcm_make_pp => fcm_make2_pp => postproc
         [[[ {{RESUB}} ]]]
-            graph = atmos_main => postproc => housekeeping
+            graph = atmos_main<ensemble> => postproc => housekeeping
         {% if PPTRANSFER %}
         [[[R1]]]
             graph = fcm_make_pptransfer => fcm_make2_pptransfer => pptransfer
@@ -106,7 +109,7 @@
 
         {% if ARCH_WALL %}
         [[[ R1//^+{{RUNLEN}}-{{RESUB}} ]]]
-            graph = atmos_main => rose_arch_wallclock
+            graph = atmos_main<ensemble> => rose_arch_wallclock
         {% endif %}
 
 # Include tests graph if required
@@ -191,11 +194,23 @@
 
     [[recon]]
         inherit = RUN_MAIN, RCF_RESOURCE, RECONFIGURE
+        [[[environment]]]
+            ASTART=$ROSE_DATA/$RUNID.astart
+            DATAM=$ROSE_DATA/{{DATAM}}
+            ENS_MEMBER=0
 
     [[atmos_main]]
         inherit = RUN_MAIN, ATMOS_RESOURCE, ATMOS
         post-script = save_wallclock.sh {{RESUB}}

####################################### Creating the atmos_maion<ensemble> task
####################################### The ASTART and DATAM files _hopefully_ contain correct file paths

+    [[atmos_main<ensemble>]]
+        inherit = RUN_MAIN, ATMOS_RESOURCE, ATMOS
+        post-script = save_wallclock.sh {{ RESUB }}
+        [[[environment]]]
+            ASTART=${ROSE_DATA}/$RUNID.astart_${CYLC_TASK_PARAM_ensemble}
+            DATAM=$ROSE_DATA/{{DATAM}}/ensemble_$CYLC_TASK_PARAM_ensemble
+            ENS_MEMBER=$CYLC_TASK_PARAM_ensemble
+            
     [[fcm_make_pp]]
         inherit = RUN_MAIN, EXTRACT_RESOURCE
     [[fcm_make2_pp]]
```

So, to summarise, the changes to the `suite.rc` file consist of:
* `ensemble = {{ range(ENSEMBLE_SIZE) | join(', ') }}` which sets up the list of ensemble members inside a `[[parameter]]` block within the `[cylc]` heading

* changing all the dependencies graphs to include a new object `atmos_main<ensemble>` where `<ensemble>` inserts the ensemble values in different instances of the `atmos_main` task

* adding in an `atmos_main<ensemble>` task which overwrites the data output location. **THIS HAS NOT BEEN CORRECTLY IMPLEMENTED ABOVE**

In addition to this, 

* the `rose-suite.conf` file needs the `ENSEMBLE_SIZE` parameter to be present to pass to suite.rc on running
* the `app/um/rose-app.conf` file needs `ENS_MEMBER=${ENS_MEMBER}` which is passed from suite.rc and replaces default value of `ENS_MEMBER=0`
* the `meta/rose-meta.conf` file should have an entry for the Ensemble_Size parameter so that it appears in the rose GUI. This is not however necessary as if this is not included the parameter will apper in the `jinja2` heading in the GUI. 

Work is still ongoing to ensure that :
1. the output data goes to the correct place
2. a perturbation can be done

With regards to the perturbation, the UKESM suite with Ensembles uses a python script as an app which does perturbations. This structure can be reused so that the graph would look like:

```
                     fcm_make_um
                          |
                    fcm_make2_um --¬  install_ancils
                     |    \   \     \/   /      /| 
                     |     \   \   / \ /      /  |
                     |      \   \/   /\     /    | 
                     |       \ / recon \  /      |
                     |       /\/  |   \ X        |
                     |     / / \  |   /\ \       |
                     |   / /    \ | /   \ \      |
                  perturb0    perturb1    perturb2
                     |           |           | 
                     |           |           | 
                  atmos_       atmos_      atmos_
                  main_        main_       main_
                  ensemble0    ensemble1   ensemble2
                        \         |         /
                         \        |       /
                          \       |     /      
                            housekeeping
```


## Current state as of 2019/02/21

Perturbed ensembles can now be generated thanks to a combination of edits to the suite.rc file and the addition of a python script as an app. The important bits to allow this are :
1. A parameterised cylc task using a `[[parameter]]` block in the runtime section, which then allows ensembles of tasks with the identifier `taskname<ens>` for an ensembling parameter `ens`
2. Use of the `ROSE_APP_OPT_CONF_KEYS` variable to select namelist additions which can overwrite previous entries in the namelist.

Descriptions of how this is done can be found in the folder "Ensemble Generation"
