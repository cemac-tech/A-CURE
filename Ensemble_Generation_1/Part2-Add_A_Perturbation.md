
# Add a perturbation

By including the section `=> perturb<ensemble> => atmos_main<ensemble>` in place of any instance of the atmos_main call in the cylc graph, you can call a new rose task called "_perturb_".

Work has currently been done integrating this into both the UKCA release job 5 suite copy (u-bf737) and a new suite released yesterday which is set up for the UKESM 1.0 (u-bf910), which will be the program being used for the actual PPe system, and also setting up ensembles on this.

One plan is to change the value of the ACURE variables in the `<<SHARED>>` file (path `/work/19880901T0000Z/atmos_main_ensemble0/SHARED`) using the perturb script. The problem with this is that the SHARED file does not exist on archer until the first `atmos_main` job is queued. The values instead are copied over from the `app/um/rose-app.conf` file. 

A possible solution would be to create a copy of rose-app.conf (possibly rose-app.conf_ensemble0) which is then read in by the ensemble member. This is almost vital, as according to a recent test, the run conditions for consecutive cycles are taken from the rose-app.conf file - **NOT** the SHARED file from the previous cycle. This was seen when changing the value in the rose-app.conf file of the A-CURE test variable and seeing the following diff
```diff
< run_target_end=0,0,1,0,0,0,
---
> run_target_end=0,0,2,0,0,0,
177c177
< ukca_aeros_volc_so2=21.7,
---
> ukca_aeros_volc_so2=29.1,
```
This means that individual rose-app.conf files will be needed for the different ensembles, and the ensembles will have to somehow know to read these files. 

_______________________________________________________________________

After further investigation, it seems that changing the um rose-app.conf file is the way forward. In-built settings allow this through putting a specific configuration file into the folder `/app/um/opt` which is then appended onto the end of the master rose-app.conf file. This can then override matching values in the main configuration file (you can find details about optional configurations [here](https://metomi.github.io/rose/doc/html/api/configuration/rose-configuration-format.html#optional-configuration) and [here](https://metomi.github.io/rose/doc/html/tutorial/rose/furthertopics/optional-configurations.html))

A python app is called in `u-bf737` which takes as arguments the location of a file with the namelist variable in it, changes it, and then places the resulting file in the `app/um/opt` folder. This python script is called by cylc.

Changes to the various configuration files are given below which are needed to implement the ensembles properly:

|File   |Original                       |With Ensemble|
--------|-------------------------------|-------------------------------------------------------
|./suite.rc|<pre> [cylc]<br>    UTC mode = True<br>    [[events]]<br>        mail events = shutdown</pre>|<pre> [cylc]<br>    UTC mode = True<br>    [[events]]<br>        mail events = shutdown<br>    [[parameters]]<br>        ens = {{ range(ENSEMBLE_SIZE) \| join(', ') }}</pre>|
|       |`graph = recon => atmos_main`    |`graph = recon => perturb<ens> => atmos<ens>`|
|       |`graph = {{FCMUM_LAST}} => atmos_main` | `graph = {{FCMUM_LAST}} => perturb<ens> => atmos<ens>`|
|       |`graph = install_ancil => atmos_main`| `graph = install_ancil => perturb<ens> => atmos<ens>`|
|       |`graph = atmos_main => housekeeping` | `graph = atmos<ens>[-{{RESUB}}] => atmos<ens> => housekeeping`|
|       |`graph = atmos_main => postproc => housekeeping` | `graph = atmos<ens>[-{{RESUB}}] => atmos<ens> => postproc => housekeeping`|
|       |`graph = atmos_main => rose_arch_wallclock` | `graph = atmos<ens> => rose_arch_wallclock`|
|       |`ROSE_APP_OPT_CONF_KEYS = ({{HORIZ}}) ({{CALENDAR}}) {{BITCOMP_NRUN_OPT}} {{UM_OPT_KEYS}}`| This line is commented out as the optional configuration keys are included later on in the `atmos<ens>` block|
|       |<pre>[[RUN_MAIN]]<br>    [[[environment]]]<br>        DATAM = =$ROSE_DATA/{{DATAM}}</pre>| The environment block is commented out as `DATAM` is defined later on|
|       |<pre>[[recon]]<br>    inherit = RUN_MAIN, RCF_RESOURCE, RECONFIGURE<br>    [[[environment]]]<br>            DATAM=$ROSE_DATA/{{DATAM}}</pre>|<pre>[[recon]]<br>    inherit = RUN_MAIN, RCF_RESOURCE, RECONFIGURE<br>    [[[environment]]]<br>        ASTART=$ROSE_DATA/$RUNID.astart<br>        DATAM=$ROSE_DATA/{{DATAM}}<br>        ENS_MEMBER=0</pre>|
|       |<pre>[[atmos_main]]<br>    inherit = RUN_MAIN, ATMOS_RESOURCE, ATMOS<br>    post-script = save_wallclock.sh {{RESUB}}</pre>|<pre>[[atmos\<ens>]]<br>    inherit = RUN_MAIN, ATMOS_RESOURCE, ATMOS<br>    post-script = save_wallclock.sh {{ RESUB }}<br>    [[[environment]]]<br>        ASTART=${ROSE_DATA}/$RUNID.astart<br>        DATAM=$ROSE_DATA/{{DATAM}}/ens_${CYLC_TASK_PARAM_ens}<br>        ENS_MEMBER=${CYLC_TASK_PARAM_ens}<br>        ROSE_APP_OPT_CONF_KEYS = ens${ENS_MEMBER} ({{HORIZ}}) ({{CALENDAR}}) {{BITCOMP_NRUN_OPT}} {{UM_OPT_KEYS}}</pre>|
|       | There is no block here as the perturb block is a wholly new addition|<pre>[[perturb<ens>]]<br>    inherit = PERTURB_RESOURCE<br>    script = "rose task-run --app-key=perturb --verbose"<br>    [[[environment]]]<br>        TMPL_LOC=${CYLC_SUITE_RUN_DIR}/app/perturb/bin/ukca_acure.conf<br>        CONF_LOC=${CYLC_SUITE_RUN_DIR}/app/um<br>        LABEL=${CYLC_TASK_PARAM_ens}</pre>|
|./meta/rose-meta.conf| No ensemble related entry |<pre> [jinja2:suite.rc=ENSEMBLE_SIZE]<br> compulsory=true<br> description=Size of the A-CURE ensemble<br> help=This is the ACURE ensemble variable<br>     =<br>     =This variable is currently only present for testing the rose suite,<br>    =and does not have any impact on the values of parameters<br> ns=cycle<br> range=0:50<br> sort-key=999<br> title=A-CURE Ensemble Size<br> type=integer</pre>|
|./rose-suite.conf| No ensemble variable exists|`ENSEMBLE_SIZE=5`|
|./site/archer.rc | No information related to perturb included |<pre>[[PERTURB_RESOURCE]]<br>    inherit = HPC_SERIAL<br>    pre-script = """<br>                 module load anaconda<br>                 export PYTHONPATH=$PYTHONPATH:$UMDIR/lib/python2.7<br>                 """<br>    [[[directives]]]<br>        -l walltime=00:05:00<br>|

In addition to these changes, a new app folder is created in `./app/` called `perturb`. This folder contains a file `rose-app.conf` and a two further files in the `bin` folder, named `perturb_ini` and `ukca_acure.conf`. All three of these files are shown in this folder, with the appropriate folder structure preserved.

