# Add a perturbation

By including the section `=> perturb<ensemble> => atmos_main<ensemble>` in place of any instance of the atmos_main call in the cylc graph, you can call a new rose task called "_perturb_".

Work has currently been done integrating this into both the UKCA release job 5 suite copy (u-bf737) and a new suite released yesterday which is set up for the UKESM 1.0 (u-bf910), which will be the program being used for the actual PPe system, and also setting up ensembles on this.

One plan is to change the value of the ACURE variables in the <<SHARED>> file (path `/work/19880901T0000Z/atmos_main_ensemble0/SHARED`) using the perturb script. The problem with this is that the SHARED file does not exist on archer until the first `atmos_main` job is queued. The values instead are copied over from the `app/um/rose-app.conf` file. 

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
This means that individual rose-app.conf files will be needed for the different ensembles, and the ensembles will have to somehow know to read these files. I honestly have no idea how this can be done yet. A different possibility would be for the perturb script to change the rose-app.conf file and for there to be an exclusive run condition for the `perturb=>atmos_main` groups 
