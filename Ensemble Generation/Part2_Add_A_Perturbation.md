# Add a perturbation

By including the section `=> perturb<ensemble> => atmos_main<ensemble>` in place of any instance of the atmos_main call in the cylc graph, you can call a new rose task called "_perturb_".

Work has currently been done integrating this into both the UKCA release job 5 suite copy (u-bf737) and a new suite released yesterday which is set up for the UKESM 1.0 (u-bf910), which will be the program being used for the actual PPe system, and also setting up ensembles on this.

One plan is to change the value of the ACURE variables in the <<SHARED>> file (path `/work/19880901T0000Z/atmos_main_ensemble0/SHARED`) using the perturb script. The problem with this is that the SHARED file does not exist on archer until the first `atmos_main` job is queued. The values instead are copied over from the `app/um/rose-app.conf` file. 

A possible solution would be to create a copy of rose-app.conf (possibly rose-app.conf_ensemble0) which is then read in by the ensemble member
