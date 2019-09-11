## [Previous ACURE Work](https://github.com/cemac/A-CURE) ##


**CEMAC:**
James developed 2 python scripts to extract and condense information from PPE
pp file to nc files aresol optical depth: (AOD) 550nm and Cloud Dropled number
Concentration (CDNC). The second script collates AOD 3 hourly into monthly averages.

[collocate_UM_AERONET.py](https://github.com/cemac/A-CURE/blob/master/Python/collocate_UM_AERONET.py) provides an example of using CIS tools (may be useful to look at for CRESCENDO).

**Version 8.4**

In the older version he and Masaru developed scripts that would generate
additional experiments for the ensemble members after the UMUI had
processed and created the first experiment job.

I think the perturbations are in the number of order (200+) and
effectively tweak the namelist entries for particular aerosol parameters.

<hr>

## NCAS CMS INPUT ##

Mark spoke to Grenville and Dave Sexton from NCAS and Dave has some ensemble runs. Using Rosie Go
these can easily be searched for: run u-bd149. (CS copy has been made at u-bf163)

>""This is a UKESM suite, but if you look at this, you
can hopefully see the logic behind the cylc graph. For example, in the
suite.rc there are parameters:

>     [[parameters]]
         ensemble = {{ range(ENSEMBLE_SIZE) | join(', ') }}

>so <ensemble> can be used as a variable when handling the work. You can
see lines such as:

>'recon => perturb<ensemble> => coupled<ensemble>'

>in the graph, so there will be an app called perturb which will tweak
the start files in some way, and the coupled<ensemble> will run that job
in the appropriate directory. What I imagine that you will be able to do
would be to keep the logic from this suite whilst running your own
model, but choose some appropriate script to run at the perturb stage.
Currently, if you look at app/perturb/rose-app.conf , you will see that:
[command]
>default=module load anaconda && perturb_ini.py ${ASTART} -o
${ASTART}_${LABEL}

>but you will have to replace the perturb_ini.py script, and set this up
to copy files, or run your own script as you feel appropriate. LABEL is
set in the suite.rc as LABEL= $CYLC_TASK_PARAM_ensemble

>This is just an example of something which is similar to what I think
that you are doing, and I don't know if it's the easiest or best way,
but if you have a look at it I hope that it'll help. Let me know if it
does(n't),""

<hr>

## Meeting 20th Dec ##

### Overview ###

* 200 member ensemble, 50 parameters being edited.
* Control runs are a median member
* 2 suties for 2 different emmision scenarios
* 3.5 months CEMAC time

### Configuration information ###

* UKESM atmosphere only release job
 * coming soon to ARCHER
 * **beware** known bugs
* Large amount of data produced, needs intermediate house keeping.
* No chemistry
* Extreme perturbations are likely to crash
  * must allow the ensemble to continue
* Ultimate destination of data will be to JASMIN

### Other Useful information ###

* Leighton has provided previous scripts in this folder
* Dave Sexton (MetOffice) previously has done hand edits
* UKESM and UKCA release job do not work the same scientifically - be aware of.

### Actions ###

* Once UKESM on Archer can begin the process
* Leighton has estimated 25 days for this part.
