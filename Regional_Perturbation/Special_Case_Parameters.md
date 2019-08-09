 # Special Case Parameters and Upgrades to the [pertub_ppe.py](Regional_Perturbations/perturb/bin/perturb_ppe.py) script #

 In this round of revisions, the perturb script was largely revised, adding comments and cleaning up the code for readability, while also adding code to deal with special cases and parameters that are treated differently.

 ## Recap and Overview of function of the perturb_ppe.py script and related files ##

The perturbation script takes as its input a set of four arguments:
* The location of a template configuration file, exported with the name TMPL_LOC by Cylc. The template file is to contain every parameter to be perturbed as a logical set to false with the exception of the regional parameters, as each parameter has an associated logical in the um namelists.
* The location of the optional configuration files folder, typically `app/um/opt` in the Cylc run directory, exported with the name CONF_LOC by Cylc. This directory is where the output files from the script will be placed, for later reading by the atmos_main task in the UM run.
* The ensemble number list, referred to as ENS_LST by Cylc. This takes the form of a space separated list of ensemble members which can be itterated through to create the required specific optional configuration files.
* The location of the PPE dataframe. This dataframe is used by the script to provide the
PPE variables their values and logical switches. A discussion of the required layout of the dataframe file is included later.

In its function, the perturb script iterates through the elements of the ensemble list. Each ensemble member corresponds to a row in the dataframe. For each ensemble member, the script reads the relevant row, and for variables that have a value other than "x" (used to indicate no perturbation) the replacement value is written to the output configuration file on the line following the corresponding logical switch, ie the section of the template file which reads:
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
l_acure_dry_dep_acc=.false.
acure_dry_dep_acc=0.26
l_acure_dry_dep_so2=.false.
acure_dry_dep_so2=17.54
l_acure_kappa_oc=.false.
l_acure_sig_w=.false.
acure_sig_w=0.001
```

The parameter values written to this output file overwrite the corresponding parameters in the `app/um/rose-app.conf` file from which the UM obtains values input in the GUI - as such if ensembling is enabled the **values from the dataframe completely supercede those from the GUI**.

## Input File Formats ##

There are two input files for the perturbation script - the *template.conf* file and the *ppe_dataframe* file. Both of these files have format rules described below.

### Template File ###

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

### PPE dataframe file

* The PPE dataframe file should be a csv file
* There should be one row per ensemble member
* If a parameter is not required as part of an ensemble members parameter set, the value in the dataframe should be "x", otherwise a real number is expected.
* There should be the same number of columns in the dataframe as there are logical switches in the template file, plus an additional 22 rows for the four regional parameters.
* Each column should have a header field with the name of the *real* parameter being perturbed. These parameters should match the name of the associated logical, omitting the `l_` characters, i.e. `l_acure_bl_nuc` in the template file would have a corresponding column in the dataframe with the header `acure_bl_nuc`. The columns need not be in the same order as in the template file however.
* The regional perturbation logical switches have their own corresponding real value in the dataframe, but *only* in the dataframe - there is no variable called `acure_anth_so2` or `acure_carb_ff_ems` in the code itself. These columns are used only as an indicator to direct the script to add the desired regional fields.
* **There is no verification to ensure that all the regional perturbation fields have values when one has values.** This means that there is no check that prevents  an optional configuration file from reading:
```
            :
            :
l_acure_carb_ff_ems=.false.
l_acure_carb_bb_ems=.true.
acure_carb_bb_ems_sam=0.65
acure_carb_bb_ems_naf=0.004
acure_carb_bb_ems_saf=x
acure_carb_bb_ems_bnh=x
acure_carb_bb_ems_rnh=0.01
acure_carb_bb_ems_rsh=2.88
l_acure_carb_res_ems=.false.
l_acure_anth_so2=.false.
            :
            :
```
which would cause a fatal error. As such **it is the responsibility of the dataframe creator** to ensure that the regional parameters are enabled as a block, and that if they are enabled then there is a value for the non-regional version of that parameter on that row.

##  Special Cases ##

There are some special cases that have been added to the perturbation script. These consist of
1. Regional Perturbations
2. Parameters with a new logical but no associated variable
3. Perturbation of parameters that already exist in the model
4. Perturbation parameters which affect more than one preexisting variable

Below the code to deal with these special cases, along with the standard case where new acure parameters are inserted in *logical-real* pairs. The code is extensively commented to give insight. Before this code block, the variables to be perturbed and their associated values are written into a dictionary `ens_dict` where the keys are the variable names, and `num_replaced` is the number of logical switches in the template file whose values were changed from `.false.` to `.true.`. The template file with the logical switches changed is stored in a variable called `outlst` which at the end is written out to form the entire optional configuration file.

```python
    if num_replaced == len(ens_dict):

        # One-to-one relationship between logicals and dataframe columns, meaning
        # compatibility between template conf file and ensemble parameter set.
        # No regional perturbations, as these have 5 or 6 parameters controlled
        # by a single logical switch.

        for (name,val) in ens_dict.items():
            teststr="_"+name+"="               # name is the variable name
            for i, line in enumerate(outlst):
                if teststr in line:
                    if (name=="acure_pcalc_index"):
                        # Special case where a pcalc file is selected by modifying the string for ukcaprec to include the pcalc_index (special case 2 above)
                        outlst.insert(i+1,"ukcaperc=/work/n02/n02/lre/pcalc_UKESM_11_1_default/RADAER_pcalc_"+val+".ukca\n")
                    else:
                        outlst.insert(i+1,name+"="+val+"\n")
                        #adds new value for parameter to output conf file below the
                        #corresponding logical switch
                        if (name=="m_ci" or name=="a_ent_1_rp"):
                            # Special cases where min and max values are set to the same as the parameter (special case 4 above)
                            outlst.insert(i+2,name+"_min="+val+"\n")
                            outlst.insert(i+3,name+"_max="+val+"\n")
        # Write out the output opt conf file
        with open(dest, 'w') as fout:
            for line in outlst:
                fout.write("%s" % line)

    else:

        # Regional perturbations may exist, throwing off the compatibility between
        # template conf file and ensemble parameter set.

        if (any("_ems_" in key for key in ens_dict) or
           any("_so2_" in key for key in ens_dict)):
            # used as indicator for region specific parameters
            noreg = [key for key,value in ens_dict.items() if ((not "_so2_" in key) and (not "_ems_" in key))]
            # Builds up a list of variables which are not region specific
            noreg_dict = {}
            for k in noreg: noreg_dict.update({k:ens_dict[k]})
            # Set a dict of parameters without the regional parameter values
            # Allowing compatibility check between template conf file and
            # the ensemble member parameter set
            if num_replaced == len(noreg_dict):
                for (name,val) in noreg_dict.items():
                    teststr="_"+name+"="
                    for i, line in enumerate(outlst):
                        if teststr in line:
                            # Proceed as before if not one of the regional perturbation
                            # special case variables
                            if (name!="acure_carb_ff_ems" and
                                name!="acure_carb_bb_ems" and
                                name!="acure_carb_res_ems" and
                                name!="acure_anth_so2"):
                                if (name=="acure_pcalc_index"):
                                    outlst.insert(i+1,"ukcaperc=/work/n02/n02/lre/pcalc_UKESM_11_1_default/RADAER_pcalc_"+val+".ukca\n")
                                else:
                                    outlst.insert(i+1,name+"="+val+"\n")
                                    if (name=="m_ci" or name=="a_ent_1_rp"):
                                        outlst.insert(i+2,name+"_min="+val+"\n")
                                        outlst.insert(i+3,name+"_max="+val+"\n")

                            else:
                                # Add the regional perturbation values rather than the indicator
                                # variable value
                                temp_dict={}
                                j=1
                                regkeys=[key for key,value in ens_dict.items() if name+"_" in key]
                                #List of region specific variables created
                                for k in regkeys: temp_dict.update({k:ens_dict[k]})
                                for (name2,val2) in temp_dict.items():
                                    outlst.insert(i+j,name2+"="+val2+"\n")
                                    j +=1

                with open(dest, 'w') as fout:
                    for line in outlst:
                        fout.write("%s" % line)
            else:
                # Handle a mismatch between the template conf file and ensemble member
                # parameter list
                for (name,val) in ens_dict.items():
                    print("%s=%s\n" % (name,val))
                raise FileError('One or more of the expected parameters were not replaced.\n'
                                + 'Expected to replace '+str(len(ens_dict))+' but replaced '
                                + str(num_replaced)+" in " + dest + "\n")
        else:
            # Handle a mismatch between the template conf file and ensemble member
            # parameter list
            for (name,val) in ens_dict.items():
                print("%s=%s\n" % (name,val))
            raise FileError('One or more of the expected parameters were not replaced.\n'
                            + 'Expected to replace '+str(len(ens_dict))+' but replaced '
                            + str(num_replaced)+" in " + dest + "\n")
```

The above code block handles special cases 1, 2, and 4, however special case 3 is handled differently by allowing variables from the dataframe without the word "acure" in the name to be paired to logicals that do, such as `bparam <--> l_acure_bparam`. An example of this being done is below, where the dictionary `rep` is created and then compiled into a regular expression to change the values of the logicals
```python
    for name in ens_dict:
        if re.search("_ems_",name) is not None:
            #Skip the regional variables as no one-to-one associated logical
            continue
        elif re.search("_so2_",name) is not None:
            #Skip the regional variables as no one-to-one associated logical
            continue
        elif re.search("acure",name) is not None:
            # Standard acure only variables
            rep["l_"+name+"=.false."]="l_"+name+"=.true."
        else:
            # Pre-existing variables
            rep["l_acure_"+name+"=.false."]="l_acure_"+name+"=.true."

    pattern = re.compile("|".join(rep.keys()))
```  
