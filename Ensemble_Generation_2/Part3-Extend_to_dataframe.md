#Extending the Perturbation System to use a Dataframe

The python script used to generate the ensemble configuration files has been expanded significantly. This has involved changing the way in which the new values for variables are determined - replacing a calculated value with values read in from a csv file which is converted into a pandas dataframe. The script now takes 4 arguments:
1. The input template location, which contains logical switches for all of the PPE variables
2. The Output location for the optional configuration files, typically `app/um/opt`
3. The dataframe location
4. A space separated list of the ensemble indices.

To be able to pass the correct values to the python script and control the ensembling, a set of modifications were needed for the suite.rc, rose-suite.conf and meta/rose-meta.conf files to allow for new parameters.

The changes in the rose-suite.conf file allow the introduction of new rose variables:
```
ENS_CSV='app/perturb/bin/ppe_dataframe.csv' # Location of the CSV dataframe file
ENS_L='1','4','37-41'  # a list of individual ensemble members to run. Ranges enabled with hyphens
ENS_N=5   # Number of ensemble members, used when autopopulating a range
ENS_RNGE=true  # Logical for whether ensembles are itemised or ranged (f/t respectively)
ENS_STRT=3    ~ Starting ensemble member for range
ENS_TMPL='app/perturb/bin/template.conf'  # Location of PPE template file
```
In the meta/rose-meta.conf file, these variables are populated in a new pane under the "suite conf" heading in the rose GUI,

```ini
[jinja2:suite.rc=ENS_RNGE]
ns          = Ensembles
compulsory  = false
help        =
title       = Run a range of members?
type        = boolean
trigger     = jinja2:suite.rc=ENS_STRT: true;
            = jinja2:suite.rc=ENS_N: true;
            = jinja2:suite.rc=ENS_L: false
sort-key    = 5a

 [jinja2:suite.rc=ENS_STRT]
ns          = Ensembles
compulsory  = false
help        =
title       = First ensemble member in range
type        = integer
sort-key    = 5b

[jinja2:suite.rc=ENS_N]
ns          = Ensembles
compulsory  = false
help        =
title       = Number of ensemble members
type        = integer
sort-key    = 5c

[jinja2:suite.rc=ENS_L]
ns          = Ensembles
compulsory  = false
title       = Individual ensemble members
help        = A range of members may be described using a hyphen.
type        = character
length      = :
sort-key    = 5d

[jinja2:suite.rc=ENS_CSV]
ns          = Ensembles
compulsory  = true
title       = Ensemble csv file
help        = A csv file that defines parameter perturbations.
type        = character
length      = 9999
sort-key    = 4

[jinja2:suite.rc=ENS_TMPL]
ns          = Ensembles
compulsory  = true
title       = Ensemble template namelist file
help        = A configuration file with a list of logical
            = switches for each of the parameters in the PPE.
type        = character
length      = 9999
sort-key    = 3
```

Finally, the suite.rc file is modified to allow non-sequential ensemble numbers as below:
```ini
{% set ENS_COMB = [] %}
{% for mem in (range(ENS_STRT,ENS_STRT+ENS_N)|list if ENS_RNGE else ENS_L) %}
    {% if '-' in mem|string %}
        {% for msub in range(mem.split('-')[0]|int, mem.split('-')[1]|int+1) %}
            {% do ENS_COMB.append(msub) %}
        {% endfor %}
    {% elif mem != '' %}
        {% do ENS_COMB.append(mem) %}
    {% endif %}
{% endfor %}



[cylc]
    UTC mode = True
    [[events]]
        mail events = shutdown

    [[parameters]]
        ens = {{ ENS_COMB | join(', ') }}
```
and the perturb task is also modified such that it now reads
```ini
    [[perturb]]
        inherit = PERTURB_RESOURCE
        script = "rose task-run --app-key=perturb --verbose"
        [[[environment]]]
            ENS_LST={{ ENS_COMB | join(' ') }}
            TMPL_LOC=${CYLC_SUITE_RUN_DIR}/{{ENS_TMPL}}
            CONF_LOC=${CYLC_SUITE_RUN_DIR}/app/um/opt
            DF_LOC=${CYLC_SUITE_RUN_DIR}/{{ENS_CSV}}
```

The important difference between this and the previous version of this task is that the perturb task is only called once at the beginning, rather than once per ensemble member as it was previously. As such it now takes a space separated list as the ENS_LST variable, which is passed into the python script.
