# Regional Perturbations #

A set of regional perturbations has been requested for the A-CURE project. This can be achieved by applying a regional mask and perturbing dependant on the region. To do this a set of ancillaries, one for each of the regionally perturbed variables, needs to be created, and then read into the UM.

## Creating the ancillaries ##

The ancillaries were made based on the N96 resolution land-sea mask. The mask was taken from the UM (at location `/work/y07/y07/umshared/ancil/atmos/n96e/orca1/land_sea_mask/etop01/v2/qrparm.mask` on Archer), and converted to NetCDF using xconv. Following this, the map grid was extracted from the netCDF file and reconfigured to be displayed as a 192 x 144 grid that could be input into a csv file using the [python script](./gridbuild.py) included in this folder. These grids were manually modified to select the regions with integer values denoting each region as described below:

Anth SO2:
1. North America
2. Europe
3. China
4. Rest of Asia
5. Rest of Globe

FF Emiss:
1. Europe
2. North America
3. China
4. Rest of Asia
5. Marine
6. Rest of Globe

Res Emiss:
1. China
2. Rest of Asia
3. Africa
4. Latin America
5. Rest of Globe

BB Emiss:
1. Southern Hemisphere South America (ARCD & SARC)
2. Northern Hemisphere Africa (NHAF)
3. Southern Hemisphere Africa (SHAF)
4. Boreal Norhtern Hemisphere (BONA-W, BONA-E, BOAS)
5. Rest of Northern Hemisphere (TENA, CEAM, NHSA, EURO, CEAS, SEAS, MIDE, Northern Hemisphere Marine)
6. Rest of Southern Hemisphere (EQAS, AUST, Southern Hemisphere Marine)

The definitions of the regions for the first three files correspond to the mapping of regions given in table S2 of doi [10.5194/gmd-11-369-2018](http://dx.doi.org/10.5194/gmd-11-369-2018), and the regions for BB Emiss correspond to the regions shown in figure 2 of doi [10.5194/gmd-10-3329-2017](http://dx.doi.org/10.5194/gmd-10-3329-2017). The csv files with the modified regions are included in this folder.

After setting the regions, the grids were refactored again into NetCDF-compatible formats (included in this folder as the * .cdl files) and then converted back into netCDF using `ncgen`. Finally, the netCDF files were converted into UM ancil files using XANCIL on Archer.

## Inputting the new ancils into the UM ##

Input of ancillaries is achieved through STASH. I was able to incorporate the ancillaries using the instructions [here](https://code.metoffice.gov.uk/trac/um/wiki/NewProgAncil) and [here](https://code.metoffice.gov.uk/doc/um/latest/um-training/stashmaster.html#method2)

The first step was to add a STASH record in the STASHmaster_A file. The STASH codes 301-320 have been reserved for single level user ancillaries with real values in the main field, so the codes 301-304 were used for the ancillaries (note: these STASH codes need to be used in the xancil conversion also). By using the 301-304 STASH codes, no extra code changes were needed to get the ancil values into the UM.

The STASHmaster_A records have the following form (read UMDP C04 for details):

```
1|    1 |    0 |  301 | ACURE ANTH SO2 REGIONAL MASK       |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    2 |    0 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |  395 |    0 |  129 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  302 | ACURE BB EMISS REGIONAL MASK       |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    2 |    0 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |  395 |    0 |  129 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  303 | ACURE FF EMISS REGIONAL MASK       |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    2 |    0 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |  395 |    0 |  129 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  304 | ACURE RES EMISS REGIONAL MASK      |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    2 |    0 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |  395 |    0 |  129 |    0 |    0 |    0 |    0 |    0 |
```

Corresponding entries are required in the STASHmaster-meta.conf file:

```
[stashmaster:code(301)]
description=ACURE ANTH SO2 REGIONAL MASK
help= Perturbation regions for Anthropogenic SO2 used in A-CURE project
    =
    = Regions
    = 1: North America
    = 2: Europe
    = 3: China
    = 4: Rest of Asia
    = 5: Rest of Globe

[stashmaster:code(302)]
description=ACURE BB EMISS REGIONAL MASK
help= Perturbation regions for Biomass Burning Emissions used in A-CURE project
    = Modelled after regions used for GFED4 emission inventory
    =
    = Regions
    = 1: Southern Hemisphere South America (ARCD, SARC)
    = 2: Northern Hemisphere Africa (NHAF)
    = 3: Southern Hemisphere Africa (SHAF)
    = 4: Boreal Northern Hemisphere (BONA, BOAS)
    = 5: Rest of Northern Hemisphere (TENA, CEAM, NHSA, EURO, CEAS, SEAS, MIDE, NH-Marine)
    = 5: Rest of Southern Hemisphere (EQAS, AUST, SH-Marine)

[stashmaster:code(303)]
description=ACURE FF EMISS REGIONAL MASK
help= Perturbation regions for Fossil Fuel Emissions used in A-CURE project
    =
    = Regions
    = 1: Europe
    = 2: North America
    = 3: China
    = 4: Rest of Asia
    = 5: Marine
    = 5: Rest of Globe

[stashmaster:code(304)]
description=ACURE RES EMISS REGIONAL MASK
help= Perturbation regions for Residential Emissions used in A-CURE project
    =
    = Regions
    = 1: China
    = 2: Rest of Asia
    = 3: Africa
    = 4: Latin America
    = 5: Rest of Globe
```

The ancils are added to the STASH list and the ITEMS namelist in the `um/namelist/Reconfiguration and Ancillary Control/Configure ancils and initialise dump fields` panel in rose edit.

To ensure that the STASHmaster_A file is picked up, it can be copied over to the rose suite app/um/file folder, ensuring that it gets copied over to Archer from puma. This should be done alongside the setting of the "STASHMASTER" environment variable to "." in `um/env/Runtime Controls`

Following these steps should ensure that the ancil data is picked up and incorporated into the main D1 array within the UM. To access these within the UKCA, follow the instructions in UMDP 304.
