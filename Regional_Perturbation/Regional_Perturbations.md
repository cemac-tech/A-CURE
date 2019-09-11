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
### Part 1 - Addition to D1 array ###

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
    = 6: Rest of Southern Hemisphere (EQAS, AUST, SH-Marine)

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
    = 6: Rest of Globe

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

Following these steps should ensure that the ancil data is picked up and incorporated into the main D1 array within the UM. To access these within the UKCA, follow the instructions in UMDP 304, summarised in the following section in specific for the ACURE PPE regional ancillaries.

### Part 2 - Interfacing with the UKCA ###

The UKCA portion of the UM accesses data from the D1 array through a subset array called UkcaD1codes. The construction and population of this array is achieved across two different files:

* `ukca_um_interf.F90` is the file containing subroutines which interface with the UM D1 array into which the ancillaries were included in the previous step.
* `ukca_setd1defs.F90` allocates and populates the UkcaD1codes array

In the UMDP 304 documentation, the arrays are declared in the `ukca_main1_ukca_main1.F90` file, but this is no longer the case - this is now done in `ukca_um_interf.F90`.

#### Stage 1 - Declare the new arrays

To declare the arrays that the UkcaD1codes array will use you need to add the following lines to the preamble section of `ukca_um_interf.F90`

```
REAL, ALLOCATABLE, PUBLIC :: acure_anth_so2_region (:,:)
REAL, ALLOCATABLE, PUBLIC :: acure_ff_emiss_region (:,:)
REAL, ALLOCATABLE, PUBLIC :: acure_res_emiss_region (:,:)
REAL, ALLOCATABLE, PUBLIC :: acure_bb_emiss_region (:,:)
```
These variable declarations should be in the order they are allocated later, ie in stash code order, so these lines being related to STASH codes 301-304 should be added on line 106 between the declarations for `cloud_liq_frac(:,:,:)` and `biogenic(:,:,:)`

#### Stage 2 - Allocate the arrays and match to a STASH code

In the same file but further down, in the `getd1flds` subroutine, you need to add code to allocate the arrays and link them to a STASH code. This is done by adding the following lines into the correct place in the prognostic block (ie again between `cloud_liq_frac` and `biogenic`):

```
  CASE (301)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(acure_anth_so2_region(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,acure_anth_so2_region)
  CASE (302)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(acure_bb_emiss_region(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,acure_bb_emiss_region)
  CASE (303)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(acure_ff_emiss_region(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,acure_ff_emiss_region)
  CASE (304)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(acure_res_emiss_region(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,acure_res_emiss_region)
```

#### Stage 3 - Deallocate the arrays

Without explicit deallocation, a new array will be created every timestep leading to a memopry leak. To prevent this, the following lines need to be added to the end of the same `ukca_um_interf.F90` file, somewhere before the call to ntp_dealloc at the end of the `ukca_um_dealloc` subroutine:

```
IF (ALLOCATED(acure_anth_so2_region)) DEALLOCATE(acure_anth_so2_region)
IF (ALLOCATED(acure_bb_emiss_region)) DEALLOCATE(acure_bb_emiss_region)
IF (ALLOCATED(acure_ff_emiss_region)) DEALLOCATE(acure_ff_emiss_region)
IF (ALLOCATED(acure_res_emiss_region)) DEALLOCATE(acure_res_emiss_region)
```

#### Stage 4 - Set the new variables and their dimensions

In the `ukca_setd1defs` file, there are a few different additions needed. The first is to ensure that the logicals that determine whether or not the regional perturbations are needed are passed over from ukca_option_mod by expanding the list of parameters passed over to include them as below:
```
USE ukca_option_mod, ONLY: L_ukca_rado3, L_ukca_radch4,              &
                           L_ukca_mode, L_ukca_chem, L_ukca_dust,    &
                           L_ukca_qch4inter, i_ukca_photol,          &
                           L_ukca_arg_act, L_ukca_radn2o,            &
                           L_ukca_ageair, L_ukca_achem,              &
                           L_ukca_aerchem, L_ukca_tropisop,          &
                           L_ukca_trop, l_ukca_prim_moc,             &
                           L_ukca_stratcfc, L_ukca_trophet,          &
                           L_ukca_strattrop, L_ukca_raq,             &
                           L_ukca_strat, jpctr, l_ukca,              &
                           tr_ukca_a,                                &
                           l_ukca_offline, l_ukca_offline_be,        &
                           l_ukca_nr_aqchem,                         &
                           l_ukca_classic_hetchem,                   &
                           l_ukca_raqaero, l_ukca_chem_aero,         &
                           l_acure_anth_so2, l_acure_carb_bb_ems,    &
                           l_acure_carb_ff_ems, l_acure_carb_res_ems
```

The next modification is to increase the number of prognostic fields needed in the UkcaD1codes array. This is in defined in a code block starting on line 427 in the unmodified file:
```
!n_in_progs     =  40     ! max no of progs reqd other than tracers/ems !!!ORIGINAL LINE
n_in_progs     =  44     ! max no of progs reqd other than tracers/ems  !!!WITH NEW FIELDS
n_in_diags0    =   4     ! max no of diags (sect 0) reqd
n_in_diags1    =   4     ! max no of diags (sect 1) reqd
n_in_diags2    =   1     ! max no of diags (sect 2) reqd
n_in_diags3    =  16     ! max no of diags (sect 3) reqd
n_in_diags4    =  10     ! max no of diags (sect 4) reqd
n_in_diags5    =   4     ! max no of diags (sect 5) reqd
n_in_diags8    =   1     ! max no of diags (sect 8) reqd
n_in_diags15   =   1     ! max no of diags (sect 15) reqd
n_in_diags30   =   1     ! max no of diags (sect 30) reqd
```

As seen above the number of prognostics is changed from 40 to 44. In the UMDP304 documentation, the older version of the UM described therein has 41 prognostics.

Now that the correct number of prognostics is added, the items can be added to the UkcaD1codes array. Starting after the `j+40` index, the following can be added:
```
UkcaD1Codes(j+41)%item=301           ! ACURE Anth SO2 regions - CCS ADDITION
UkcaD1Codes(j+41)%len_dim1=row_length
UkcaD1Codes(j+41)%len_dim2=rows
IF (.NOT. l_acure_anth_so2) UkcaD1Codes(j+41)%required=.FALSE.
UkcaD1Codes(j+42)%item=302           ! ACURE BB Emissions Regions - CCS ADDITION
UkcaD1Codes(j+42)%len_dim1=row_length
UkcaD1Codes(j+42)%len_dim2=rows
IF (.NOT. l_acure_carb_bb_ems) UkcaD1Codes(j+42)%required=.FALSE.
UkcaD1Codes(j+43)%item=303           ! ACURE FF Emissions Regions - CCS ADDITION
UkcaD1Codes(j+43)%len_dim1=row_length
UkcaD1Codes(j+43)%len_dim2=rows
IF (.NOT. l_acure_carb_ff_ems) UkcaD1Codes(j+43)%required=.FALSE.
UkcaD1Codes(j+44)%item=304           ! ACURE RES Emissions Regions - CCS ADDITION
UkcaD1Codes(j+44)%len_dim1=row_length
UkcaD1Codes(j+44)%len_dim2=rows
IF (.NOT. l_acure_carb_res_ems) UkcaD1Codes(j+44)%required=.FALSE.
```

This will add the regions to the UkcaD1Codes array at index `j=dtn_dim + n_dust_emissions + 41`. The presence of the prognostics can be confirmed by setting the PRINT_STATUS in `um/env/Runtime Controls/Atmosphere only` in the rose GUI to "Extra diagnostic messages" and examining the output from the cylc run in the `cylc-run/<suite-id>/log/job/<date stamp>/atmos_main` folder on puma. It should be noted however that an issue caused by a write statement in EasyAerosol and triggered at more verbose levels of output will cause the atmos_main task to crash. This can be fixed in the file `src/atmosphere/radiation_control/easyaerosol_read_input_mod.F90:1497` with the following change:

```
--          DO k = 1, dimsize(4)

++          DO k = 1, dimsize(3)
```
