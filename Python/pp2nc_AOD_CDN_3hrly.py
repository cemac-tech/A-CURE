#!/usr/bin/env python2.7
"""
Script name: pp2nc_3hrly.py
Author: James O'Neill. Based on very helpful python scripts provided by Masaru Yoshioka as well as
        UCKA example scripts pointed to by Kirsty Pringle
Date: March 2018
Purpose: Extract and condense information from 3-hourly pp files belonging to the UKCA26AER
         perturbed parameter ensemble (PPE) set into netCDF format. One set of netCDF files
         contains aerosol optical depth (AOD) at 550nm, and another set contain column-integrated
         cloud droplet number concentration (CDNC) fields for all 235 PPE members for one day
         (8 timesteps per file) on a coarsened grid (N96 down to N48)
Usage: ./pp2nc_AOD_CDN_3hrly.py <ppRoot> <orogFile> <ncRef> <ncRoot> <startDate> <endDate>
        <ppRoot> - path (either relative to the current directory, or full) to the root directory
                   containing the pp files
        <orogFile> - path (relative or full) to ancilliary UM file containing the orgraphy data
                     (file typically called 'qrparm.orog')
        <ncRef> - path (relative or full) to a reference netCDF file whose coordinate system is at
                  the desired coarsened resolution (N48) onto which the original high-resolution (N96)
                  data should be regridded
        <ncRoot> - path (relative or full) to the desired output directory
        <startDate> - in the format YYYYMMDD and refers to the first day to be processed
        <endDate> - in the format YYYYMMDD and refers to the last day (inclusive) to be processed
Output: Daily netCDF files (8 timesteps per file) with the naming convention 
        <ncRoot>/aod_tebaa-tebiz_teafw_pbYYYYMMDD_N48.nc (for AOD), and
        <ncRoot>/cdn_tebaa-tebiz_teafw_pbYYYYMMDD_N48.nc (for CDN)
"""

def readAOD():
    """
    Try reading in all 6 cubes relating to AOD, check dimensions are as expected, and sum to get total AOD.
    Return with error flag if something goes wrong, so that this file can be skipped by the main routine.
    """
    AODcube=0
    err=False
    only7timesteps=False
    for s,stashid in enumerate(stash_AOD):
        #Try loading in each mode's cube:
        try:
            modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stashid))
        except:
            err=True
            flog.write('Warning: Could not load AOD data from input file '+ppFile+'. File has been skipped\n')
            flog.close
            break
        #Save the shape of the first mode, or check shape is the same as the first mode
        if s==0:
            mode1shape=modeCube.shape
        else:
            if modeCube.shape != mode1shape:
                err=True
                flog.write('Warning: Inconsistent AOD mode dimensions in input file '+ppFile+'. File has been skipped\n')
                flog.close
                break
        #Sum over all modes:
        if s==0:
            AODcube=modeCube
        else:
            AODcube+=modeCube
    #Check the dimensions are as expected:
    if not err and AODcube.shape != (nTimes,N96res[0],N96res[1]):
        #First check the case where there are only 7 rather than 8 timesteps.
        #Masaru said that this sometimes happens on the first day of the month,
        #and that in this case, we should specify the first timestep as missing
        #but fill the other 7 timesteps:
        if AODcube.shape == (nTimes-1,N96res[0],N96res[1]):
            only7timesteps=True
            flog.write('Warning: only 7 rather than 8 AOD timesteps in file '+ppFile+ \
                   '. Assuming first timestep is missing and processing remaining 7 timesteps\n')
            flog.close
        else:
            err=True
            flog.write('Warning: Unexpected AOD dimensions in input file '+ppFile+ \
              '. Expected '+str((nTimes,N96res[0],N96res[1]))+ \
              ', got '+str(AODcube.shape)+'. File has been skipped\n')
            flog.close
    #return to main routine:
    return(AODcube,err,only7timesteps)

import argparse
import os
#####READ IN COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser(description="Script to extract and condense information from 3-hourly pp files belonging to \
         the UKCA26AER perturbed parameter ensemble set into netCDF format. Each resulting netCDF file contains aerosol optical \
         depth at 550nm and column-integrated cloud droplet number concentration fields for all 235 PPE members for one day \
         (8 timesteps per file) on a coarsened grid (N96 down to N48)")
parser.add_argument("ppRoot",help="absolute/relative path to root directory containing the pp files",type=str)
parser.add_argument("orogFile",help="absolute/relative path to orography data file",type=str) 
parser.add_argument("ncRef",help="absolute/relative path to reference nc file at the desired coarsened (N48) resolution",type=str)
parser.add_argument("ncRoot",help="absolute/relative path to desired output directory",type=str) 
parser.add_argument("startDate",help="First day to be processed, in format YYYYMMDD, e.g. 20080701",type=str)
parser.add_argument("endDate",help="Last day to be processed, in format YYYYMMDD, e.g. 20080710",type=str)
args = parser.parse_args()
ppRoot=args.ppRoot
orogFile=args.orogFile
ncRef=args.ncRef
ncRoot=args.ncRoot
startDate=args.startDate
endDate=args.endDate
#####

#####CHECK FILES/DIRS exist
assert os.path.exists(ppRoot), "ppRoot does not exist"
assert os.path.exists(orogFile), "orogFile does not exist"
assert os.path.exists(ncRef), "ncRef does not exist"
assert os.path.exists(ncRoot), "ncRoot does not exist"
#####

import iris
import numpy as np
import string
import datetime as dt
import netCDF4
from dateutil.parser import parse

#####PARAMETERS
stash_CDNC = 'm01s38i479' #STASH code for CDNC field in pp files
stash_AOD = ['m01s02i500','m01s02i501','m01s02i502','m01s02i503','m01s02i504','m01s02i505'] #STASH codes for 6 'modes' in pp files. AOD_550 is sum over all these modes
stash_orog = 'm01s00i033' #STASH code for orography data
maxDays = 31 #Maximum number of days that can be processed in one go
N96res=(145,192) #Expected number of lats/lons in finer resolution data
N48res=(73,96) #Expected number of lats/lons in coarser resolution data
nTimes=8 #Expected number of times in pp files (one day at 3-hrly resolution per file)
nLevels=52 #Expected number of vertical levels in pp files
missing=np.nan #Missing value written to nc file if field is missing
######

#####GENERATE JOBID ARRAY
#In this particular PPE (named 'UKCA26AER'), there are 182 'training' jobs,
#one 'median' job and 52 'validation' jobs, making 235 jobs (or simulations
#or members) altogether. The jobs must be in this order. The 'training' and
#'validation' jobs all have job ids that start with the string 'teb' and the
#last two characters loop through the alphabet (final character loops the
#fastest) starting from 'aa' and ending at 'iz'. The median job id is 'teafw'
jobids=[]
#Training jobs first:
for c1 in string.ascii_lowercase[0:7]:
    for c2 in string.ascii_lowercase:
      jobids.append('teb'+c1+c2)
#median job next:
jobids.append('teafw')
#validation jobs last:
for c1 in string.ascii_lowercase[7:9]:
    for c2 in string.ascii_lowercase:
      jobids.append('teb'+c1+c2)
#####

#####GENERATE ARRAY WITH DATE PART OF FILENAME CONVENTION
#Expected filename format under root directory is: <jobid>/<jobid>a.pbYYYYDDMM.pp
#The 'pb' indicates 3-hrly output runs ('pm', 'pa' and 'pc' indicate monthly,
#daily and 1hrly, respectively).
try:
    startDateObj=parse(startDate)
except:
    print("Error: startDate input parameter in unexpected format")
    raise
try:
    endDateObj=parse(endDate)
except:
    print("Error: endDate input parameter in unexpected format")
    raise
numDays=(endDateObj-startDateObj).days + 1 #Add one to include both extremeties
assert numDays > 0, "startDate must come before endDate"
assert numDays <= maxDays, "Can only process up to 31 days of data at a time"
fileDates=[]
for date in [startDateObj + dt.timedelta(n) for n in range(numDays)]:
    #The date.strftime() functions below ensure e.g. that the day is written as '01' and not '1'
    fileDates.append(date.strftime('%Y')+date.strftime('%m')+date.strftime('%d'))

#####LOAD AND STORE OROGRAPHY DATA
try:
    orog=iris.load_cube(orogFile,iris.AttributeConstraint(STASH=stash_orog))
    #Create auxiliary coordinate for orography:
    auxcoord=iris.coords.AuxCoord(orog.data,standard_name=str(orog.standard_name),
                      long_name="orography",var_name="orog",units=orog.units)
except:
    print("Error: Could not load in orography file "+orogFile)
    raise
assert orog.shape == N96res, "Unexpected dimensions in orography data file"
#####

#####LOAD AND STORE REFERENCE NC FILE DATA
try:
    refCube=iris.load(ncRef)[0]
    lons = refCube.coord('longitude').points
    lats = refCube.coord('latitude').points
except:
    print("Error: Couldn't load or extract coordinates from reference data file")
    raise
assert len(lats) == N48res[0] and len(lons) == N48res[1], "Unexpected dimensions in reference data file"
#####

#####FOR EACH DAY, LOOP OVER ALL JOBS, STORE DATA IN ONE BIG ARRAY, AND OUTPUT TO NC FILE
nValidFiles=0
flog=open(ncRoot+'/logfile.log','w')
flog.close
#Loop over days:
for day in fileDates:
    #open and close log file to allow reading of it during execution
    flog=open(ncRoot+'/logfile.log','a')
    flog.write("---------------------\n")
    flog.write("Processing day: "+day+'\n')
    flog.close
    outFileAOD=ncRoot+'/aod550_tebaa-tebiz_teafw_pb'+day+'_N48.nc' #path to AOD output nc file
    outFileCDN=ncRoot+'/cdn_tebaa-tebiz_teafw_pb'+day+'_N48.nc' #path to CDN output nc file
    #Initialise big arrays (all jobs for current day):
    CDNpm2_all=np.full(shape=(nTimes,N48res[0],N48res[1],len(jobids)),fill_value=missing,dtype=np.float64)
    AOD_all=np.full(shape=(nTimes,N48res[0],N48res[1],len(jobids)),fill_value=missing,dtype=np.float64)
    #Loop over jobs:
    for j,jobid in enumerate(jobids):
        #Open log file here, close it at end of loop to allow reading of it during execution
        flog=open(ncRoot+'/logfile.log','a')
        flog.write("Processing job: "+jobid+" ("+str(j+1)+" of "+str(len(jobids))+")\n")
        #Generate full PP file path
        ppFile=jobid+'a.pb'+day+'.pp' #specific pp file name
        ppPath=ppRoot+'/'+jobid+'/'+ppFile #specific pp file path
        #Check input file exists, skip if not
        if not os.path.exists(ppPath):
            flog.write('Warning: Input file '+ppFile+'was not found. File has been skipped\n')
            flog.close
            continue
        #Try reading in CDNC cube and check dimensions are as expected
        try:
            CDNCcube = iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stash_CDNC))
        except:
            flog.write('Warning: Could not load CDNC data from input file '+ppFile+'. File has been skipped\n')
            flog.close
            continue
        if CDNCcube.shape != (nTimes,nLevels,N96res[0],N96res[1]):
            flog.write('Warning: Unexpected CDNC dimensions in input file '+ppFile+ \
              '. Expected '+str((nTimes,nLevels,N96res[0],N96res[1]))+ \
              ', got '+str(CDNCcube.shape)+'. File has been skipped\n')
            flog.close
            continue
        #Try reading in all 6 cubes relating to AOD, check dimensions are as expected, and sum to get total AOD
        AODcube,err,only7timesteps = readAOD()
        if err==True: continue
        if not only7timesteps: nValidFiles+=1
        #If first file read in with 8 timesteps, store time information:
        if nValidFiles==1:
            try:
                CDNCtimes = CDNCcube.coord('time').points
            except:
                print("Error: Unexpected coordinates for CDNC in pp file "+ppFile)
                raise
            try:
                modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stash_AOD[0]))
                AODtimes = modeCube.coord('time').points
            except:
                print("Error: Unexpected time coordinates for AOD in pp file "+ppFile)
                raise
                
        ###Calculate column-integrated CDNC using IRIS method:
        CDNCcube.add_aux_coord(auxcoord,(2,3,)) #Add auxiliary orography coordinate to cube on lat/lon grid
        #Create a hybrid-height coordinate factory with the formula z = a + b * orog:
        factory=iris.aux_factory.HybridHeightFactory(delta=CDNCcube.coord("level_height"),
                sigma=CDNCcube.coord("sigma"),orography=CDNCcube.coord("surface_altitude")) 
        CDNCcube.add_aux_factory(factory) #Add factory to cube to create altitude derived coordinate
        #Calculate depth of each grid cell as difference between altitude bounds:
        depths = CDNCcube.coord('altitude').bounds[:,:,:,1] - CDNCcube.coord('altitude').bounds[:,:,:,0]
        #Multiply each CDNC data point by the depth of its cell to give CDN per m^2 through that cell:
        CDNCcube.data = CDNCcube.data * depths
        #Sum all data from each column to give CDN per m^2 through entire atmosphere
        CDNCcube_int=CDNCcube.collapsed('model_level_number',iris.analysis.SUM)
        
        #Regrid onto coarser resolution (N48):
        #Note: Currently using IRIS's bilinear regridding method. Discussion is ongoing as to whether this 
        #is the most appropriate method.
        CNDpm2_N48 = CDNCcube_int.regrid(refCube, iris.analysis.Linear())
        AOD_N48 = AODcube.regrid(refCube, iris.analysis.Linear())
        #Update big arrays (all jobs for current day)
        CDNpm2_all[:,:,:,j]=CNDpm2_N48.data
        if not only7timesteps:
            AOD_all[:,:,:,j]=AOD_N48.data
        else:
            AOD_all[1:nTimes,:,:,j]=AOD_N48.data
        #Close log file:
        flog.close

    assert nValidFiles > 0, "No processed data"
    flog=open(ncRoot+'/logfile.log','a')
    ####Create AOD nc file:
    flog.write("Writing AOD fields to nc file: "+outFileAOD+'\n')
    ncfile = netCDF4.Dataset(outFileAOD,mode='w',format='NETCDF4_CLASSIC')
    #Add time dimension/variable:
    tAOD = ncfile.createDimension('time',len(AODtimes))
    tAOD_out = ncfile.createVariable('time', np.float64, ('time'))
    tAOD_out[:] = AODtimes
    tAOD_out.long_name = "time"
    tAOD_out.standard_name = "time"
    tAOD_out.units = "hours since 1970-01-01 00:00:00"
    tAOD_out.calendar = "gregorian"
    #Add lat/lon dimensions/variables:
    lon_dim = ncfile.createDimension('longitude', len(lons))
    lat_dim = ncfile.createDimension('latitude', len(lats))
    lon_out = ncfile.createVariable('longitude', np.float64, ('longitude'))
    lat_out = ncfile.createVariable('latitude', np.float64, ('latitude'))
    lon_out[:] = lons
    lat_out[:] = lats
    lon_out.long_name="longitude"
    lat_out.long_name="latitude"
    lon_out.standard_name="longitude"
    lat_out.standard_name="latitude"
    lon_out.units="degrees_east"
    lat_out.units="degrees_north"
    #Add job dimension/variable. Note that multi-character strings must be handled using stringtochar:
    job_dim = ncfile.createDimension('job',len(jobids))
    nchar_dim = ncfile.createDimension('nchar',5)
    str_out = netCDF4.stringtochar(np.array(jobids,'S5'))
    job_out = ncfile.createVariable('job', 'S1', ('job','nchar'))
    job_out[:,:] = str_out
    job_out.long_name = "job index of UKCA26AER PPE member"
    job_out.units = "none"
    #Add AOD550 variable:
    aod_out = ncfile.createVariable('aod550', np.float64, ('time','latitude','longitude','job'),fill_value=missing)
    aod_out[:,:,:,:] = AOD_all
    aod_out.long_name="Aerosol optical depth at 550nm"
    aod_out.units = "none"
    aod_out.stash_codes="summation of m01s02i500, m01s02i501, m01s02i502, m01s02i503, m01s02i504, m01s02i505"
    aod_out.stash_names="AITKENMODESOLUBLEOPTDEPNONADV, ACCUMMODESOLUBLEOPTDEPNONADV, COARSEMODESOLUBLEOPTDEPNONADV, \
    AITKENMODEINSOLOPTDEPNONADV, ACCUMMODEINSOLOPTDEPNONADV, COARSEMODEINSOLOPTDEPNONADV"
    #close netCDF file
    ncfile.close()
    
    ####Create CDN nc file:
    flog.write("Writing CDN fields to nc file: "+outFileCDN+'\n')
    ncfile = netCDF4.Dataset(outFileCDN,mode='w',format='NETCDF4_CLASSIC')
    #Add time dimension/variables:
    tCDN_dim = ncfile.createDimension('time',len(CDNCtimes))
    tCDN_out = ncfile.createVariable('time', np.float64, ('time'))
    tCDN_out[:] = CDNCtimes
    tCDN_out.long_name = "time"
    tCDN_out.standard_name = "time"
    tCDN_out.units = "hours since 1970-01-01 00:00:00"
    tCDN_out.calendar = "gregorian"
    #Add lat/lon dimensions/variables:
    lon_dim = ncfile.createDimension('longitude', len(lons))
    lat_dim = ncfile.createDimension('latitude', len(lats))
    lon_out = ncfile.createVariable('longitude', np.float64, ('longitude'))
    lat_out = ncfile.createVariable('latitude', np.float64, ('latitude'))
    lon_out[:] = lons
    lat_out[:] = lats
    lon_out.long_name="longitude"
    lat_out.long_name="latitude"
    lon_out.standard_name="longitude"
    lat_out.standard_name="latitude"
    lon_out.units="degrees_east"
    lat_out.units="degrees_north"
    #Add job dimension/variable. Note that multi-character strings must be handled using stringtochar:
    job_dim = ncfile.createDimension('job',len(jobids))
    nchar_dim = ncfile.createDimension('nchar',5)
    str_out = netCDF4.stringtochar(np.array(jobids,'S5'))
    job_out = ncfile.createVariable('job', 'S1', ('job','nchar'))
    job_out[:,:] = str_out
    job_out.long_name = "job index of UKCA26AER PPE member"
    job_out.units = "none"
    #Add CDN per m^2 variable:
    cdn_out = ncfile.createVariable('cdn', np.float64, ('time','latitude','longitude','job'),fill_value=missing)
    cdn_out[:,:,:,:] = CDNpm2_all
    cdn_out.long_name="Column integrated cloud droplet number concentration (i.e. cloud droplet number per m^2)"
    cdn_out.units = "m-2"
    cdn_out.stash_code="column-integrated m01s38i479 per m^2"
    #close netCDF file
    ncfile.close()
    flog.close
