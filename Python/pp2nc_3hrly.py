#!/usr/bin/env python2.7
"""
Script name: pp2nc_3hrly.py
Author: James O'Neill. Based on very helpful python scripts provided by Masaru Yoshioka.
Date: March 2018
Purpose: Extract and condense information from 3-hourly pp files belonging to the UKCA26AER
         perturbed parameter ensemble (PPE) set into netCDF format. Each netCDF file contains
         aerosol optical depth (AOD) at 550nm and column-integrated cloud droplet number
         concentration (CDNC) fields for all 235 PPE members for one day (8 timesteps per file)
         on a coarsened grid (N96 down to N48)
Usage: ./pp2nc_3hrly.py <ppRoot> <ncRoot> <ncRef> <startDate> <endDate>
        <ppRoot> - path (either relative to the current directory, or full) to the root directory containing the pp files
        <ncRoot> - path (relative or full) to the desired output directory
        <ncRef> - path (relative or full) to a reference netCDF file whose coordinate system is at the desired coarsened resolution (N48) onto which the original high-resolution (N96) data should be regridded
        <startDate> - in the format YYYYMMDD and refers to the first day to be processed
        <endDate> - in the format YYYYMMDD and refers to the last day (inclusive) to be processed
Output: Daily netCDF files (8 timesteps per file) with the naming convention <ncRoot>/pbYYYYMMDD.nc
"""

def readAOD():
    """
    Try reading in all 6 cubes relating to AOD, check dimensions are as expected, and sum to get total AOD.
    Return with error flag if something goes wrong, so that this file can be skipped by the main routine.
    """
    AODcube=0
    err=False
    for s,stashid in enumerate(stash_AOD):
        try:
            modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stashid))
        except:
            err=True
            flog.write('Warning: Could not load AOD data from input file '+ppFile+'. File has been skipped\n')
            flog.close
            break
        if modeCube.shape != (nTimes,N96res[0],N96res[1]):
            err=True
            flog.write('Warning: Unexpected AOD dimensions in input file '+ppFile+ \
              '. Expected '+str((nTimes,N96res[0],N96res[1]))+ \
              ', got '+str(modeCube.shape)+'. File has been skipped\n')
            flog.close
            break
        if s==0:
            AODcube=modeCube
        else:
            AODcube+=modeCube
    return(AODcube,err)


import argparse
import os
#####READ IN COMMAND LINE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument("ppRoot",help="absolute/relative path to root directory containing the pp files",type=str)
parser.add_argument("ncRoot",help="absolute/relative path to desired output directory",type=str)
parser.add_argument("ncRef",help="absolute/relative path to reference nc file at the desired coarsened (N48) resolution",type=str)  
parser.add_argument("startDate",help="First day to be processed, in format YYYYMMDD, e.g. 20080701",type=str)
parser.add_argument("endDate",help="Last day to be processed, in format YYYYMMDD, e.g. 20080710",type=str)
args = parser.parse_args()
ppRoot=args.ppRoot
ncRoot=args.ncRoot
ncRef=args.ncRef
startDate=args.startDate
endDate=args.endDate
#####

#####CHECK FILES/DIRS exist
assert os.path.exists(ppRoot), "ppRoot does not exist"
assert os.path.exists(ncRoot), "ncRoot does not exist"
assert os.path.exists(ncRef),"ncRef does not exist"
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
maxDays = 31 #Maximum number of days that can be processed in one go
N96res=(145,192) #Expected number of lats/lons in finer resolution data
N48res=(73,96) #Expected number of lats/lons in coarser resolution data
nTimes=8 #Expected number of times in pp files (one day at 3-hrly resolution per file)
nLevels=52 #Expected number vertical levels in pp files
missing=-999 #Missing value written to nc file if pp file is skipped
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
assert numDays > 0, "Error: startDate must come before endDate"
assert numDays <= maxDays, "Error: Can only process up to 31 days of data at a time"
fileDates=[]
for date in [startDateObj + dt.timedelta(n) for n in range(numDays)]:
    #The date.strftime() functions below ensure e.g. that the day is written as '01' and not '1'
    fileDates.append(date.strftime('%Y')+date.strftime('%m')+date.strftime('%d'))

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
nFilesRead=0
flog=open(ncRoot+'/logfile.log','w')
flog.close
#Loop over days:
for day in fileDates:
    #open and close log file to allow reading of it during execution
    flog=open(ncRoot+'/logfile.log','a')
    flog.write("---------------------\n")
    flog.write("Processing day: "+day+'\n')
    flog.close
    outFile=ncRoot+'/pb'+day+'.nc' #path to output nc file
    #Initialise big arrays (all jobs for current day):
    CDNC_all=np.full(shape=(nTimes,N48res[0],N48res[1],len(jobids)),fill_value=missing,dtype=np.float64)
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
        AODcube,err = readAOD()
        if err==True: continue
        nFilesRead+=1
        #If first file read in, store time and level information:
        if nFilesRead==1:
            try:
                CDNCtimes = CDNCcube.coord('time').points
                hgts = CDNCcube.coord('level_height').points
                sigmas = CDNCcube.coord('sigma').points
            except:
                print("Error: Unexpected coordinates for CDNC in pp file "+ppFile)
                raise
            try:
                modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stash_AOD[0]))
                AODtimes = modeCube.coord('time').points
            except:
                print("Error: Unexpected time coordinates for AOD in pp file "+ppFile)
                raise
        #***FOR NOW*** just take arithmetic mean of CDNC over all heights
        for h in range(nLevels):
            if h==0:
                CDNCsum=CDNCcube[:,h,:,:]
            else:
                CDNCsum+=CDNCcube[:,h,:,:]
        CDNCave=CDNCsum/float(nLevels)
        #Regrid onto coarser resolution (N48):
        CNDC_N48 = CDNCave.regrid(refCube, iris.analysis.Linear())
        AOD_N48 = AODcube.regrid(refCube, iris.analysis.Linear())
        #Update big arrays (all jobs for current day)
        CDNC_all[:,:,:,j]=CNDC_N48.data
        AOD_all[:,:,:,j]=AOD_N48.data
        flog.close

    #Create nc file:
    assert not (np.max(CDNC_all) == missing and np.min(CDNC_all) == missing), "No processed data"
    flog=open(ncRoot+'/logfile.log','a')
    flog.write("Writing to nc file: "+outFile+'\n')
    ncfile = netCDF4.Dataset(outFile,mode='w',format='NETCDF4_CLASSIC')
    #Add time dimension/variables (times are potentially different for CDNC and AOD):
    tCDNC_dim = ncfile.createDimension('time_CDNC',len(CDNCtimes))
    tCDNC_out = ncfile.createVariable('time_CDNC', np.float64, ('time_CDNC'))
    tCDNC_out[:] = CDNCtimes
    tCDNC_out.units = "hours since 1970-01-01 00:00:00"
    tCDNC_out.calendar = "gregorian"
    #
    tAOD_dim = ncfile.createDimension('time_AOD',len(AODtimes))
    tAOD_out = ncfile.createVariable('time_AOD', np.float64, ('time_AOD'))
    tAOD_out[:] = AODtimes
    tAOD_out.units = "hours since 1970-01-01 00:00:00"
    tAOD_out.calendar = "gregorian"
    #Add lat/lon dimensions/variables:
    lon_dim = ncfile.createDimension('longitude', len(lons))
    lat_dim = ncfile.createDimension('latitude', len(lats))
    lon_out = ncfile.createVariable('longitude', np.float64, ('longitude'))
    lat_out = ncfile.createVariable('latitude', np.float64, ('latitude'))
    lon_out[:] = lons
    lat_out[:] = lats
    lon_out.units="degrees_east"
    lat_out.units="degrees_north"
    #Add job dimension/variable. Note that multi-character strings must be handled using stringtochar:
    job_dim = ncfile.createDimension('job',len(jobids))
    nchar_dim = ncfile.createDimension('nchar',5)
    str_out = netCDF4.stringtochar(np.array(jobids,'S5'))
    job_out = ncfile.createVariable('job', 'S1', ('job','nchar'))
    job_out[:,:] = str_out
    job_out.units = "none"
    #Add AOD550 variable:
    aod_out = ncfile.createVariable('AOD550_total', np.float64, ('time_AOD','latitude','longitude','job'))
    aod_out[:,:,:,:] = AOD_all
    aod_out.units = "1"
    #Add CDNC variable:
    cdnc_out = ncfile.createVariable('CDNC', np.float64, ('time_CDNC','latitude','longitude','job'))
    cdnc_out[:,:,:,:] = CDNC_all
    cdnc_out.units = "number per metre cubed"
    #close netCDF file
    ncfile.close()
    flog.close
