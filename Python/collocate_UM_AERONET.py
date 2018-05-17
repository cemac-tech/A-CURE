#!/usr/bin/env python2.7
import string
from calendar import monthrange
import datetime as dt
import os
import iris
import numpy as np
import netCDF4

#####PARAMETERS
month='200807' #Month in format YYYYMM
ppRoot='/group_workspaces/jasmin2/ukca/vol1/myoshioka/um/Dumps' #pp root directory
ncRoot='/home/users/earjjo/GASSP_WS/aod550_total_jobid_pbYYYYMMDD_nc' #Output directory for pp->nc files
stash_AOD = ['m01s02i500','m01s02i501','m01s02i502','m01s02i503','m01s02i504','m01s02i505'] #STASH codes for 6 'modes' in pp files. AOD_550 is sum over all these modes
missing=np.nan #Missing value written to nc file if field is missing
run_pp2nc=True
#####

#####GENERATE ARRAY WITH MIDDLE PART OF FILENAME CONVENTION
#Expected filename format under root directory is: <jobid>/<jobid>a.pbYYYYDDMM.pp
numDaysInMon=monthrange(int(month[0:4]),int(month[4:6]))[1]
fileDates=[]
#####TEMPORARILY ONLY RUN FOR 2 DAYS RATHER THAN WHOLE MONTH:#####
#for date in [dt.date(int(month[0:4]),int(month[4:6]),1) + dt.timedelta(n) for n in range(numDaysInMon)]:
for date in [dt.date(int(month[0:4]),int(month[4:6]),1) + dt.timedelta(n) for n in range(2)]:
    #The date.strftime() functions below ensure e.g. that the day is written as '01' and not '1'
    fileDates.append('pb'+date.strftime('%Y')+date.strftime('%m')+date.strftime('%d'))
#####

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

#####GENERATE NC FILES FROM PP FILES
if(run_pp2nc):
    for date in fileDates:
        for j,jobid in enumerate(jobids):
            outFile=ncRoot+'/aod550_total_'+jobid+'_'+date+'.nc' #path to output nc file
            #Generate full PP file path
            ppFile=jobid+'a.'+date+'.pp' #specific pp file name
            ppPath=ppRoot+'/'+jobid+'/'+ppFile #specific pp file path
            #Check input file exists, skip if not
            if not os.path.exists(ppPath):
                print('Warning: Input file '+ppFile+' was not found. File has been skipped')
                continue
            #Read in all 6 cubes relating to AOD; sum to get total AOD
            AODcube=0
            for s,stashid in enumerate(stash_AOD):
                modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stashid))
                if s==0:
                    AODcube=modeCube
                    AODtimes=modeCube.coord('time').points
                    AODlats=modeCube.coord('latitude').points
                    AODlons=modeCube.coord('longitude').points
                else:
                    AODcube+=modeCube
            ####Create AOD nc file:
            print("Writing AOD fields to nc file: "+outFile)
            ncfile = netCDF4.Dataset(outFile,mode='w',format='NETCDF4_CLASSIC')
            #Add time dimension/variable:
            tAOD = ncfile.createDimension('time',len(AODtimes))
            tAOD_out = ncfile.createVariable('time', np.float64, ('time'))
            tAOD_out[:] = AODtimes
            tAOD_out.long_name = "time"
            tAOD_out.standard_name = "time"
            tAOD_out.units = "hours since 1970-01-01 00:00:00"
            tAOD_out.calendar = "gregorian"
            #Add lat/lon dimensions/variables:
            lon_dim = ncfile.createDimension('longitude', len(AODlons))
            lat_dim = ncfile.createDimension('latitude', len(AODlats))
            lon_out = ncfile.createVariable('longitude', np.float64, ('longitude'))
            lat_out = ncfile.createVariable('latitude', np.float64, ('latitude'))
            lon_out[:] = AODlons
            lat_out[:] = AODlats
            lon_out.long_name="longitude"
            lat_out.long_name="latitude"
            lon_out.standard_name="longitude"
            lat_out.standard_name="latitude"
            lon_out.units="degrees_east"
            lat_out.units="degrees_north"
            #Add AOD550 variable:
            aod_out = ncfile.createVariable('aod550', np.float64, ('time','latitude','longitude'),fill_value=missing)
            aod_out[:,:,:] = AODcube.data
            aod_out.long_name="Aerosol optical depth at 550nm"
            aod_out.units = "none"
            aod_out.stash_codes="summation of m01s02i500, m01s02i501, m01s02i502, m01s02i503, m01s02i504, m01s02i505"
            aod_out.stash_names="AITKENMODESOLUBLEOPTDEPNONADV, ACCUMMODESOLUBLEOPTDEPNONADV, COARSEMODESOLUBLEOPTDEPNONADV, \
            AITKENMODEINSOLOPTDEPNONADV, ACCUMMODEINSOLOPTDEPNONADV, COARSEMODEINSOLOPTDEPNONADV"
            #close netCDF file
            ncfile.close()
#####