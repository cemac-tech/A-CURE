#!/usr/bin/env python2.7
import string
from calendar import monthrange
import datetime as dt
import os
import iris
import numpy as np
import netCDF4
import glob
import cis
import re
import pandas as pd

assert "CIS_PLUGIN_HOME" in os.environ, "Environment variable CIS_PLUGIN_HOME not set"
assert os.path.exists(os.path.join(os.getenv("CIS_PLUGIN_HOME"),"cis_plugin_AERONETv3nc.py")), "Cannot find cis_plugin_AERONETv3nc.py in CIS_PLUGIN_HOME directory"

#####PARAMETERS
YYYYMM='200807' #Month in format YYYYMM
ppRoot='/group_workspaces/jasmin2/ukca/vol1/myoshioka/um/Dumps' #pp root directory
AERONETRoot='/group_workspaces/jasmin2/crescendo/Data/AERONET/AOT/ver3/LEV20/Monthly' #AERONET root directory
UMRoot='/home/users/earjjo/GASSP_WS/aod550_total_jobid_pbYYYYMMDD_nc' #Output directory for pp->nc files
colRoot='/home/users/earjjo/GASSP_WS/aod550_total_jobid_pbYYYYMM_station_nc' #Output directory for collocated nc files
mavRoot='/home/users/earjjo/GASSP_WS/aod550_total_jobid_pbYYYYMM_station_mav_nc' #Output directory for monthly averaged nc files
outRoot='/home/users/earjjo/GASSP_WS/aod550_total_jobid_pbYYYYMM_mav_nc' #Output directory for concatenated files
stash_AOD = ['m01s02i500','m01s02i501','m01s02i502','m01s02i503','m01s02i504','m01s02i505'] #STASH codes for 6 'modes' in pp files. AOD_550 is sum over all these modes
missing=np.nan #Missing value written to nc file if field is missing
run_pp2nc=False
#####

#####GENERATE ARRAY WITH MIDDLE PART OF FILENAME CONVENTION
#Expected filename format under root directory is: <jobid>/<jobid>a.pbYYYYDDMM.pp
year=int(YYYYMM[0:4])
mon=int(YYYYMM[4:6])
numDaysInMon=monthrange(year,mon)[1]
numDaysInMon=2 ####TEMPORARY CODE
monStart=dt.datetime(year,mon,1)
monEnd=monStart+dt.timedelta(days=numDaysInMon,seconds=-1)
fileDates=[]
for date in [monStart + dt.timedelta(n) for n in range(numDaysInMon)]:
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
jobids=jobids[0:1] ####TEMPORARY CODE
#####

for j,jobid in enumerate(jobids):
    print('=====PROCESSING JOB '+jobid+'=====')
    
#####GENERATE NC FILES FROM PP FILES
    if(run_pp2nc):
        for date in fileDates:
            outFile=os.path.join(UMRoot,'aod550_total_'+jobid+'_'+date+'.nc') #path to output nc file
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
            print("Extracting AOD_550 from pp file and writing to nc file: "+outFile)
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
            aod_out.units = "1"
            aod_out.stash_codes="summation of m01s02i500, m01s02i501, m01s02i502, m01s02i503, m01s02i504, m01s02i505"
            aod_out.stash_names="AITKENMODESOLUBLEOPTDEPNONADV, ACCUMMODESOLUBLEOPTDEPNONADV, COARSEMODESOLUBLEOPTDEPNONADV, \
            AITKENMODEINSOLOPTDEPNONADV, ACCUMMODEINSOLOPTDEPNONADV, COARSEMODEINSOLOPTDEPNONADV"
            #close netCDF file
            ncfile.close()
#####

#####LOOP THROUGH AERONET FILES (STATIONS); USE CIS TO SUBSET AERONET DATA, COLLOCATE UM DATA AND GET MONTHLY AVERAGES           
    AERONETFilePtn="AOD_440_*_"+YYYYMM[0:4]+"-"+YYYYMM[4:6]+"_v3.nc"
    AERONETPaths=glob.glob(os.path.join(AERONETRoot,AERONETFilePtn))
    AERONETFiles=[os.path.basename(x) for x in AERONETPaths]
    AERONETFiles=AERONETFiles[0:10] ###TEMPORARY CODE
    monAveDF=pd.DataFrame(missing,index=range(0,len(AERONETFiles)),columns=['stn','lat','lon','ave','std','num'])
    monAveDF['stn'] = monAveDF['stn'].astype(str)
    monAveDF['num'] = 0
    monAveDF['num'] = monAveDF['num'].astype(int)
    for i,AERONETFile in enumerate(AERONETFiles):
        underscores=[m.start() for m in re.finditer(r"_",AERONETFile)]
        station=AERONETFile[(underscores[1]+1):underscores[-2]]
        print("Collocating UM data with AERONET data and calculating monthly averages for station "+station)
        #Read in AERONET data with new plugin:
        AERONETData=cis.read_data(os.path.join(AERONETRoot,AERONETFile),"AOD_440",product="AERONETv3nc")
        #Subset AERONET data in time to ensure all points are within UM data time bounds:
        AERONETsubset=AERONETData.subset(time=[monStart,monEnd])
        #Read in relevant UM files:
        UMFilePtn='aod550_total_'+jobid+'_pb'+YYYYMM[0:4]+YYYYMM[4:6]+'??.nc'
        UMFiles=glob.glob(os.path.join(UMRoot,UMFilePtn))
        UMData=cis.read_data(UMFiles,"aod550")
        #Collocate UM data onto AERONET data:
        colData=UMData.collocated_onto(AERONETData,how="lin",var_name="collocated_AOD550",var_units="1")
        #colData.save_data(os.path.join(colRoot,'aod550_total_'+jobid+'_pb'+YYYYMM+'_'+station+'.nc'))
        #Calculate monthly average of collocated aod550:
        aggData=colData.aggregate(how='moments',t=[monStart,monEnd,dt.timedelta(days=numDaysInMon)])
        #aggData.save_data(os.path.join(mavRoot,'aod550_total_'+jobid+'_pb'+YYYYMM+'_'+station+'_mav.nc'))
        monAveDF.at[i,'stn']=station
        monAveDF.at[i,'lat']=aggData[0].get_all_points()[0].latitude
        monAveDF.at[i,'lon']=aggData[0].get_all_points()[0].longitude
        monAveDF.at[i,'ave']=aggData[0].get_all_points()[0].val[0]
        monAveDF.at[i,'std']=aggData[1].get_all_points()[0].val[0]
        monAveDF.at[i,'num']=aggData[2].get_all_points()[0].val[0]
#####SAVE ALL STATION MONTHLY AVERAGES TOGETHER INTO ONE ARRAY
    print(monAveDF)

######SAVE ARRAY TO UNGRIDDED NC FILE
#    ncfilename = os.path.join(outRoot,'aod550_total_'+jobid+'_pb'+YYYYMM+'_mav.nc')
#    ncfile = netCDF4.Dataset(ncfilename,mode='w',format='NETCDF4_CLASSIC')
#    #Add station dimension/variable. Note that multi-character strings must be handled using stringtochar:
#    stn_dim = ncfile.createDimension('station',len(stations))
#    nchars=len(max(stations,key=len))
#    nchar_dim = ncfile.createDimension('nchar',nchars)
#    str_out = netCDF4.stringtochar(np.array(stations,'S'+str(nchars)))
#    stn_out = ncfile.createVariable('station', 'S1', ('station','nchar'))
#    stn_out[:] = str_out
#    stn_out.long_name = "AERONET station name"
#    stn_out.units = "none"
#
#    #Add lat/lon dimensions/variables:
#    lon_dim = ncfile.createDimension('longitude', len(lons))
#    lat_dim = ncfile.createDimension('latitude', len(lats))
#    lon_out = ncfile.createVariable('longitude', np.float64, ('longitude'))
#    lat_out = ncfile.createVariable('latitude', np.float64, ('latitude'))
#    lon_out[:] = lons
#    lat_out[:] = lats
#    lon_out.long_name="longitude"
#    lat_out.long_name="latitude"
#    lon_out.standard_name="longitude"
#    lat_out.standard_name="latitude"
#    lon_out.units="degrees_east"
#    lat_out.units="degrees_north"
#    #Add job dimension/variable. Note that multi-character strings must be handled using stringtochar:
#    job_dim = ncfile.createDimension('job',len(jobids))
#    nchar_dim = ncfile.createDimension('nchar',5)
#    str_out = netCDF4.stringtochar(np.array(jobids,'S5'))
#    job_out = ncfile.createVariable('job', 'S1', ('job','nchar'))
#    job_out[:,:] = str_out
#    job_out.long_name = "job index of UKCA26AER PPE member"
#    job_out.units = "none"
#    #Add AOD550 variable:
#    aod_out = ncfile.createVariable('aod550', np.float64, ('time','latitude','longitude','job'),fill_value=missing)
#    aod_out[:,:,:,:] = AOD_all
#    aod_out.long_name="Aerosol optical depth at 550nm"
#    aod_out.units = "none"
#    aod_out.stash_codes="summation of m01s02i500, m01s02i501, m01s02i502, m01s02i503, m01s02i504, m01s02i505"
#    aod_out.stash_names="AITKENMODESOLUBLEOPTDEPNONADV, ACCUMMODESOLUBLEOPTDEPNONADV, COARSEMODESOLUBLEOPTDEPNONADV, \
#    AITKENMODEINSOLOPTDEPNONADV, ACCUMMODEINSOLOPTDEPNONADV, COARSEMODEINSOLOPTDEPNONADV"
#    #close netCDF file
#    ncfile.close()