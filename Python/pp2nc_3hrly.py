import iris
import numpy as np
import string
import datetime as dt
import netCDF4
import os
import argparse

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
#ppRoot = '/home/users/myoshioka/ET/Retrieved/group_workspaces/jasmin2/ukca/vol1/myoshioka/um/Dumps'
#ncRoot = '/group_workspaces/jasmin2/gassp/earjjo/NC_3hrly'
#ncRef = '/group_workspaces/jasmin2/ukca/vol1/myoshioka/um/Dumps/Links2NC/teafw/aod550_total_teafw_pm2008jan_N48_dim_att_j.nc'
#startDate = '20080101'
#endDate = '20080101'
#####

#####CHECK FILES/DIRS exist
assert os.path.exists(ppRoot), "ppRoot does not exist"
assert os.path.exists(ncRoot), "ncRoot does not exist"
assert os.path.exists(ncRef),"ncRef does not exist"
#####

#####PARAMETERS
stash_CDNC = 'm01s38i479' #STASH code for CDNC field in pp files
stash_AOD = ['m01s02i500','m01s02i501','m01s02i502','m01s02i503','m01s02i504','m01s02i505'] #STASH codes for 6 'modes' in pp files. AOD_550 is sum over all these modes
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
    startYr=startDate[0:4]
    startMon=startDate[4:6]
    startDay=startDate[6:8]
    startDateObj=dt.date(int(startYr),int(startMon),int(startDay))
except:
    print("Error: startDate input parameter in unexpected format")
    raise
try:
    endYr=endDate[0:4]
    endMon=endDate[4:6]
    endDay=endDate[6:8]
    endDateObj=dt.date(int(endYr),int(endMon),int(endDay))
except:
    print("Error: endDate input parameter in unexpected format")
    raise
numDays=(endDateObj-startDateObj).days + 1 #Add one to include both extremeties
assert numDays > 0, "Error: startDate must come before endDate"
assert numDays < 32, "Error: Can only process up to one month of data at a time"
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
assert len(lats) == 73 and len(lons) == 96, "Unexpected dimensions in reference data file"
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
    missing=-999
    CDNC_all=np.full(shape=(8,73,96,235),fill_value=missing)
    AOD_all=np.full(shape=(8,73,96,235),fill_value=missing)
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
        #Try reading in CDNC cube and check dimensions are as expected (assume AOD cubes are OK by association):
        try:
            CDNCcube = iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stash_CDNC))
        except:
            flog.write('Warning: Could not load CDNC data from input file '+ppFile+'. File has been skipped\n')
            flog.close
            continue
        if CDNCcube.shape != (8,52,145,192):
            flog.write('Warning: Unexpected dimensions in input file '+ppFile+'. File has been skipped\n')
            flog.close
            continue
        nFilesRead+=1
        #If first file read in, store time information:
        if nFilesRead==1:
            try:
                times = CDNCcube.coord('time').points
            except:
                print("Error: No time coordinate in pp file")
                raise
        #Read in all 6 cubes relating to AOD and sum to get total AOD:
        try:
            for s,stashid in enumerate(stash_AOD):
                modeCube=iris.load_cube(ppPath,iris.AttributeConstraint(STASH=stashid))
                if s==0:
                    AODcube=modeCube
                else:
                    AODcube+=modeCube
        except:
            flog.write('Warning: Could not load AOD data from input file '+ppFile+'. File has been skipped\n')
            flog.close
            continue
        #***FOR NOW*** just take arithmetic mean of CDNC over all heights
        for h in range(52):
            if h==0:
                CDNCsum=CDNCcube[:,h,:,:]
            else:
                CDNCsum+=CDNCcube[:,h,:,:]
        CDNCave=CDNCsum/52.0
        #Regrid onto coarser resolution (N48):
        CNDC_N48 = CDNCave.regrid(refCube, iris.analysis.Linear())
        AOD_N48 = AODcube.regrid(refCube, iris.analysis.Linear())
        #Update big arrays (all jobs for current day)
        CDNC_all[:,:,:,j]=CNDC_N48.data
        AOD_all[:,:,:,j]=AOD_N48.data
        flog.close

    #Create nc file:
    assert not (np.max(CDNC_all) == missing and np.min(CDNC_all) == missing), "No processed data"
    ncfile = netCDF4.Dataset(outFile,mode='w',format='NETCDF4_CLASSIC')
    #Add time dimension/variable:
    t_dim = ncfile.createDimension('time',len(times))
    t_out = ncfile.createVariable('time', np.float64, ('time'))
    t_out[:] = times
    t_out.units = "hours since 1970-01-01 00:00:00"
    t_out.calendar = "gregorian"
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
    aod_out = ncfile.createVariable('AOD550_total', np.float64, ('time','latitude','longitude','job'))
    aod_out[:,:,:,:] = AOD_all
    aod_out.units = "1"
    #Add CDNC variable:
    cdnc_out = ncfile.createVariable('CDNC', np.float64, ('time','latitude','longitude','job'))
    cdnc_out[:,:,:,:] = CDNC_all
    cdnc_out.units = "number per metre cubed"
    #close netCDF file
    ncfile.close()
