#!/bin/bash
#
#Creates and run ed scripts
# - modify_INITFILEENV_start_dump.ed
#   that changes the name and location of start dump to read in INITFILEENV.
#and
# - modify_CNTLALL.ed
#   that changes model basis time and jobid in CNTLALL.
#
#This script creates ed scripts and executes them.
#
#This replaces the IDL program modify_INITFILEENV_start_dump.pro
#
#Written by Masaru Yoshioka based on modify_INITFILEENV_start_dump.pro
#
#below, set up the run by checking and changing in the part 'SET RUN' and 'CHECK&CHANGE'.
#
# LRE 28thSep2015. I've changed the name of the Dump folder, used to indicate where start file is.
# I considered altering .astart (for 1st NRUN only) to *_00 restart file. 
# This needs to happen for resubmission of NRUN ensembles!!! 
# I put a copy of base_job's completed restart simulation into individual job folders (renamed) for use in the CRUN case.
# Also copy all other (non-restart) files, including history folder, which may/may not be used for CRUN.
#
# Next step is ensemble_archer_mod_MY.sh or "run"
#
#Run: ./mod_INITFILEENV_CNTLALL_ed_LRE.sh


ICOMMANDLINE=0  #never comment out. never copy&paste to command line

#comp=MONSOON
comp=ARCHER

srun=n #do not run jobs with this script #SELECT THIS TO SUBMIT AS ENSEMBLE
#srun=y #run jobs immediately with this script
#srun=a #ask if I want to run each job

#SET RUN>>>>>>>>>>>>>>>>>>>>




## LRE 29thDec2015 - !! PI !!

############################################################################################
## LRE 26thNov2015 - Revised and repeated 8thDec2015 b/c VERBOSE_OFF required.
## UKCA_27_PPE_1850 - 1st wave submission, so isamejobid4std=0, to use .astart from base job.
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=342105344 ##330144355
#xpmt=xmdg; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
##xpmt=xmdh; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
##xpmt=xmdi; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
##xpmt=xmdj; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
##xpmt=xmdk; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
##xpmt=xmdl; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
##xpmt=xmdm; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
##xpmt=xmdn; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
##xpmt=xmdo; ja=( a b c d e f g h i ) ## 9/9

#jobidSTD0=xmcnl ## LRE26thNov2015 - correct.   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=0  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=.astart                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 3 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################


############################################################################################
## LRE 26thFeb2016 - 2nd year 3rd, final, wave.
## LRE 24thFeb2016 - 2nd year, 2nd wave.
## LRE 10thFeb2016 - 2nd year resubmission preparation.
## LRE 4thJan2016 - 4th wave
## LRE 2ndJan2016 - 3rd wave
## LRE 29thDec2015 - 2nd wave
## UKCA_27_PPE_1850, isamejobid4std=1, to use 20060201 restart files in individual folders.
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=342105344 ##330144355
#xpmt=xmdg; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#xpmt=xmdh; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
#xpmt=xmdi; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmdj; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmdk; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
#xpmt=xmdl; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmdm; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmdn; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmdo; ja=( a b c d e f g h i ) ## 9/9
######jobidSTD0=xmcnk   #jobid for start dump #effective only if isamejobid4std=0
#
# LRE 24thFeb2016.
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20071001               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 4, 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
# LRE 24thFeb2016.
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070601               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 06 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 7, 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
# LRE 10thFeb2016.
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070201               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 02 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 11, 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
# LRE 4thJan2016:
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20061001               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
############################################################################################
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20060601               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 6 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 7 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#########################################################################################
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20060201               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 2 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 11 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################






## LRE 29thDec2015 !! PD !!

############################################################################################
## LRE 26thNov2015
## LRE 29thFeb2016 - 2nd year 1st submission.
## UKCA_27_PPE_2006_8 - 1st wave submission, so isamejobid4std=0, to use .astart from base job.
## LRE -2ndDec2015 - updated base job with VERBOSE_OFF. 
## LRE -2ndDec2015 - need to reduce to 4 monthly output b/c 40Tb limit.
## LRE 23rdDec2015 - Forced to repeat because the manual umui folder change didn't work.
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
## OLD: runnum=329145744 ## 336122950  
## LRE2ndDec2015. Updated VERBOSE OFF, but runnum used only for ensemble folders. Base job jobidSTD0 (below) used only for name replacement within ensemble folders.
#runnum=345170551
##xpmt=xmct; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
##xpmt=xmcu; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
##xpmt=xmcv; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
##xpmt=xmcw; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
##xpmt=xmcx; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
##xpmt=xmcy; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
##xpmt=xmcz; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
##xpmt=xmda; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmdb; ja=( a b c d e f g h i ) ## 9/9
##
## LRE 29thFeb2016 - 2nd year, first submission.
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070601               #Use individual restart files.
#resub=" 0, 8, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 6 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
##
## LRE 29thFeb2016 - 2nd year, first submission.
##isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070201               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 2 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 11 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
#jobidSTD0=xmcnk   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=0  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=.astart                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 3 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################

############################################################################################
## 4th Submission (completes a year) - running this script 3rdJan2016
## 3rd submission completed.
## 2nd PPE_2006_2008 submission. Using Feb01 restart files. - LRE 29thDec - Success!!
## UKCA_27_PPE_2006_8 - 1st wave submission, so isamejobid4std=0, to use .astart from base job.
## LRE -2ndDec2015 - updated base job with VERBOSE_OFF. 
## LRE -2ndDec2015 - need to reduce to 4 monthly output b/c 40Tb limit.
## LRE 23rdDec2015 - Forced to repeat because the manual umui folder change didn't work.
#
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
####runnum=329145744 ## 336122950  ## LRE2ndDec2015. Updated VERBOSE OFF, but runnum used only for ensemble folders. Base job jobidSTD0 (below) used only for name replacement within ensemble folders.
#runnum=345170551
#xpmt=xmct; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#
#xpmt=xmcu; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
### LRE 23rdDec2015 - NOTE 'J' CAN'T BE REMOVED HERE!!! 
### LRE 27thDec2015 - but job fails b/c restart file not available!!
### LRE 27thDec2015. CAN skip for 2nd submission? There's no change to parameters, so expect to be fine.
### LRE 29thDec2015 - Removed the check later in this script for the restart file existence.
### so that all jobs can be included here.
#xpmt=xmcv; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmcw; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmcx; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
### LRE 27thDec2015 - Job xmcy.a failed at runtime 17/11/05
### LRE 30thDec2015 - Above comment relates to 20060201 resubmission only. Parameters in other files unchanged when 'a' removed from list. Checked.
#xpmt=xmcy; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmcz; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmda; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmdb; ja=( a b c d e f g h i ) # 9/9
#jobidSTD0=xmcnk   #jobid for start dump #effective only if isamejobid4std=0
# LRE 3rdJan2016:
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20061001               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
#
############################################################################################



## LRE 29thDec2015 - !! 1978 !!

############################################################################################
## LRE 11thJan2016 - First submission (.astart)
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=011095052
#xpmt=xmja; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#xpmt=xmjb; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
#xpmt=xmjc; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmjd; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmje; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
#xpmt=xmjf; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmjg; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmjh; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmji; ja=( a b c d e f g h i ) ## 9/9

#jobidSTD0=xmcnm ## LRE26thNov2015 - correct.   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=0  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=.astart                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 3, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 3 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################

############################################################################################
## LRE 20thJan2016 - resubmission 2.
## LRE 22ndJan2016 - resubmission 3.
## LRE 25thJan2016 - resubmission 4.
## LRE 28thJan2016 - 2nd YEAR.
## LRE 4thFeb2016  - 2nd YEAR, resubmission 2.
## LRE 10thFeb2016 - Resubmitting December b/c space crashed ensemble.
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=011095052
#xpmt=xmja; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#xpmt=xmjb; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
#xpmt=xmjc; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmjd; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmje; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
#xpmt=xmjf; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmjg; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmjh; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmji; ja=( a b c d e f g h i ) ## 9/9

###jobidSTD0=xmcnm ## LRE20thJan2016. Now using restarts within folders.   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20071201                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 1, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 12 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 1 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################




############################################################################################
## LRE 16thJun2016 - !! 1998 !!
############################################################################################

## LRE 11thJan2016 - First submission (.astart)
## LRE 21stJan2016 2nd submission
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=165100427
#xpmt=xmtn; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#xpmt=xmto; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
#xpmt=xmtp; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmtq; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmtr; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
#xpmt=xmts; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmtt; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmtu; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmtv; ja=( a b c d e f g h i ) ## 9/9


#jobidSTD0=xmcnm ## jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20060701  ##.astart             #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 6, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 7 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 6 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
##resub=" 0, 5, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
##bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
##tgend=" 3 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)

############################################################################################

############################################################################################


### LRE 21stOct2016 - Starting 2nd year
### LRE 17thNov2016 - started 2nd year again b/c of space issues
### LRE 22ndNov2016 - again restarting b/c of space. 10Tb added.
### causing jobs to stop at multiple places with half-made files
### Also /karthee filepath deleted for ancillary files.

### LRE 28thNov2016 - finally a successful batch of runs. Rinse/repeat
### LRE 30thNov2016 - last batch on last day available!

#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=165100427
##### TEST: xpmt=xmtn; ja=( a )
#xpmt=xmtn; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#xpmt=xmto; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 2/9 LRE 26thNov2015
#xpmt=xmtp; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 3/9
#xpmt=xmtq; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 4/9
#xpmt=xmtr; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 5/9
#xpmt=xmts; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 6/9
#xpmt=xmtt; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 7/9
#xpmt=xmtu; ja=( a b c d e f g h i j k l m n o p q r s t u v w x y z ) ## 8/9
#xpmt=xmtv; ja=( a b c d e f g h i ) ## 9/9

###jobidSTD0=xmcnn ## LRE20thJan2016. Now using restarts within folders.   #jobid for start dump #effective only if isamejobid4std=0

#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070901                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 9 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 5 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)

#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20070101                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2007 , 1 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 1 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)

############################################################################################











#
# LRE 3rdJan2016 - ARCHIVE:
#
#jobidSTD0=xmcnk   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20060601               #Use individual restart files.
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 6 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 2 , 7 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################









## LRE 29thDec2015 - CRUN ARCHIVE:

############################################################################################
## LRE 28thSep2015. CRUN OAT PPE !!!!!!!!!   Will .astart work for CRUN??  (NO)    !!!!!!!!!

#-----------COPY A BLOCK AND MODFY IT----------------
## LRE 29thSept2015 - THIS JOB DIDN'T WORK B/C OF THE .astart when using CRUN
###iextra_check=0
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=271085348
#xpmt=xlwo; ja=( b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
#jobidSTD0=xlwoa   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=0  #0=jobid for start dump set here; 1=same as current jobid
### LRE This means if using .astart, set this to 1. 0 if using restart
#yyyymmdd=.astart                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 3, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)

############################################################################################
## LRE 28thSep2015. CRUN OAT PPE  - COPY OF ABOVE, now using RESTART
## What tells the um code that this is a CRUN job?
#-----------COPY A BLOCK AND MODFY IT----------------
###iextra_check=0
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=271085348
#xpmt=xlwo; ja=( b c d e f g h i j k l m n o p q r s t u v w x y z ) #tgt xpmt & job IDs
##jobidSTD0=xlwoa   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
### LRE This means set to 1 if using .astart, set this 0 if using restart.
#yyyymmdd=20060101               #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 3, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 1 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)

############################################################################################
## LRE 28thSep2015. 2 member CRUN test. Have run CRUN to completion.
## and copied numerous folder and files into work.
#-----------COPY A BLOCK AND MODFY IT----------------
###iextra_check=0
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=288141815
#xpmt=xlwv; ja=( b c ) #tgt xpmt & job IDs
##jobidSTD0=xlwva   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#### LRE This means set to 1 if using .astart, set this 0 if using restart.
#yyyymmdd=20060701               #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2006 , 7 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 11 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################

############################################################################################
## LRE 30thSep2015
## NRUN OAT PPE - returned to perform 3 extra OAT jobs with appropriate mparwtr values.
##iextra_check=0
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=267174400
#xpmt=xlwk; ja=( e f g ) #tgt xpmt & job IDs
#jobidSTD0=xlwga   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=0  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=.astart                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 4, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 1 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 1 , 3 , 1 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################

############################################################################################
## LRE 29thOct2015
## CRUN test 2day regular Q w script changes
##iextra_check=0
#ixarr=0; irarr=0 #experiment id & run number NOT given as an array
#runnum=302094120
#xpmt=xlwv; ja=( n o ) #tgt xpmt & job IDs
#jobidSTD0=xlwva   #jobid for start dump #effective only if isamejobid4std=0
#isamejobid4std=1  #0=jobid for start dump set here; 1=same as current jobid
#yyyymmdd=20051003                     #Use *.astart file!! CAREFUL 4 START DATE!!
#resub=" 0, 2, 0, 0, 0, 0,"          #set RUN_RESUBMIT_INC = resubmission pattern
#bstim=" 2005 , 10 , 3 , 0 , 0 , 0 ," #set MODEL_BASIS_TIME #HAS 2 MATCH yyyymmdd above
#tgend=" 0 , 0 , 8 , 0 , 0 , 0 , "    #set RUN_TARGET_END: y,m,d,0,0,0 (run length)
############################################################################################

#-----------COPY A BLOCK AND MODFY IT----------------

#CAUTION: You may be tempted to use the very last dump available
# (eg, jobida.da20080721_00 stored in /work/n02/n02/$USER/um/$jobid/)
# instead of that at the first day of the last month
# (eg, jobida.da20080701_00 stored in /work/n02/n02/$USER/Dumps/$jobid/).
# The model will run OK with that, but bear in mind that THE MONTHLY OUTPUT
# FOR THAT MONTH (eg, jobida.pm2008jul) WILL NOT BE CREATED. Monthly outputs
# will be created for the second month and later.

#SET RUN<<<<<<<<<<<<<<<<<<<<

#./mod_INITFILEENV_CNTLALL_ed.sh

#SET DIRECTORIES <=============================================
if [ $comp == ARCHER ]; then
  dirHome=/home/n02/n02/$USER
  dirWork=/work/n02/n02/$USER
 #dirDump=/nerc/n02/n02/$USER/Dumps #nerc dir not accesible at run time
  dirum=/work/n02/n02/$USER/um      # so put dump file in work dir
  dirDump=/work/n02/n02/$USER/Dumps  ## LRE 27thDec2015. This is only used when yyyymmdd NOT_EQ .astart.
  ## Only completed month restarts are stored in Dumps, so this is the safer location to use for restarting jobs.
elif [ $comp == MONSOON ]; then
  dirHome=/home/$USER
  dirWork=/projects/ukca/$USER
  dirDump=/nerc/ukca/$USER
else
  echo "comp = "$comp" is not correct."
  if [ "$ICOMMANDLINE" != "1" ]; then exit; fi
fi

dirRUN0=$dirHome/umui_runs
fnmINI=INITFILEENV
fnmINI0=INITFILEENV_ini
fnmCTL=CNTLALL
fnmCTL0=CNTLALL_ini

dirSCR=$dirHome/umui_runs/scripts
dirED=$dirHome/umui_runs/scripts/ed
fnmINIed=modify_INITFILEENV_start_dump.ed
fnmCTLed=modify_CNTLALL.ed

fnmSUB=umuisubmit_run

dirSTD0=$dirDump

ans1=undef
ans2=undef
i=0
for j in ${ja[@]}; do
  echo "**************************************"

  if [ $ixarr == 1 ]; then xpmt=${xpmta[$i]}; fi
  echo xpmt: $xpmt

  jobid=$xpmt${ja[$i]}
  echo jobid: $jobid

  if [ $irarr == 1 ]; then runnum=${runnuma[$i]}; fi
  runid=$jobid-$runnum
  echo runid: $runid

  dirRUN=$dirRUN0/$runid ; ls -ld $dirRUN

##INITFILENV
  ls -l $dirRUN/$fnmINI
  ls -l $dirRUN/$fnmINI0

  if [ ! -e $dirRUN/$fnmINI0 ]; then
    echo "File "$dirRUN/$fnmINI0" does not exist--to be created."
    cp $dirRUN/$fnmINI $dirRUN/$fnmINI0
    ls -l $dirRUN/$fnmINI0
  fi
  echo

  if [ $isamejobid4std == 1 ]; then
    jobidSTD=$jobid
  else
    jobidSTD=$jobidSTD0
  fi

  if [ "$yyyymmdd" == ".astart" ]; then
    dirSTD=$dirum/$jobidSTD
    fnmSTD=${jobidSTD}${yyyymmdd}
  else
    dirSTD=$dirSTD0/$jobidSTD
    fnmSTD=${jobidSTD}a.da${yyyymmdd}_00
  fi

## Deleting this check. I check for all .astart/restart files manually and delete some jobs from submission,
## but for ease want to include all jobs here. Avoids repeated running of this script.

#  ls -l $dirSTD/$fnmSTD
#  if [ ! -e $dirSTD/$fnmSTD ]; then
#    echo Specified start dump file does not exist.
#    if [ "$icmdline" != "1" ]; then exit; fi
#  fi

 #CREATE ED SCRIPT
  echo
  echo Creating: $dirED/$fnmINIed

  cat << endfile > $dirED/$fnmINIed
ed $dirRUN/$fnmINI <<\EOF
/export ASTART=/
d
i
export ASTART=$dirSTD/$fnmSTD
.
w
q
EOF
endfile
  echo

 #CHECK ED SCRIPT
  ls -l $dirED/$fnmINIed
  echo =======================================
  cat $dirED/$fnmINIed
  echo =======================================

  if [ "$ans1" != "yes4all" ]; then
    echo "Ed script is created as above. OK to execute it?"
    echo "  <Type y for yes, yes4all for yes for all jobs, or anything else to exit>"
    read -a ans1
  fi
  if [ "$ans1" != "y" ] && [ "$ans1" != "yes4all" ]; then
    if [ "$icmdline" != "1" ]; then exit; fi
  fi

  chmod +x $dirED/$fnmINIed
  echo Doing: $dirED/$fnmINIed
  $dirED/$fnmINIed

  echo
  echo diff $dirRUN/$fnmINI $dirRUN/$fnmINI0
  echo =======================================
  diff $dirRUN/$fnmINI $dirRUN/$fnmINI0
  echo =======================================
  echo
  sleep 2

##CNTLALL
  ls -l $dirRUN/$fnmCTL
  ls -l $dirRUN/$fnmCTL0

  if [ ! -e $dirRUN/$fnmCTL0 ]; then
    echo "File "$dirRUN/$fnmCTL0" does not exist--to be created."
    cp $dirRUN/$fnmCTL $dirRUN/$fnmCTL0
    ls -l $dirRUN/$fnmCTL0
  fi

# if [ imodjobid == 1 ]; then
#   xid=$xid0; jid=$jid0
# else
    xid=$xpmt; jid=$j
# fi

  echo Creating: $dirED/$fnmCTLed

#  cat << endfile > $dirED/$fnmCTLed
#ed $dirRUN/$fnmCTL <<\EOF
#/ EXPT_ID=/
#d
#i
# EXPT_ID='$xid',
#.
#/ JOB_ID=/
#d
#i
# JOB_ID='$jid',
#.
#/ MODEL_BASIS_TIME=/
#d
#i
# MODEL_BASIS_TIME=$bstim
#.
#w
#q
#EOF
#endfile

  cat << endfile > $dirED/$fnmCTLed
ed $dirRUN/$fnmCTL <<\EOF
/ RUN_RESUBMIT_INC=/
d
i
 RUN_RESUBMIT_INC=$resub
.
/ MODEL_BASIS_TIME=/
d
i
 MODEL_BASIS_TIME=$bstim
.
/ RUN_TARGET_END=/
d
i
 RUN_TARGET_END=$tgend
.
w
q
EOF
endfile

##CHECK ED SCRIPT
  ls -l $dirED/$fnmCTLed
  echo =======================================
  cat $dirED/$fnmCTLed
  echo =======================================

  if [ "$ans2" != "yes4all" ]; then
    echo "Ed script is created as above. OK to execute it?"
    echo "  <Type y for yes, yes4all for yes for all jobs, or anything else to exit>"
    read -a ans2
  fi
  if [ "$ans2" != "y" ] && [ "$ans2" != "yes4all" ]; then
    if [ "$ICOMMANDLINE" != "1" ]; then exit; fi
  fi

##EXECUTE ED SCRIPT
  chmod +x $dirED/$fnmCTLed
  echo Doing: $dirED/$fnmCTLed
  $dirED/$fnmCTLed

##CHECK THE CHANGES MADE
  echo
  echo diff $dirRUN/$fnmCTL $dirRUN/$fnmCTL0
  echo =======================================
  diff $dirRUN/$fnmCTL $dirRUN/$fnmCTL0
  echo =======================================
  echo
  sleep 2

#####SUBMIT THE JOB TO RUN!!!!
  ansrun=n
  if [ $srun == a ]; then

    if [ $comp == MONSOON ]; then
      echo -n "Submit job $jobid to run? This has to be done on IBM."
    fi
    if [ $comp == ARCHER ]; then
      echo -n "Submit job $jobid to run?"
    fi
    echo "  <Type y for yes, yes4all for yes for all jobs, e to exit, or anything else for no>"
    read -a ansrun
  fi
  if [ $ansrun == yes4all ]; then
    srun=y
  fi
  if [ $srun == y ] || [ $ansrun == y ]; then
    if [ $comp == MONSOON ]; then
      echo llsubmit $dirRUN/$fnmSUB
      llsubmit $dirRUN/$fnmSUB
      echo llq -u $USER
      llq -u $USER
    fi
    if [ $comp == ARCHER ]; then
      echo qsub $dirRUN/$fnmSUB
      qsub $dirRUN/$fnmSUB
      echo qstat -u $USER
      qstat -u $USER
    fi
    sleep 1
  fi
  if [ $ansrun == e ]; then
    echo "Exiting the script."
    if [ $ICOMMANDLINE == 0 ]; then exit; fi
  fi
  echo

  i=$(( i+1 ))
done



