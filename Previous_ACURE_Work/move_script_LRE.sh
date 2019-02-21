#!/bin/bash

#Modified by Masaru Yoshioka as marked with 'MY'.

# copies CNTLATM from jobid-####.x1 etc to tdxxa-#### etc.
# CNTLATM contains parameters and their values.
#
## LRE 24thSep2015. I'm adding to this script so that
## INITFILEENV and SHARED are also moved from 
## jobid-####.x1 to jobid-####
#
# set up the parts 'GIVE VALUES HERE' and 'CHECK&CHANGE' below.
#
# this file (move_script_MY.sh) is stored in umui_runs/scripts/
# but it is (and should be) linked to umui_runs/move_script.sh .
#
# To run this script, go to umui_runs/ and type;
# ./scripts/move_script_LRE.sh


ICOMMANDLINE=0

#############COPY A BLOCK & MODIFY VALUES#############
#TRUNK_EXPERIMENT_ID=xlvva ## LRE 24thSep2015. Note the 5 letter code. ## LRE 24thSep2015 - This definition is replaced below!!!
#TRUNK_EXPERIMENT_NUMBER=265084941

#RUN_EXPERIMENT_ID_ONE=tduo #1st set from  x1 to x26 (or 25 if starting xxxxb)
#RUN_EXPERIMENT_ID_ONE=xlvv #1st set from  x1 to x26 (or 25 if starting xxxxb)
#RUN_EXPERIMENT_ID_TWO=xlvw #2nd set from x27 to x48 #dummy if #ensemble=<26

##NUMBER_OF_RUNS=2 ## LRE 25thSep2015 - this is no longer used. 
#############COPY A BLOCK & MODIFY VALUES#############
##TRUNK_EXPERIMENT_ID=tdugt
##TRUNK_EXPERIMENT_NUMBER=018170001

##RUN_EXPERIMENT_ID_ONE=tdvj #1st set from  x1 to x26 (or 25 if starting xxxxb)
##RUN_EXPERIMENT_ID_TWO=tdvk #2nd set from x27 to x48 #dummy if #ensemble=<26

#
##NUMBER_OF_RUNS=31
#NUMBER_OF_RUNS=1
#############COPY A BLOCK & MODIFY VALUES#############
##TRUNK_EXPERIMENT_ID=tdugq
##TRUNK_EXPERIMENT_NUMBER=017224319

#RUN_EXPERIMENT_ID_ONE=tdvh #1st set from  x1 to x26 (or 25 if starting xxxxb)
##RUN_EXPERIMENT_ID_ONE=tdvi #1st set from  x1 to x26 (or 25 if starting xxxxb)
##RUN_EXPERIMENT_ID_TWO=tdvi #2nd set from x27 to x48 #dummy if #ensemble=<26

#NUMBER_OF_RUNS=31
##NUMBER_OF_RUNS=1
#############COPY A BLOCK & MODIFY VALUES#############

######CAUTION: THE WAY TO GIVE VALUES CHANGED BELOW THIS LINE###############################

#############COPY A BLOCK & MODIFY VALUES#############
## LRE 25thSep2015 - CRUN FULL OAT
#TRUNK_RUNID=xlwoa-271085348
##
#RUN_EXPERIMENT_ID_ONE=xlwo #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( b c d e f g h i j k l m n o p q r s t u v w x y z )
#iset2=1                    #=1 2nd set necessary, =0 not necessary
#RUN_EXPERIMENT_ID_TWO=xlwn #2nd set #dummy if iset2=0
#ja2=( a b c )          #dummy if iset2=0
######################################################

#############COPY A BLOCK & MODIFY VALUES#############
## LRE 28thSep2015 - FULL OAT NRUN.
#TRUNK_RUNID=xlwga-267174400
##
#RUN_EXPERIMENT_ID_ONE=xlwg #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( b c d e f g h i j k l m n o p q r s t u v w x y z )
#iset2=1                    #=1 2nd set necessary, =0 not necessary
#RUN_EXPERIMENT_ID_TWO=xlwk #2nd set #dummy if iset2=0
#ja2=( a b c )          #dummy if iset2=0
######################################################

#############COPY A BLOCK & MODIFY VALUES#############
## LRE 28thSep2015 - FULL OAT NRUN - LRE 12thOct2015. Repeated call of base job
#TRUNK_RUNID=xlwga-267174400
##
#RUN_EXPERIMENT_ID_ONE=xlwk #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( d e )
#iset2=0                    #=1 2nd set necessary, =0 not necessary
##RUN_EXPERIMENT_ID_TWO=xlwk #2nd set #dummy if iset2=0
##ja2=( a b c )          #dummy if iset2=0
######################################################


#############COPY A BLOCK & MODIFY VALUES#############
## LRE 25thSep2015 - 2 day CRUN test
#TRUNK_RUNID=xlvza-267173247
#
#RUN_EXPERIMENT_ID_ONE=xlvz #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( b c d e f g h i j k l m n o p q r s t u v w x y z )
#iset2=1                    #=1 2nd set necessary, =0 not necessary
#RUN_EXPERIMENT_ID_TWO=xlwa #2nd set #dummy if iset2=0
#ja2=( a b c )          #dummy if iset2=0
#############COPY A BLOCK & MODIFY VALUES#############
## LRE 25thSep2015 - 2 day NRUN test
#TRUNK_RUNID=xlvva-265084941
#
#RUN_EXPERIMENT_ID_ONE=xlvv #1st set from  x1 to x26 (or 25 if starting from xxxxb)
##ja1=( b c d e f g h i j k l m n o p q r s t u v w x y z )
#ja1=( b c d e f g h i j k l m n o p q r s t u v w x y z )
#iset2=1                    #=1 2nd set necessary, =0 not necessary
#RUN_EXPERIMENT_ID_TWO=xlvw #2nd set #dummy if iset2=0
#ja2=( a b c )          #dummy if iset2=0
#############COPY A BLOCK & MODIFY VALUES#############

## Masaru's:

##TRUNK_RUNID=tdwnq-100140808

##RUN_EXPERIMENT_ID_ONE=tdwn #1st set from  x1 to x26 (or 25 if starting from xxxxb)
##ja1=( r s t )

##iset2=0                    #=1 2nd set necessary, =0 not necessary
##RUN_EXPERIMENT_ID_TWO=tdwn #2nd set #dummy if iset2=0
##ja2=( a b c d e )          #dummy if iset2=0
#############COPY A BLOCK & MODIFY VALUES#############
##TRUNK_RUNID=tdyod-226152252

##RUN_EXPERIMENT_ID_ONE=tdyo #1st set from  x1 to x26 (or 25 if starting from xxxxb)
##ja1=( e f )

##iset2=0                    #=1 2nd set necessary, =0 not necessary
##RUN_EXPERIMENT_ID_TWO=tdxl #dummy if iset2=0
##ja2=( a b c d e )          #dummy if iset2=0
#############COPY A BLOCK & MODIFY VALUES#############
##TRUNK_RUNID=tdyog-231110756

##RUN_EXPERIMENT_ID_ONE=tdyr #1st set from  x1 to x26 (or 25 if starting from xxxxb)
##ja1=( a b c d e f g h i j k l m n o p q r s t u v w x y z )

##iset2=0                    #=1 2nd set necessary, =0 not necessary
##RUN_EXPERIMENT_ID_TWO=tdxl #dummy if iset2=0
##ja2=( a b c d e )          #dummy if iset2=0


#############COPY A BLOCK & MODIFY VALUES#############
## LRE 19thOct2015 - 2 job test of CRUN
#TRUNK_RUNID=xlwva-288141815
##
#RUN_EXPERIMENT_ID_ONE=xlwv #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( b c )
#iset2=0                    #=1 2nd set necessary, =0 not necessary
#RUN_EXPERIMENT_ID_TWO=xlwk #2nd set #dummy if iset2=0
#ja2=( a b c )          #dummy if iset2=0
######################################################


#############COPY A BLOCK & MODIFY VALUES#############
## LRE 20thOct2015. Returned to OAT job to implement 3 tests with mparwtr
#TRUNK_RUNID=xlwga-267174400
##
#RUN_EXPERIMENT_ID_ONE=xlwk #1st set from  x1 to x26 (or 25 if starting from xxxxb)
#ja1=( e f g )
#iset2=0                    #=1 2nd set necessary, =0 not necessary
###RUN_EXPERIMENT_ID_TWO=xlwk #2nd set #dummy if iset2=0
###ja2=( a b c )          #dummy if iset2=0
######################################################

#############COPY A BLOCK & MODIFY VALUES#############
## LRE 20thOct2015. CRUN 2day regular Q test w script changes
TRUNK_RUNID=xlwvm-302094120
#
RUN_EXPERIMENT_ID_ONE=xlwv #1st set from  x1 to x26 (or 25 if starting from xxxxb)
ja1=( n o )
iset2=0                    #=1 2nd set necessary, =0 not necessary
##RUN_EXPERIMENT_ID_TWO=xlwk #2nd set #dummy if iset2=0
##ja2=( a b c )          #dummy if iset2=0
######################################################





# To run this script, go to umui_runs/ and type;
# move_script_LRE.sh

strarray=(${TRUNK_RUNID//-/ })
#${string//a/b} replaces a with b. as b=' ', this forms (str1 str2) which defines an array!!!

TRUNK_EXPERIMENT_ID=${strarray[0]}
TRUNK_EXPERIMENT_NUMBER=${strarray[1]}

ARCHIVE_NAMELIST_FOLDER=$TRUNK_EXPERIMENT_ID_store_namelists

######## LOOP OVER THE FIRST SET OF JOB IDs

count=0

#MY: FIRST SET ################################################################
for letter in ${ja1[@]}; do
  count=$((count+1))

  echo $count::$letter

  ls -l $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM
  ls -ld $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/

  echo cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/.
  cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/.

  echo diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/CNTLATM $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/CNTLATM
  diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/CNTLATM $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/CNTLATM

  ## LRE 25thSep2015. Changes to INITFILEENV and SHARED included here:

  echo diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV
  diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV

  ls -l $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED
  ls -ld $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/

  echo cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/.
  cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/.

  echo diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/SHARED $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/SHARED
  diff $RUN_EXPERIMENT_ID_ONE$letter-$TRUNK_EXPERIMENT_NUMBER/SHARED $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/SHARED

  #MY: echo $RUN_EXPERIMENT_ID_ONE${letter}a.dak7c10

  echo count = $count
  echo
done
#MY: FIRST SET ################################################################

#MY: SECOND SET ###############################################################
if [ $iset2 != 1 ]; then
  echo move_script finished.
  if [ $ICOMMANDLINE == 0 ]; then exit; fi
fi

for letter in ${ja2[@]}; do

  count=$((count+1))

  echo $count::$letter

  ls -l $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM
  ls -ld $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/

  echo cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/.
  cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/CNTLATM $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/.

  echo diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/CNTLATM $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/CNTLATM
  diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/CNTLATM $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/CNTLATM

  ## LRE 25thSep2015. INITFILEENV and SHARED changes:

  echo diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV
  diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/INITFILEENV

  ls -l $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED
  ls -ld $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/

  echo cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/.
  cp -f $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER.x$count/SHARED $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/.

  echo diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/SHARED $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/SHARED
  diff $RUN_EXPERIMENT_ID_TWO$letter-$TRUNK_EXPERIMENT_NUMBER/SHARED $TRUNK_EXPERIMENT_ID-$TRUNK_EXPERIMENT_NUMBER/SHARED

  echo count = $count
  echo
done
#MY: SECOND SET ###############################################################

####MY: COPY THE BLOCK ABOVE AS MANY TIMES AS NECESSARY TO ACCOMMODATE ALL JOBS
####MY: ALL YOU WILL NEED TO CHANGE FOR 3RD SET IS '2'. THAT NEEDS TO BE CHANGED INTO '3' FOR EXAMPLE.

echo move_script finished.
