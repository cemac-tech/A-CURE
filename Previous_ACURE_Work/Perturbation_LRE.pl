#!/usr/bin/perl -w
#
# This Script Generates Perturbations To A UM Job
#
# Written By: Steven Pickering
# Last Modified: 4/2/2013
#
# Modified by Masaru Yoshioka. Modifications are marked with 'MY'.
#
## LRE 23rdSep2015. All MY alterations previously removed by MY.
## I'm adding Namelist and Control_file options for atmospheric parameters.
#
# As the first step of creating a new set of PPE jobs
# - Modify Variables.list as necessary,
# - Modify Parameters.list as necessary,
# - Run this script: modifying this file is not usually necessary. 
#
# Run: ./Perturbation_LRE.pl

#----------------#
# Use Statements #
#----------------#

use strict;
use warnings;
use File::Copy;

#------------------#
# Global Variables #
#------------------#

# Control Files

my @Valid_Control_File_Names = ("CNTLALL", "CNTLATM", "CNTLGEN", "SHARED");
my @Valid_All_Control_Namelist_Names = ("NLSTCALL", "NLSTWRITEDATA");
my @Valid_Atmosphere_Control_Namelist_Names = ("NLSTCATM", "RUN_TFILT", "RUN_BL", "RUN_BLICE", "RUN_BLVEG", "RUN_PFT", "RUN_LAND", "RUN_Precip", "RUN_Stochastic", "RUN_Cloud", "RUN_Convection", "RUN_Radiation", "RUN_GWD", "RUN_GWD", "RUN_UKCA", "RUN_Dyn", "RUN_SL", "RUN_Diffusion", "RUN_RIVERS", "RADFCDIA", "PPRINTXN", "ANCILCTA", "INTFCNSTA", "R2SWNCAL", "R2SWCLNL", "R2LWNCAL", "R2LWCLNL", "CLMCHFCG");
my @Valid_General_Control_Namelist_Names = ("NLSTCGEN", "NLST_MPP", "BOUNCNST");

# Variables

my $Running_On_MONSOON = 0;
my $Number_Of_Variables = 0;
my $Number_Of_Perturbations = 0;
my $Host_Name = "";
my %Variables_Offset = ();
my @Permutation_Tokens;

# UM Job

my $UM_Directory_Name = "";
my $UM_Work_Directory_Name = "";
my $Perturbation_Directory_Name = "";
my $Perturbation_Work_Directory_Name = "";
my $UM_Job_Name = "";

#-------------#
# Subroutines #
#-------------#

#---------------------#
# Ask For UM Job Name #
#---------------------#

sub Ask_For_UM_Job_Name()

{

# Asks For The Name Of The UM Job

print("\n");
print("Enter The Name Of The UM Job.   E.g. xhuna-313092737\n");
print("\n");

# The User Enters The Value

$UM_Job_Name = <STDIN>;
chomp($UM_Job_Name);

# Display A Blank Line

print("\n");

# Store The Name Of The Host

$Host_Name = system("hostname");

# Determine Which Host You Are Running On

if ($Host_Name =~ m/hector/)

   {

   $Running_On_MONSOON = 0;

   }

# Calculate The Name Of The UM Directory

$UM_Directory_Name = $ENV{'HOME'} . "/umui_runs/" . $UM_Job_Name;

# Calculate The Name Of The UM Work Directory

## LRE 23rdSep2015. Printing the switch for RUNNING_ON_MONSOON

print("Monsoon on/off index " . $Running_On_MONSOON . "\n");


if ($Running_On_MONSOON == 1)

   {

   $UM_Work_Directory_Name = $ENV{'HOME'} . "/data_output/" . substr($UM_Job_Name, 0, index($UM_Job_Name, '-'));

   }

else

   {

## LRE 23rdSep2015. Problem with /bin copying.
## Must be on /work, so include a symbolic link in $home
## called 'work' which points to /work/n02/n02/lre/

##    $UM_Work_Directory_Name = $ENV{'HOME'} . "/um/" . substr($UM_Job_Name, 0, index($UM_Job_Name, '-')); 
## LRE 23rdSep2015 Adding print statement for work directory:
# ##This is just /home/n02/n02/lre/um/xlvva which already exists
#   print("Created UM_Work_Directory " . $UM_Work_Directory_Name . "\n");
#
# ## LRE removed /work  This now fails b/c 'bin' does not exist on /home version. Replaced 23rdSep2015 b/c dir non-existent.
### Simply making the /bin directory to check what happens.
#    $UM_Work_Directory_Name = $ENV{'HOME'} . "/work/" . substr($UM_Job_Name, 0, index($UM_Job_Name, '-'));
#   print("Created alternative UM_Work_Directory " . $UM_Work_Directory_Name . "\n");
#    $UM_Work_Directory_Name = "/work/n02/n02/lre/um/" . substr($UM_Job_Name, 0, index($UM_Job_Name, '-'));
#   print("Created !another! alternative UM_Work_Directory " . $UM_Work_Directory_Name . "\n");
    $UM_Work_Directory_Name = $ENV{'HOME'} . "/work/um/" . substr($UM_Job_Name, 0, index($UM_Job_Name, '-'));
   print("reCreated original UM_Work_Directory " . $UM_Work_Directory_Name . "\n");

   }

}

#----------------#
# Copy Directory #
#----------------#

sub Copy_Directory()

{

# Creates A Copy Of The Specified Directory


my $File_Name = "";
my $Source_Directory_Name = "";
my $Destination_Directory_Name = "";

# Store The Passed Parameter

$Source_Directory_Name = $_[0];
$Destination_Directory_Name = $_[1];

print("Source_Directory_Name = " . $Source_Directory_Name . "\n" );
print("Destination_Directory_Name = " . $Destination_Directory_Name . "\n" );

# The Destination Directory Already Exists

# Open The Directory

opendir(DIR, $Source_Directory_Name) or die $!;

# Loop Over The Files In The Directory

while ($File_Name = readdir(DIR))

      {

# Skip The Files Starting With A Dot

   if (($File_Name !~ m/^\./) && (-f $Source_Directory_Name . "/" . $File_Name))

      {

# Make A Backup Copy Of The File

      &Copy_File($Source_Directory_Name . "/" . $File_Name, $Destination_Directory_Name . "/" . $File_Name);

      }

   }

# Close The Directory

closedir(DIR);

}

#-----------#
# Copy File #
#-----------#

sub Copy_File()

{

# Makes A Copy Of The Specified File

my $Count = -1;
my $Source_File_Name = "";
my $Destination_File_Name = "";

# Store The Passed Parameter

$Source_File_Name = $_[0];
$Destination_File_Name = $_[1];

# Copy The File

copy($Source_File_Name, $Destination_File_Name);

# Change The Permissions On The File

$Count = chmod 0755, $Destination_File_Name;

}

#----------------------------#
# Is Valid Control File Name #
#----------------------------#

sub Is_Valid_Control_File_Name()

{

# Checks Whether The Specified Name Is That Of A Valid Control File

my $Loop = -1;
my $Found = 0;
my $File_Name_1 = "";

# Store The Passed Parameter

$File_Name_1 = $_[0];

# Check To See If It Is A Valid Control File Name

for ($Loop = 0; $Loop < @Valid_Control_File_Names; $Loop++)

   {

# Compare The Names

   if ($File_Name_1 eq $Valid_Control_File_Names[$Loop])

      {

# Show That A Match Was Found

      $Found = 1;

      }

   }

# Return The Result Of The Function

return($Found);

}

#------------------------------------#
# Is Valid All Control Namelist Name #
#------------------------------------#

sub Is_Valid_All_Control_Namelist_Name()

{

# Checks Whether The Specified Name Is That Of A Valid All Control Name List

my $Loop_1 = -1;
my $Found_1 = 0;
my $Name_List = "";

# Store The Passed Parameter

$Name_List = $_[0];

# Check To See If It Is A Valid All Control Name List

for ($Loop_1 = 0; $Loop_1 < @Valid_All_Control_Namelist_Names; $Loop_1++)

   {

# Compare The Names

   if ($Name_List eq $Valid_All_Control_Namelist_Names[$Loop_1])

      {

# Show That A Match Was Found

      $Found_1 = 1;

      }

   }

# Return The Result Of The Function

return($Found_1);

}

#-------------------------------------------#
# Is Valid Atmosphere Control Namelist Name #
#-------------------------------------------#

sub Is_Valid_Atmosphere_Control_Namelist_Name()

{

# Checks Whether The Specified Name Is That Of A Valid Atmosphere Control Name List

my $Loop_2 = -1;
my $Found_2 = 0;
my $Name_List_2 = "";

# Store The Passed Parameter

$Name_List_2 = $_[0];

# Check To See If It Is A Valid Atmosphere Control Name List

for ($Loop_2 = 0; $Loop_2 < @Valid_Atmosphere_Control_Namelist_Names; $Loop_2++)

   {

# Compare The Names

   if ($Name_List_2 eq $Valid_Atmosphere_Control_Namelist_Names[$Loop_2])

      {

# Show That A Match Was Found

      $Found_2 = 1;

      }

   }

# Return The Result Of The Function

return($Found_2);

}

#----------------------------------------#
# Is Valid General Control Namelist Name #
#----------------------------------------#

sub Is_Valid_General_Control_Namelist_Name()

{

# Checks Whether The Specified Name Is That Of A Valid General Control Name List

my $Loop_3 = -1;
my $Found_3 = 0;
my $Name_List_3 = "";

# Store The Passed Parameter

$Name_List_3 = $_[0];

# Check To See If It Is A Valid General Control Name List

for ($Loop_3 = 0; $Loop_3 < @Valid_General_Control_Namelist_Names; $Loop_3++)

   {

# Compare The Names

   if ($Name_List_3 eq $Valid_General_Control_Namelist_Names[$Loop_3])

      {

# Show That A Match Was Found

      $Found_3 = 1;

      }

   }

# Return The Result Of The Function

return($Found_3);

}

#-------------------------#
# Modify All Control File #
#-------------------------#

sub Modify_All_Control_File()

{

# Modifies The All Control File

my $Loop_4 = -1;
my $Index = -1;
my $Permutation_Name = "";
my $First_Character = "";
my $Line = "";
my $Original_Line = "";
my $Key = "";
my $Value = "";
my $Namelist_Name = "";
my $Variable_Name = "";
my $File_Name_4 = "";
my @Lines;
my @Values;
my @Tokens;

# Store The Passed Parameter

$Permutation_Name = $_[0];

# Store The File Name

$File_Name_4 = $Permutation_Name . "/CNTLALL";

# Open The Control File

open(CONTROL_FILE, "<" . $File_Name_4);

# Read In All Of The Lines

@Lines = <CONTROL_FILE>;

# Close The Control File

close(CONTROL_FILE);

# Create The Modified Control File

open(MODIFIED_CONTROL_FILE, ">" . $File_Name_4 . ".modified");

# Process The Lines Of Text

for ($Loop_4 = 0; $Loop_4 < @Lines; $Loop_4++)

   {

# Extract A Line

   $Original_Line = $Lines[$Loop_4];
   chomp($Original_Line);

# Save The Line

   $Line = $Original_Line;

# Remove All Spaces

   $Line =~ s/ //g;

# Extract The First Character

   $First_Character = substr($Line, 0, 1);

# Check If It Is A Name List

   if ($First_Character eq "&")

      {

# Store The Name List Name

      $Namelist_Name = substr($Line, 1);

      }

# Check If It Is A Variable Pair

   if (index($Line, "=") != -1)

      {

# Split The Variable Pair Into Two Tokens

      @Tokens = split(/=/, $Line);

# Extract The Variable Name

      $Variable_Name = $Tokens[0];

# Create The Key

      $Key = "CNTLALL_" . $Namelist_Name . "_" . $Variable_Name;

# Check If The Key Exists

      if (exists($Variables_Offset{$Key}))

         {

# Extract The Variable's Index

         $Index = $Variables_Offset{$Key};

# Extract The Value, For That Variable

         $Value = $Permutation_Tokens[$Index - 1];

# Create The Modified Line

         $Line = $Variable_Name . "=" . $Value . ",";

# Store It Back In The Original Line

         $Original_Line = $Line;

         }

      }

# Write The Line Back To The Modified File

   print(MODIFIED_CONTROL_FILE $Original_Line . "\n");

   }

# Close The Modified Control File

close(MODIFIED_CONTROL_FILE);

# Rename The Original Control File

rename($File_Name_4 . ".modified", $File_Name_4);

}

#--------------------------------#
# Modify Atmosphere Control File #
#--------------------------------#

sub Modify_Atmosphere_Control_File()

{

# Modifies The Atmosphere Control File

my $Loop = -1;
my $Index = -1;
my $Permutation_Name = "";
my $First_Character = "";
my $Line = "";
my $Original_Line = "";
my $Key = "";
my $Value = "";
my $Namelist_Name = "";
my $Variable_Name = "";
my $File_Name = "";
my @Lines;
my @Values;
my @Tokens;

# Store The Passed Parameter

$Permutation_Name = $_[0];

# Store The File Name

$File_Name = $Permutation_Name . "/CNTLATM";

# Open The Control File

open(CONTROL_FILE, "<" . $File_Name);

# Read In All Of The Lines

@Lines = <CONTROL_FILE>;

# Close The Control File

close(CONTROL_FILE);

# Create The Modified Control File

open(MODIFIED_CONTROL_FILE, ">" . $File_Name . ".modified");

# Process The Lines Of Text

for ($Loop = 0; $Loop < @Lines; $Loop++)

   {

# Extract A Line

   $Original_Line = $Lines[$Loop];
   chomp($Original_Line);

# Save The Line

   $Line = $Original_Line;

# Remove All Spaces

   $Line =~ s/ //g;

# Extract The First Character

   $First_Character = substr($Line, 0, 1);

# Check If It Is A Name List

   if ($First_Character eq "&")

      {

# Store The Name List Name

      $Namelist_Name = substr($Line, 1);

      }

# Check If It Is A Variable Pair

   if (index($Line, "=") != -1)

      {

# Split The Variable Pair Into Three Tokens

      @Tokens = split(/=/, $Line);

# Extract The Variable Name

      $Variable_Name = $Tokens[0];

# Create The Key

      $Key = "CNTLATM_" . $Namelist_Name . "_" . $Variable_Name;

# Check If The Key Exists

      if (exists($Variables_Offset{$Key}))

         {

# Extract The Variable's Index

         $Index = $Variables_Offset{$Key};

# Extract The Value, For That Variable

         $Value = $Permutation_Tokens[$Index - 1];

# Create The Modified Line

         $Line = $Variable_Name . "=" . $Value . ",";

# Store It Back In The Original Line

         $Original_Line = $Line;

         }

      }

# Write The Line Back To The Modified File

   print(MODIFIED_CONTROL_FILE $Original_Line . "\n");

   }

# Close The Modified Control File

close(MODIFIED_CONTROL_FILE);

# Rename The Original Control File

rename($File_Name . ".modified", $File_Name);

}


#--------------------------------#
# Modify Shared Control File     #
#--------------------------------#

sub Modify_Shared_Control_File()

{

# Modifies The Shared File (for cloud erosion rate dbsdtbs_turb_0 )

my $Loop = -1;
my $Index = -1;
my $Permutation_Name = "";
my $First_Character = "";
my $Line = "";
my $Original_Line = "";
my $Key = "";
my $Value = "";
my $Namelist_Name = "";
my $Variable_Name = "";
my $File_Name = "";
my @Lines;
my @Values;
my @Tokens;

# Store The Passed Parameter

$Permutation_Name = $_[0];

# Store The File Name

$File_Name = $Permutation_Name . "/SHARED";

# Open The Control File

open(CONTROL_FILE, "<" . $File_Name);

# Read In All Of The Lines

@Lines = <CONTROL_FILE>;

# Close The Control File

close(CONTROL_FILE);

# Create The Modified Control File

open(MODIFIED_CONTROL_FILE, ">" . $File_Name . ".modified");

# Process The Lines Of Text

for ($Loop = 0; $Loop < @Lines; $Loop++)

   {

# Extract A Line

   $Original_Line = $Lines[$Loop];
   chomp($Original_Line);

# Save The Line

   $Line = $Original_Line;

# Remove All Spaces

   $Line =~ s/ //g;

# Extract The First Character

   $First_Character = substr($Line, 0, 1);

# Check If It Is A Name List

   if ($First_Character eq "&")

      {

# Store The Name List Name

      $Namelist_Name = substr($Line, 1);

      }

# Check If It Is A Variable Pair

   if (index($Line, "=") != -1)

      {

# Split The Variable Pair Into Three Tokens

      @Tokens = split(/=/, $Line);

# Extract The Variable Name

      $Variable_Name = $Tokens[0];

# Create The Key

      $Key = "SHARED_" . $Namelist_Name . "_" . $Variable_Name;

# Check If The Key Exists

      if (exists($Variables_Offset{$Key}))

         {

# Extract The Variable's Index

         $Index = $Variables_Offset{$Key};

# Extract The Value, For That Variable

         $Value = $Permutation_Tokens[$Index - 1];

# Create The Modified Line

         $Line = $Variable_Name . "=" . $Value . ",";

# Store It Back In The Original Line

         $Original_Line = $Line;

         }

      }

# Write The Line Back To The Modified File

   print(MODIFIED_CONTROL_FILE $Original_Line . "\n");

   }

# Close The Modified Control File

close(MODIFIED_CONTROL_FILE);

# Rename The Original Control File

rename($File_Name . ".modified", $File_Name);

}



#-----------------------------#
# Modify General Control File #
#-----------------------------#

sub Modify_General_Control_File()

{

# Modifies The General Control File

my $Loop = -1;
my $Index = -1;
my $Permutation_Name = "";
my $First_Character = "";
my $Line = "";
my $Original_Line = "";
my $Key = "";
my $Value = "";
my $Namelist_Name = "";
my $Variable_Name = "";
my $File_Name = "";
my @Lines;
my @Values;
my @Tokens;

# Store The Passed Parameter

$Permutation_Name = $_[0];

# Store The File Name

$File_Name = $Permutation_Name . "/CNTLGEN";

# Open The Control File

open(CONTROL_FILE, "<" . $File_Name);

# Read In All Of The Lines

@Lines = <CONTROL_FILE>;

# Close The Control File

close(CONTROL_FILE);

# Create The Modified Control File

open(MODIFIED_CONTROL_FILE, ">" . $File_Name . ".modified");

# Process The Lines Of Text

for ($Loop = 0; $Loop < @Lines; $Loop++)

   {

# Extract A Line

   $Original_Line = $Lines[$Loop];
   chomp($Original_Line);

# Save The Line

   $Line = $Original_Line;

# Remove All Spaces

   $Line =~ s/ //g;

# Extract The First Character

   $First_Character = substr($Line, 0, 1);

# Check If It Is A Name List

   if ($First_Character eq "&")

      {

# Store The Name List Name

      $Namelist_Name = substr($Line, 1);

      }

# Check If It Is A Variable Pair

   if (index($Line, "=") != -1)

      {

# Split The Variable Pair Into Three Tokens

      @Tokens = split(/=/, $Line);

# Extract The Variable Name

      $Variable_Name = $Tokens[0];

# Create The Key

      $Key = "CNTLGEN_" . $Namelist_Name . "_" . $Variable_Name;

# Check If The Key Exists

      if (exists($Variables_Offset{$Key}))

         {

# Extract The Variable's Index

         $Index = $Variables_Offset{$Key};

# Extract The Value, For That Variable

         $Value = $Permutation_Tokens[$Index - 1];

# Create The Modified Line

         $Line = $Variable_Name . "=" . $Value . ",";

# Store It Back In The Original Line

         $Original_Line = $Line;

         }

      }

# Write The Line Back To The Modified File

   print(MODIFIED_CONTROL_FILE $Original_Line . "\n");

   }

# Close The Modified Control File

close(MODIFIED_CONTROL_FILE);

# Rename The Original Control File

rename($File_Name . ".modified", $File_Name);

}

#--------------------------#
# Modify Submit Check File #
#--------------------------#

sub Modify_Submit_Check_File()

{

# Modifies The Submit Check File

my $Loop = -1;
my $Permutation_Name = "";
my $Line = "";
my $File_Name = "";
my @Lines;

# Store The Passed Parameter

$Permutation_Name = $_[0];

# Store The File Name

$File_Name = $Permutation_Name . "/submitchk";

# Open The Submit File

open(SUBMIT_FILE, "<" . $File_Name);

# Read In All Of The Lines

@Lines = <SUBMIT_FILE>;

# Close The Submit File

close(SUBMIT_FILE);

# Create The Modified Submit File

open(MODIFIED_SUBMIT_FILE, ">" . $File_Name . ".modified");

# Initialise The Line Counter

$Loop = 0;

# Extract The First Line

$Line = $Lines[$Loop];

# Search For "DELJOBDIR"

while ($Line !~ /DELJOBDIR=/)

   {

# Write The Line Out Unchaged

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

   }

# The Line Now Cotains "DELJOBDIR="

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE "DELJOBDIR=false\n");

# Increment The Number Of Lines

$Loop++;

# Process The Remaining Lines

while ($Loop < @Lines)

   {

# Write Out The Remaining Line

   print(MODIFIED_SUBMIT_FILE $Lines[$Loop]);

# Increment The Number Of Lines

   $Loop++;

   }

# Close The Modified Submit File

close(MODIFIED_SUBMIT_FILE);

# Rename The Original Submit File

rename($File_Name . ".modified", $File_Name);

}

#-------------------------------#
# Modify UMUI Submit Clear File #
#-------------------------------#

sub Modify_UMUI_Submit_Clear_File()

{

# Modifies The UMUI Submit Clear File

my $Loop = -1;
my $Permutation_Index = -1;
my $Line = "";
my $Permutation_Name = "";
my $Extended_UM_Job_Name = "";
my $File_Name = "";
my @Lines;

# Store The Passed Parameter

$Permutation_Name = $_[0];
$Permutation_Index = $_[1];

# Store The File Name

$File_Name = $Permutation_Name . "/umuisubmit_clr";

# Check If The File Exists

if (-e $File_Name)

   {

# Open The Submit File

   open(SUBMIT_FILE, "<" . $File_Name);

# Read In All Of The Lines

   @Lines = <SUBMIT_FILE>;

# Close The Submit File

   close(SUBMIT_FILE);

# Create The Modified Submit File

   open(MODIFIED_SUBMIT_FILE, ">" . $File_Name . ".modified");

# Initialise The Line Counter

   $Loop = 0;

# Extract The First Line

   $Line = $Lines[$Loop];

# Search For "PBS -o" On Hector
# Search For "#@ output" On Monsoon

   if ($Running_On_MONSOON == 0)

      {

      while ($Line !~ /PBS -o/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

   else

      {

      while ($Line !~ /@ output/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

# The Line Now Contains "PBS -o" Or "@ output"

# Remove The End Of Line Symbol

   chomp($Line);

# Append The Permutation Index

   $Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Modfy Standard Error If You Are On Monsoon

   if ($Running_On_MONSOON == 1)

      {

# Copy The Line For Standard Error As Well

      $Line =~ s/output/error/;

# Write Out The Modified Line

      print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

      $Loop++;

      }

# Extract The Next Line

   $Line = $Lines[$Loop];

# Search For "export UM_RDATADIR"

   while ($Line !~ /export UM_RDATADIR/)

      {

# Write The Line Out Unchaged

      print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

      $Loop++;

# Extract The Next Line

      $Line = $Lines[$Loop];

      }

# The Line Now Contains "export UM_RDATADIR"

# Perform The Substitution

   $Line =~ s/RUNID/{RUNID}/;

# Remove The End Of Line Symbol

   chomp($Line);

# Append The Permutation Index

   $Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

# Search For "qsub" On Hector
# Search for "llsubmit" On monsoon

   if ($Running_On_MONSOON == 0)

      {

      while ($Line !~ /qsub/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

   else

      {

      while ($Line !~ /llsubmit/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

# The Line Now Contains "qsub" Or "llsubmit"

# Add The Permutation Index To The UM Job Name

   $Extended_UM_Job_Name = $UM_Job_Name . ".x" . $Permutation_Index;

# Perform The Substitution

   $Line =~ s/$UM_Job_Name/$Extended_UM_Job_Name/g;

# Write Out The Modified Line

   print(MODIFIED_SUBMIT_FILE $Line);

# Close The Modified Submit File

   close(MODIFIED_SUBMIT_FILE);

# Rename The Original Submit File

   rename($File_Name . ".modified", $File_Name);

   }

}

#---------------------------------#
# Modify UMUI Submit Compile File #
#---------------------------------#

sub Modify_UMUI_Submit_Compile_File()

{

# Modifies The UMUI Submit Compile File

my $Loop = -1;
my $Permutation_Index = -1;
my $Line = "";
my $Permutation_Name = "";
my $Extended_UM_Job_Name = "";
my $File_Name = "";
my @Lines;

# Store The Passed Parameter

$Permutation_Name = $_[0];
$Permutation_Index = $_[1];

# Store The File Name

$File_Name = $Permutation_Name . "/umuisubmit_compile";

# Check If The File Exists

if (-e $File_Name)

   {

# Open The Submit File

   open(SUBMIT_FILE, "<" . $File_Name);

# Read In All Of The Lines

   @Lines = <SUBMIT_FILE>;

# Close The Submit File

   close(SUBMIT_FILE);

# Create The Modified Submit File

   open(MODIFIED_SUBMIT_FILE, ">" . $File_Name . ".modified");

# Initialise The Line Counter

   $Loop = 0;

# Extract The First Line

   $Line = $Lines[$Loop];

# Search For "PBS -o" On Hector
# Search For "#@ output" On Monsoon

   if ($Running_On_MONSOON == 0)

      {

      while ($Line !~ /PBS -o/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

   else

      {

      while ($Line !~ /@ output/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

# The Line Now Cotains "PBS -o" Or "@ output"

# Remove The End Of Line Symbol

   chomp($Line);

# Append The Permutation Index

   $Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Modfy Standard Error If You Are On Monsoon

   if ($Running_On_MONSOON == 1)

      {

# Copy The Line For Standard Error As Well

      $Line =~ s/output/error/;

# Write Out The Modified Line

      print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

      $Loop++;

      }

#MY THE BLOCK BELOW HAS BEEN DISABLED FOLLOWING EMAIL FROM STEVEN PICKERING 06/08/2014 WHICH SAYS;
#MY   It sounds like the variable UM_RDATADIR doesn't exist in the basis file.
#MY   If it doesn't, you will need to remove those few lines, it is basically
#MY   searching for those lines but can't find them.
#MY   But be careful, that variable may be important to other versions of the UM.
#MY
#MY# Extract The Next Line
#MY
#MY   $Line = $Lines[$Loop];
#MY
#MY# Search For "export UM_RDATADIR" 
#MY
#MY   while ($Line !~ /export UM_RDATADIR/)
#MY
#MY      {
#MY
#MY# Write The Line Out Unchaged
#MY
#MY      print(MODIFIED_SUBMIT_FILE $Line);
#MY
#MY# Increment The Number Of Lines
#MY
#MY      $Loop++;
#MY
#MY# Extract The Next Line
#MY
#MY      $Line = $Lines[$Loop];
#MY
#MY      }

# The Line Now Contains "export UM_RDATADIR"

# Perform The Substitution

   $Line =~ s/RUNID/{RUNID}/;

# Remove The End Of Line Symbol

   chomp($Line);

# Append The Permutation Index

   $Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Process The Remaining Lines

   while ($Loop < @Lines)

      {

# Write Out The Remaining Line

      print(MODIFIED_SUBMIT_FILE $Lines[$Loop]);

# Increment The Number Of Lines

      $Loop++;

      }

# Close The Modified Submit File

   close(MODIFIED_SUBMIT_FILE);

# Rename The Original Submit File

   rename($File_Name . ".modified", $File_Name);

   }

}

#-----------------------------#
# Modify UMUI Submit Run File #
#-----------------------------#

sub Modify_UMUI_Submit_Run_File()

{

# Modifies The UMUI Submit Run File

my $Loop = -1;
my $Permutation_Index = -1;
my $Line = "";
my $Permutation_Name = "";
my $Extended_UM_Job_Name = "";
my $File_Name = "";
my @Lines;

# Store The Passed Parameter

$Permutation_Name = $_[0];
$Permutation_Index = $_[1];

# Store The File Name

$File_Name = $Permutation_Name . "/umuisubmit_run";

# Open The Submit File

open(SUBMIT_FILE, "<" . $File_Name);

# Read In All Of The Lines

@Lines = <SUBMIT_FILE>;

# Close The Submit File

close(SUBMIT_FILE);

# Create The Modified Submit File

open(MODIFIED_SUBMIT_FILE, ">" . $File_Name . ".modified");

# Initialise The Line Counter

$Loop = 0;

# Extract The First Line

$Line = $Lines[$Loop];

# Search For "PBS -o" On Hector
# Search For "#@ output" On Monsoon

   if ($Running_On_MONSOON == 0)

      {

      while ($Line !~ /PBS -o/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

   else

      {

      while ($Line !~ /@ output/)

         {

# Write The Line Out Unchaged

         print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

         $Loop++;

# Extract The Next Line

         $Line = $Lines[$Loop];

         }

      }

# The Line Now Cotains "PBS -o" Or "@ output"

# Remove The End Of Line Symbol

chomp($Line);

# Append The Permutation Index

$Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

$Loop++;

# Modfy Standard Error If You Are On Monsoon

   if ($Running_On_MONSOON == 1)

      {

# Copy The Line For Standard Error As Well

      $Line =~ s/output/error/;

# Write Out The Modified Line

      print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

      $Loop++;

      }

# Extract The Next Line

$Line = $Lines[$Loop];

# Search For "export UMRUN_OUTPUT"

while ($Line !~ /export UMRUN_OUTPUT/)

   {

# Write The Line Out Unchaged

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

   }

# The Line Now Contains "export UMRUN_OUTPUT"

# Remove The End Of Line Symbol

chomp($Line);

# Append The Permutation Index

$Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

$Loop++;

# Extract The Next Line

$Line = $Lines[$Loop];

# Search For "export JOBDIR"

while ($Line !~ /export JOBDIR/)

   {

# Write The Line Out Unchaged

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

   }

# The Line Now Contains "export JOBDIR"

# Remove The End Of Line Symbol

chomp($Line);

# Append The Permutation Index

$Line .= (".x" . $Permutation_Index . "\n");

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

$Loop++;

# Extract The Next Line

$Line = $Lines[$Loop];

# Search For "UM_DATAW="

while ($Line !~ /UM_DATAW=/)

   {

# Write The Line Out Unchaged

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

   }

# The Line Now Contains "UM_DATAW="

# Calculate The Extended Job Name

$Extended_UM_Job_Name = "{RUNID}.x" . $Permutation_Index;

# Perform The Substitution

$Line =~ s/RUNID/$Extended_UM_Job_Name/;

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

$Loop++;

# Extract The Next Line

$Line = $Lines[$Loop];

# Search For "UM_DATAM="

while ($Line !~ /UM_DATAM=/)

   {

# Write The Line Out Unchaged

   print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

   $Loop++;

# Extract The Next Line

   $Line = $Lines[$Loop];

   }

# The Line Now Contains "UM_DATAM="

# Calculate The Extended Job Name

$Extended_UM_Job_Name = "{RUNID}.x" . $Permutation_Index;

# Perform The Substitution

$Line =~ s/RUNID/$Extended_UM_Job_Name/;

# Write Out The Modified Line

print(MODIFIED_SUBMIT_FILE $Line);

# Increment The Number Of Lines

$Loop++;

# Process The Remaining Lines

while ($Loop < @Lines)

   {

# Write Out The Remaining Line

   print(MODIFIED_SUBMIT_FILE $Lines[$Loop]);

# Increment The Number Of Lines

   $Loop++;

   }

# Close The Modified Submit File

close(MODIFIED_SUBMIT_FILE);

# Rename The Original Submit File

rename($File_Name . ".modified", $File_Name);

}

#----------------------#
# Read Parameters List #
#----------------------#

sub Read_Parameters_List()

{

# Reads The Parameters List

my $Loop = -1;
my $Line = "";
my @Lines;

# Open The File

open(PARAM_FILE, "<Parameters.list");

# Read In All Of The Lines

@Lines = <PARAM_FILE>;

# Close The File

close(PARAM_FILE);

# Process The Lines Of Text

for ($Loop = 0; $Loop < @Lines; $Loop++)

   {

# Extract A Line

   $Line = $Lines[$Loop];
   chomp($Line);

# Split The Line Into Tokens

   @Permutation_Tokens = split(/\|/, $Line);

# Check That There Are The Correct Number Of Tokens Present
# There Should Be "Number_Of_Variables" Tokens Per Line

   if (@Permutation_Tokens != $Number_Of_Variables)

      {

# Display The Message

      print("Error In Perturbation " . ($Loop + 1) . ": " . $Line . "\n");
      print("Incorrect Number Of Tokens Found, There Should Be " . $Number_Of_Variables . " But " . @Permutation_Tokens . " Were Found!!!\n");
      print("\n");

# Exit The Script

      exit(1);

      }

# Create The Perturbation Directory Name

   $Perturbation_Directory_Name = $UM_Directory_Name . ".x" . ($Loop + 1);

# Display The Name Of The Perturbation Directory

   print("Creating Pertubation Directory " . $Perturbation_Directory_Name . "\n");

# Create The Perturbation Directory

   mkdir($Perturbation_Directory_Name);

# Copy The Files From The Original UM Job Directory

   &Copy_Directory($UM_Directory_Name, $Perturbation_Directory_Name);

# Modify The All Control File

   &Modify_All_Control_File($Perturbation_Directory_Name);

# Modify The Atompshere Control File

   &Modify_Atmosphere_Control_File($Perturbation_Directory_Name);

## LRE 23rdSep2015. Added this subroutine above.
# Modify The Shared Control File

   &Modify_Shared_Control_File($Perturbation_Directory_Name);

# Modify The General Control File

   &Modify_General_Control_File($Perturbation_Directory_Name);

# Modify The UMUI Submit Clear File

   &Modify_UMUI_Submit_Clear_File($Perturbation_Directory_Name, $Number_Of_Perturbations + 1);

# Modify The UMUI Submit Compile File

   &Modify_UMUI_Submit_Compile_File($Perturbation_Directory_Name, $Number_Of_Perturbations + 1);

# Modify The UMUI Submit Run File

   &Modify_UMUI_Submit_Run_File($Perturbation_Directory_Name, $Number_Of_Perturbations + 1);

# Create The Perturbation Work Directory Name

   $Perturbation_Work_Directory_Name = $UM_Work_Directory_Name . ".x" . ($Loop + 1);

# Display The Name Of The Perturbation Directory

   print("Creating Pertubation Work Directory " . $Perturbation_Work_Directory_Name . "\n");

# Create The Perturbation Work Directory

   mkdir($Perturbation_Work_Directory_Name);
   mkdir($Perturbation_Work_Directory_Name . "/bin");

# Copy The Files From The UM Work Directory

#MY   print("Calling Copy_Directory 2nd time!\n");
#MY   print($UM_Work_Directory_Name . "/bin\n");           #/home/myosh/data_output/tdsqa/bin
#MY   print($Perturbation_Work_Directory_Name . "/bin\n"); #/home/myosh/data_output/tdsqa.x1/bin

   &Copy_Directory($UM_Work_Directory_Name . "/bin", $Perturbation_Work_Directory_Name . "/bin");

#MY   print(" 2nd call to Copy_Directory finished\n");

# Modify The Submit Check File

   &Modify_Submit_Check_File($Perturbation_Work_Directory_Name . "/bin");

# Increment The Number Of Perturbations

   $Number_Of_Perturbations++;

   }

# Display The Number Of Perturbations

print("\n");
print($Number_Of_Perturbations . " Perturbations Found\n");
print("\n");

}

#--------------------#
# Read Variable List #
#--------------------#

sub Read_Variable_List()

{

# Reads The Variable List

my $Loop = -1;
my $Line = "";
my $Key = "";
my $Value = "";
my $Control_File = "";
my $Name_List = "";
my $Variable_Name = "";
my @Lines;
my @Tokens;

# Open The File

open(VAR_FILE, "<Variables.list");

# Read In All Of The Lines

@Lines = <VAR_FILE>;

# Close The File

close(VAR_FILE);

# Process The Lines Of Text

for ($Loop = 0; $Loop < @Lines; $Loop++)

   {

# Extract A Line

   $Line = $Lines[$Loop];
   chomp($Line);

# Split The Line Into Three Tokens

   @Tokens = split(/ /, $Line);

# Store The Tokens

   $Control_File = $Tokens[0];
   chomp($Control_File);
   $Name_List = $Tokens[1];
   chomp($Name_List);
   $Variable_Name = $Tokens[2];
   chomp($Variable_Name);

# Check That The First Token Is A Valid Control File

   if (&Is_Valid_Control_File_Name($Control_File) == 0)

      {

# Display The Message

      print($Control_File . " Isn't A Valid Control File Name!!!\n");

# Exit The Script

      exit(1);

      }

# At This Point, We Know We Have A Valid Control File Name

# Check That The Second Token Is A Valid Namelist Name For That Control File

   if (((&Is_Valid_All_Control_Namelist_Name($Name_List) == 0) && ($Control_File eq "CNTLALL")) ||
       ((&Is_Valid_Atmosphere_Control_Namelist_Name($Name_List) == 0) && ($Control_File eq "CNTLATM")) ||
       ((&Is_Valid_General_Control_Namelist_Name($Name_List) == 0) && ($Control_File eq "CNTLGEN")))

      {

# Display The Message

      print($Name_List . "Isn't A Valid Namelist Name For The Control File " . $Control_File . "!!!\n");

# Exit The Script

      exit(1);

      }

# At This Point, We Know We Have A Valid Namelist Name 

# Increment The Number Of Variables

   $Number_Of_Variables++;

# Create The Key

   $Key = $Control_File . "_" . $Name_List . "_" . $Variable_Name;

# Store The Variable's Offset

   $Variables_Offset{$Key} = ($Loop + 1);

   }

# Display The Number Of Variables Found

print($Number_Of_Variables . " Variables Found\n");
print("\n");

}

#--------------#
# Main Program #
#--------------#

# Ask For The Job Name

&Ask_For_UM_Job_Name;

# Read In The Variables List

&Read_Variable_List;

# Read In The Parameters List

&Read_Parameters_List;

# Exit The Script

exit(0);

