#!/usr/bin/env python2.7
'''

Name: perturb_ppe.py

Description:
    Overwrite a set of parameters in the um app by creating a new rose-app.conf
    file in the app/um/opt folder, or custom alternative folder. A directive in
    the suite.rc file points the atmos_main task to the relevant altered
    rose-app.conf file. Values for the parameters are determined by reading
    specific rows in a dataframe csv file

Required Arguments:
    input  : location of template configuration file, referred to as TMPL_LOC by cylc.
             Template file is to contain every parameter to be perturbed as a
             logical set to false, as each parameter has an associated logical in the
             um namelists.
    ens    : The ensemble number list, referred to as ens by cylc. The elements of the
             list are provided by cylc as a space separated list of parameters
    output : location of the optional configuration files folder
    data   : location of the PPE dataframe, arranged with variables as column headers
             and one row per ensemble member, row indices used as the ensemble identifier
             and the character "x" in cells where the variable is not used for that ensemble member

'''
import sys
import os
import re
import argparse
import pandas as pd

class ArgumentsError(Exception):
    '''
    Exception raised when there is an error detected in the argument list.
    '''
    def __init__(self, msg):
        sys.stderr.write('[FATAL ERROR] : %s' % msg )
        sys.exit(9)

class FileError(Exception):
    '''
    Exception raised when contents of files are not as expected
    '''
    def __init__(self,msg):
        sys.stderr.write('[FILE ERROR] : %s' % msg )
        sys.exit(9)

def input_val(args):
    '''
    Returns the validated variables obtained through processing the input
    arguments.
    Args:
        args (list): List of arguments from arg parser
    Returns:
        conf_in (str)         : path of input file with filename
        path_out (str)        : path of output file, directory only
        df (pandas dataframe) : dataframe of ensemble member ppe values
        ens_lst (list)        : requested ensemble members
    '''

    # Validate the input path_in and generate conf_in

    path_in = args.input
    if not os.path.exists(path_in):
        raise ArgumentsError('Path to template namelist file does not exist: '
                             + path_in)

    if os.path.isdir(path_in):
        conf_in = os.path.join(path_in, 'template.conf')
    else:
        conf_in=path_in
        path_in=os.path.split(conf_in)[0]

    if not os.path.exists(conf_in):
        raise ArgumentsError('Template namelist file does not exist: ' + conf_in +"\n")

    # Validate the output path, or generate it from default value

    if args.output==".":
        path_out=os.path.join(path_in, 'opt')
    else:
        path_out=args.output

    if not os.path.exists(path_out):
        raise ArgumentsError('Directory to write opt conf file(s) to'
                             + ' does not exist\n')

    if path_out and not os.path.isdir(path_out):
        raise ArgumentsError('Directory to write opt conf file(s) to'
                             + ' does not exist\n')

    # Validate location of PPE dataframe

    if args.data==".":
        path_data=os.path.join(path_in, 'ppe_dataframe.csv')
    else:
        path_data=args.data

    if not os.path.exists(path_data):
        raise ArgumentsError('PPE dataframe csv file does not exist: ' + path_data +"\n")

    # Load PPE csv file into pandas dataframe

    try:
        df = pd.read_csv(path_data)
    except:
        raise FileError("Cannot open the dataframe file "+path_data)

    head=list(df)
    chkset=set(head) # Ensure no duplicates in header of dataframe
    if not len(chkset) == len(head):
        raise ArgumentsError("Duplicate values detected in the dataframe\n")

    # Validate that the csv file and the template file are compatible
    # 22 columns are disregarded as these are regional values that are enabled or disabled in groups

    k=file_val(conf_in,head)
    if k != df.shape[1]-22:
        raise FileError("Number of non-regional variables in "+path_data+" file not as expected\n"+
                        str(k)+" non-regional variables found in both"+conf_in+" and csv file but "+
                        str(df.shape[1]-22)+" variables with 22 regional variables found\n")

    # Validate the list of requested ensemble members

    ens_lst = args.ens
    chkset=set(ens_lst)
    if not len(chkset) == len(ens_lst):
        raise ArgumentsError("Duplicate values detected in the ensemble list\n")

    if len(ens_lst) > df.shape[0]:
        raise ArgumentsError("More ensemble members requested than exist within dataframe\n")

    for ens in ens_lst:

        # Check each requested ensemble member exists, is an integer, and is within the bounds of the dataframe

        if not isinstance(ens, int):
            raise ArgumentsError("Non-integer ensemble indices detected in ensemble member list")
        ens_int_val=int(ens)
        if ens_int_val < 0:
            raise ArgumentsError("Negative ensemble member index detected in ensemble member list")
        if ens_int_val > df.shape[0]:
            raise ArgumentsError("Ensemble index detected which is larger than maximum")
        if not ens in df.index:
            raise ArgumentsError('Ensemble member '+str(ens)+' not found in dataframe\n')

    return path_out, conf_in, df, ens_lst

def file_val(fname,head):
    '''
    Returns number of lines in file template which are not namelist block header
    lines or blank lines. Also verifies that the template.conf file has the
    correct format and has all correct variables present, not counting regional
    perturbation variables.
    Args:
        fname (str): name and path of the template file
        head (list): header line of the ppe dataframe as read from csv file
    '''
    h,j,k,m=0,0,0,0
    headlog=[]
    for e in head:
        if re.search("_ems_",e) is not None:
            # Don't count the regional perturbations
            continue
        elif re.search("_so2_",e) is not None:
            # Don't count the regional perturbations
            continue
        elif re.search("acure",e) is not None:
            # for all acure only parameters
            headlog.append("l_"+e+"=.false.")
        else:
            # for the pre-existing parameters
            headlog.append("l_acure_"+e+"=.false.")

    pattern=re.compile('|'.join(headlog))
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            if re.match("\[namelist:([^\]]*?)\]",l): j+=1 #namelist headers
            if re.match("l_([^=]*?)=\.false\.$", l): h+=1 #logical placeholders
            if l=="\n": m+=1                              #blank lines
            if pattern.search(l): k+=1                    #lines matching csv header fields

    if not h == i-j-m:
        raise FileError(str(i-j-h-m)+" entries in the template file had wrong format.\n"
                        +"Entries should be assignment of false to logical switches"
                        +" for every PPE variable or namelist block headers.\n")

    if not h == k:
        raise FileError("Mismatch between dataframe variable values in "+path_data+
                        " and template file "+fname+" in "+str(h-k)+" entries.\n")
    return k

def subparam(ens_dict, source, dest):
    """Reads a source file and writes the destination file.
    Checks that required parameters are in source files
    if present, changes logical to true and writes entry for parameter

    Args:
        ens_dict (dict) : dict of variable names and associated new values
        source  (str): input filename
        dest (str): destination filename
    """

    outlst = []

    num_replaced = 0

    rep={}

    for name in ens_dict:
        if re.search("_ems_",name) is not None:
            #Skip the regional variables as no one-to-one associated logical
            continue
        elif re.search("_so2_",name) is not None:
            #Skip the regional variables as no one-to-one associated logical
            continue
        elif re.search("acure",name) is not None:
            # Standard acure only variables
            rep["l_"+name+"=.false."]="l_"+name+"=.true."
        else:
            # Pre-existing variables
            rep["l_acure_"+name+"=.false."]="l_acure_"+name+"=.true."

    pattern = re.compile("|".join(rep.keys()))

    with open (source, 'r') as fin:
        for line in fin:

            repl_flg=0
            if len(rep)==0:
                #To handle exception where ens_dict has no replacement values
                out=line
            else:
                # substitute logical values based on contents of ens_dict
                out = pattern.sub(lambda match: rep[match.group(0)], line)

            if not "\n" in out:
                # Prevent a line without a \n character being treated as a replacement
                out = out + "\n"
                if ".false." in out:
                    num_replaced -=1

            #Outlst used to build up the entries in the output conf file
            outlst.append(out)
            if out != line:
                num_replaced += 1
                repl_flg += 1
                print (out)

    if num_replaced == len(ens_dict):

        # Compatibility between template conf file and ensemble parameter set
        # No regional perturbations, as these have 5 or 6 parameters controlled
        # by a single logical switch.

        for (name,val) in ens_dict.items():
            teststr="_"+name+"="
            for i, line in enumerate(outlst):
                if teststr in line:
                    if (name=="acure_pcalc_index"):
                        # Special case where a pcalc file is selected by modifying the string for ukcaprec to include the pcalc_index (special case 2 above)
                        outlst.insert(i+1,"ukcaperc=/work/n02/n02/lre/pcalc_UKESM_11_1_default/RADAER_pcalc_"+val+".ukca\n")
                    else:
                        outlst.insert(i+1,name+"="+val+"\n")
                        #adds new value for parameter to output conf file
                        if (name=="m_ci" or name=="a_ent_1_rp"):
                            # Special cases where min and max values are set to the same as the parameter
                            outlst.insert(i+2,name+"_min="+val+"\n")
                            outlst.insert(i+3,name+"_max="+val+"\n")

        # Write out the output opt conf file
        with open(dest, 'w') as fout:
            for line in outlst:
                fout.write("%s" % line)

    else:

        # Regional perturbations exist, throwing off the compatibility between
        # template conf file and ensemble parameter set.

        if (any("_ems_" in key for key in ens_dict) or
           any("_so2_" in key for key in ens_dict)):
            noreg = [key for key,value in ens_dict.items() if ((not "_so2_" in key) and (not "_ems_" in key))]
            noreg_dict = {}
            for k in noreg: noreg_dict.update({k:ens_dict[k]})
            # Set a dict of parameters without the regional parameter values
            # Allowing compatibility check between template conf file and
            # the ensemble member parameter set
            if num_replaced == len(noreg_dict):
                for (name,val) in noreg_dict.items():
                    teststr="_"+name+"="
                    for i, line in enumerate(outlst):
                        if teststr in line:
                            # Proceed as before if not one of the regional perturbation
                            # special case variables
                            if (name!="acure_carb_ff_ems" and
                                name!="acure_carb_bb_ems" and
                                name!="acure_carb_res_ems" and
                                name!="acure_anth_so2"):
                                if (name=="acure_pcalc_index"):
                                    outlst.insert(i+1,"ukcaperc=/work/n02/n02/lre/pcalc_UKESM_11_1_default/RADAER_pcalc_"+val+".ukca\n")
                                else:
                                    outlst.insert(i+1,name+"="+val+"\n")
                                if (name=="m_ci" or name=="a_ent_1_rp"):
                                    outlst.insert(i+2,name+"_min="+val+"\n")
                                    outlst.insert(i+3,name+"_max="+val+"\n")

                            else:
                                # Add the regional perturbation values rather than the indicator
                                # variable value
                                temp_dict={}
                                j=1
                                regkeys=[key for key,value in ens_dict.items() if name+"_" in key]
                                for k in regkeys: temp_dict.update({k:ens_dict[k]})
                                for (name2,val2) in temp_dict.items():
                                    outlst.insert(i+j,name2+"="+val2+"\n")
                                    j +=1

                with open(dest, 'w') as fout:
                    for line in outlst:
                        fout.write("%s" % line)
            else:
                # Handle a mismatch between the template conf file and ensemble member
                # parameter list
                for (name,val) in ens_dict.items():
                    print("%s=%s\n" % (name,val))
                raise FileError('One or more of the expected parameters were not replaced.\n'
                                + 'Expected to replace '+str(len(ens_dict))+' but replaced '
                                + str(num_replaced)+" in " + dest + "\n")
        else:
            # Handle a mismatch between the template conf file and ensemble member
            # parameter list
            for (name,val) in ens_dict.items():
                print("%s=%s\n" % (name,val))
            raise FileError('One or more of the expected parameters were not replaced.\n'
                            + 'Expected to replace '+str(len(ens_dict))+' but replaced '
                            + str(num_replaced)+" in " + dest + "\n")

def main():
    ''' Main function '''

    parser = argparse.ArgumentParser(description=(
        'Allow a set of ensemble files to be created based on a  '
        'predefined dataframe'
        ))

    parser.add_argument('--input',
                        type=str,
                        help='Path to template configuration file',
                        default='.')

    parser.add_argument('--ens', '-e',
                        nargs="*",
                        type=int,
                        help='Space separated list of ensemble numbers',
                        default=[0, 7, 250])

    parser.add_argument('--output', '-o ',
                        type=str,
                        help='Path to rose-app optional conf files',
                        default='.')

    parser.add_argument('--data', '-d ',
                        type=str,
                        help='Path to PPE dataframe',
                        default=".")

    args = parser.parse_args()

    print("Input location     : %s" % args.input)
    print("Output location    : %s" % args.output)
    print("Dataframe Location : %s" % args.data)
    print("Ensemble Indices   : %s" % args.ens)

    #Check the input Arguments

    path_out, conf_in, df, ens_lst = input_val(args)

    for ens in ens_lst:

        conf_out = os.path.join(path_out, 'rose-app-ens_' + str(ens) +'.conf' )

        if os.path.exists(conf_out):
            print('Output conf file already exists, will not overwrite: '
                  + conf_out +"\n Skipping this ensemble member.\n")
            continue

        row=df.iloc[ens].to_dict()
        ens_dict={}
        for (name, newval) in row.items():
            if newval != "x":
                ens_dict[name]=newval

        subparam(ens_dict,conf_in,conf_out)

if __name__ == '__main__':
    main()
