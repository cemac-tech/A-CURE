#!/usr/bin/env python2.7
'''

Name: perturb_ppe.py

Description:
    Perturb a set parameter in the namelist by creating a new rose-app.conf
    file in the app/um/opt folder. A directive in the suite.rc file points the
    atmos_main task to the relevant altered rose-app.conf file. Values for the 
    parameter are determined by reading specific rows in a dataframe csv file

Required Arguments:
    input  : location of template configuration file, referred to as CONF_LOC by cylc.
             Template file is to contain every parameter to be perturbed as a 
             logical set to false, as each parameter has an associated logical in the 
             um namelists.
    ens    : The ensemble number list, referred to as ens by cylc. The elements of the
             list are provided by cylc as a space separated list of parameters
    output : location of the optional configuration files folder
    data   : location of the PPE dataframe, arranged with variables as column headers
             and one row per ensemble member, row indices used as the ensemble identifier
             and the character "x" in cells where the variable is not used for that ensemble

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

def file_val(fname,head,path_data):
    '''
    Returns number of lines in file which are not namelist block header lines
    Also verifies that the template.conf file has the correct format and has all
    correct variables present
    '''
    h,j,k,m=0,0,0,0
    headlog=[]
    for e in head:
        if not re.match("acure",e):
            headlog.append("l_acure_"+e+"=.false.")
        else:
            headlog.append("l_"+e+"=.false.")
    pattern=re.compile('|'.join(headlog))
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            if re.match("^\[namelist:([^\]]*?)\]",l): j+=1
            if re.match("^l_([^=]*?)=\.false\.$", l): h+=1
            if l=="\n": m+=1
            if pattern.search(l): k+=1
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
        if re.match("acure",name):
            rep["l_"+name+"=.false."]="l_"+name+"=.true."
        else:
            rep["l_acure_"+name+"=.false."]="l_acure_"+name+"=.true."

    pattern = re.compile("|".join(rep.keys()))

    with open (source, 'r') as fin:
        for line in fin:

            repl_flg=0

            out = pattern.sub(lambda match: rep[match.group(0)], line)
            if not "\n" in out:
                out = out + "\n"
                if ".false." in out:
                    num_replaced -=1

            outlst.append(out)
            if out != line:
                num_replaced += 1
                repl_flg += 1
                print (out)

    if num_replaced == len(ens_dict):
        for (name,val) in ens_dict.items():
            teststr="_"+name+"="
            for i, line in enumerate(outlst):
                if teststr in line:
                    outlst.insert(i+1,name+"="+val+"\n")
        with open(dest, 'w') as fout:
            for line in outlst:
                fout.write("%s" % line)
    else:
        for (name,val) in ens_dict.items():
            print("%s=%s\n" % (name,val))
        raise FileError('One or more of the expected parameters were not replaced.\n'
                        + 'Expected to replace '+str(len(ens_dict))+' but replaced '
                        + str(num_replaced)+" in " + dest + "\n")

def check_int(s):
    s = str(s)
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()  


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
                        default=[3, 63])

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

    """
    Check the input Arguments
    """

    path_in = args.input
    if not os.path.exists(path_in):
        raise ArgumentsError('Path to template namelist file does not exist: '
                             + path_in)

    if os.path.isdir(path_in):
        conf_in = os.path.join(path_in, 'template.conf')
    else:
        conf_in=path_in

    if not os.path.exists(conf_in):
        raise ArgumentsError('Template namelist file does not exist: ' + conf_in +"\n")
        
    path_out = args.output

    if args.output==".":
        path_out=os.path.join(path_in, 'opt')
    else:
        path_out=args.output

    if not os.path.exists(path_out):
        raise ArgumentsError('Directory to write opt conf file to'
                             + ' does not exist\n')
    if args.data==".":
        path_data=os.path.join(path_in, 'ppe_dataframe.csv')
    else:
        path_data=args.data

    if not os.path.exists(path_data):
        raise ArgumentsError('PPE dataframe csv file does not exist: ' + path_data +"\n")

    try:
        df = pd.read_csv(path_data)
    except:
        raise FileError("Cannot open the dataframe file "+path_data)

    head=list(df)
    chkset=set(head)
    if not len(chkset) == len(head):
        raise ArgumentsError("Duplicate values detected in the dataframe\n")

    if file_val(conf_in,head,path_data) != df.shape[1]:
        raise FileError("Mismatch between dataframe variable number in "+path_data+
                        " and template file "+conf_in +"\n")
           
    ens_lst = args.ens
    chkset=set(ens_lst)
    if not len(chkset) == len(ens_lst):
        raise ArgumentsError("Duplicate values detected in the ensemble list\n")

    if len(ens_lst) > df.shape[0]:
        raise ArgumentsError("More ensemble members requested than exist within dataframe\n")

    for ens in ens_lst:
        
        if not check_int(ens):
            raise ArgumentsError("Non-integer ensemble indices detected in ensemble member list")
        ens_int_val=int(ens)
        if ens_int_val < 0:
            raise ArgumentsError("Negative ensemble member index detected in ensemble member list")
        if ens_int_val > df.shape[0]:
            raise ArgumentsError("Ensemble index detected which is larger than maximum")

        conf_out = os.path.join(path_out, 'rose-app-ens_' + str(ens) +'.conf' )

        if os.path.exists(conf_out):
            print('Output conf file already exists, will not overwrite: '
                  + conf_out +"\n Skipping this ensemble member.\n")
            continue
        else:
            if path_out and not os.path.isdir(path_out):
                raise ArgumentsError('Directory to write ' + conf_out
                                     + ' does not exist')

        if not ens in df.index:
            raise ArgumentsError('Ensemble member '+str(ens)+' not found in dataframe\n')

        row=df.iloc[ens].to_dict()
        ens_dict={}
        for (name, newval) in row.items():
            if newval != "x":
                ens_dict[name]=newval


        subparam(ens_dict,conf_in,conf_out)



if __name__ == '__main__':
    main()
