#!/usr/bin/env python2.7
'''

Name: perturn_ini.py

Description:
    Perturb a certain parameter in the namelist by creating a new rose-app.conf
    file in the app/um/opt folder. A directive in the suite.rc file points the
    atmos_main task to the relevant altered rose-app.conf file.

Required Arguments:
    path_in: location of configuration file, referred to as CONF_LOC by cylc
    ens: The ensemble number, referred to as LABEL or as CYLC_TASK_PARAM_x by
         cylc

'''

import sys
import os
import re
import shutil
from tempfile import mkstemp

class ArgumentsError(Exception):
    '''
    Exception raised when there is an error detected in the argument list.
    '''
    def __init__(self, msg):
        print ('[FATAL ERROR] : %s' % msg )
        raise SystemExit

def sed(pattern, replace, source, dest=None, count=0):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.
    """

    fin = open(source, 'r')
    num_replaced = count

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')

    for line in fin:
        out = re.sub(pattern, replace, line)
        fout.write(out)

        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()

    if not dest:
        shutil.move(name, source)


def main():
    ''' Main function '''

    conf_in = sys.argv[1]
    if not os.path.exists(conf_in):
        raise ArgumentsError('Template namelist file does not exist: ' + conf_in)

    path_out = sys.argv[2]
    ens = sys.argv[3]
    conf_out = os.path.join(path_out, 'opt/rose-app-ens' + ens +'.conf' )

    if os.path.exists(conf_out):
        raise ArgumentsError('Output conf file already exists, will not overwrite: '
                             + conf_out)
    else:
        if path_out and not os.path.isdir(path_out):
            raise ArgumentsError('Directory to write ' + conf_out
                                 + ' does not exist')

    newval=str(float(ens) + 25.3)

    sed("^ukca_aeros_volc_so2=[^\n]*", "ukca_aeros_volc_so2="+newval, conf_in, conf_out)

    print ("Changed value of ukca_aeros_volc_so2 to " + newval)


if __name__ == '__main__':
    main()
