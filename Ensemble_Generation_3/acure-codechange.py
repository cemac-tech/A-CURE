#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:26:04 2019

@author: chmcsy
"""


var_lst=[]

with open("./acure.txt") as f:
    for line in f:
        line=line.strip()
        var_lst.append(line)

j=850

out=open("./acure_meta.conf","w")

for var in var_lst:
    out.write("[namelist:run_ukca=l_%s]\n" % var)
    out.write ("compulsory=true\n")
    out.write ("description=Turn on the ACURE %s variable\n" % var)
    out.write ("help=This can be used in the absence of a PPE dataframe file\n")
    out.write ("    =to enable ACURE variables\n")
    out.write ("ns=namelist/UM Science Settings/section34/A-CURE\n")
    out.write ("sort-key=xx%d\n" % j)
    out.write ("trigger=namelist:run_ukca=%s: .true.;\n" % var)
    out.write ("type=logical\n\n")
    j+=1

for var in var_lst:
    out.write ("[namelist:run_ukca=%s]\n" % var)
    out.write ("compulsory=true\n")
    out.write ("description=Value for the ACURE %s variable\n" % var)
    out.write ("help=This can be used in the absence of a PPE dataframe file\n")
    out.write ("    =to set values for ACURE variables\n")
    out.write ("ns=namelist/UM Science Settings/section34/A-CURE/PPE Variables\n")
    out.write ("range=0.0:100.0\n")
    out.write ("sort-key=xx%d\n" % j)
    out.write ("type=real\n\n")
    j+=1

out.close()

out2=open("./acure_macro.py","w")

out2.write('class vn112_t4751(rose.upgrade.MacroUpgrade):\n')

out2.write('    """Upgrade macro for ticket #4751 by christophersymonds."""\n\n')

out2.write('    BEFORE_TAG = "vn11.2"\n')
out2.write('    AFTER_TAG = "vn11.2_t4751"\n\n')

out2.write('    def upgrade(self, config, meta_config=None):\n')
out2.write('        \"\"\"Upgrade a UM runtime app configuration with ACURE PPE variables.\"\"\"\n\n')

for var in var_lst:
    out2.write('        self.add_setting(config,[\"namelist:run_ukca\",\n')
    out2.write('                                 \"l_%s"],".false.")\n\n' % var)

    out2.write('        self.add_setting(config,[\"namelist:run_ukca\",\n')
    out2.write('                                 \"%s\"],\"0.0\")\n\n' % var)

out2.write('        return config, self.reports')

out2.close()

out3=open("./acure_optionsmod.F90","w")

out3.write('!******************ACURE LOGICALS*************************!\n\n')

for var in var_lst:
    out3.write('LOGICAL :: l_%s           = .FALSE.\n' % var)

out3.write('\n\n!*****************End ACURE LOGICALS**********************!\n\n')

out3.write('!******************ACURE REALS*************************!\n\n')

for var in var_lst:
    out3.write('REAL :: %s           = rmdi\n' % var)

out3.write('\n\n!*****************End ACURE REALS**********************!\n\n')

out3.write('Insert the following at the end of the namelist block starting \"NAMELIST/run_ukca/\"\n')

for var in var_lst:
    out3.write("         l_%s, %s,                  &\n" % (var,var))

out3.write('\n\nInsert the following before the line \"CALL umPrint(\'- - - - - - end of namelist - - - - - -\', &)\"\n\n')

for var in var_lst:
    out3.write('WRITE(lineBuffer,\'(A,L1)\')\'  l_%s = \', l_%s\n' % (var,var))
    out3.write('CALL umPrint(lineBuffer,src=\'ukca_option_mod\')\n')
    out3.write('WRITE(lineBuffer,\'(A,F16.4)\')\' %s = \',   %s\n' % (var,var))
    out3.write('CALL umPrint(lineBuffer,src=\'ukca_option_mod\'))\n')

out3.write('\n\nInsert the following in \"TYPE my_namelist"\n\n')

for var in var_lst:
    out3.write('  REAL :: %s\n' % var)

out3.write('\n\n')

for var in var_lst:
    out3.write('  LOGICAL :: l_%s\n' % var)

out3.write('\nAfter \"IF (mype == 0) THEN\", insert\n\n')

for var in var_lst:
    out3.write('  my_nml %% %s = %s\n' % (var,var))

out3.write('\n\n')

for var in var_lst:
    out3.write('  my_nml %% l_%s = l_%s\n' % (var,var))

out3.write('\nAfter \"IF (mype /= 0) THEN\", insert\n\n')

for var in var_lst:
    out3.write('  %s = my_nml %% %s\n' % (var,var))

out3.write('\n\n')

for var in var_lst:
    out3.write('  %s = my_nml %% l_%s\n' % (var,var))

out3.close()
