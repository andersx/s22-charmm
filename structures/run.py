import sys
import os
from ref import energies, names

import numpy as np


def get_energy(s22_id, footer):


    f = os.popen("./run %s %i" % (footer, s22_id))

    lines = f.readlines()

    e_monomer1 = float(lines[0].split()[2])
    e_monomer2 = float(lines[1].split()[2])
    e_complex  = float(lines[2].split()[2])

    e_interaction  = e_complex - e_monomer1 - e_monomer2

    return e_interaction


# footer = "foot_test.inp"
# footer = "foot_dftb3.inp"
# footer = "foot_dftb3-d3.inp"
# footer = "foot_cpe-dftb3.inp"
# footer = "foot_cpe-dftb3-d3.inp"
footer = "foot_test-cpe.inp"

e_list = []
de_list = []

for s22_id in range(81, 102):

    e_interaction = get_energy(s22_id, footer)
    de_interaction = e_interaction - energies[s22_id]

    e_list.append(e_interaction)
    de_list.append(de_interaction)

    print "%4i %-32s  %6.2f  %6.2f  %6.2f" % \
        (s22_id, names[s22_id][:-4], 
        energies[s22_id], e_interaction, de_interaction)

e_list = np.array(e_list)
de_list = np.array(de_list)

print "Max AE  [kcal/mol]  %6.2f" % np.amax(abs(de_list))
print "RMSD    [kcal/mol]  %6.2f" % np.sqrt(np.sum(de_list**2.0)/len(de_list))
print "Mean     [kcal/mol]  %6.2f" % np.average(de_list)
print "MAD     [kcal/mol]  %6.2f" % np.average(abs(de_list))
print "MAD     [kcal/mol]  %6.2f" % np.mean(abs(de_list))



 








