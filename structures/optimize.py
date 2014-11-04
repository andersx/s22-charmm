import sys
import os

import multiprocessing as mp

from numpy.random import normal, choice, seed

import sys
import os
from ref import energies, names

import numpy as np

CPE_FILE = "/home/andersx/dev/charmm_cpe/source/sccdftbint/sccdftbsrc/cpe.src"

case_line = 274
line_id = dict()
line_id[0] = case_line + 1
line_id[1] = line_id[0] + 9
line_id[2] = line_id[1] + 9
line_id[3] = line_id[2] + 9

atom_type_names = dict()
atom_type_names[0] = "H"
atom_type_names[1] = "C"
atom_type_names[2] = "O"
atom_type_names[3] = "N"

parameter_names = dict()
parameter_names[0] = "Z"
parameter_names[1] = "B"
parameter_names[2] = "Ru"
parameter_names[3] = "Rl"

footer = "foot_test-cpe.inp"

def get_energy(s22_id, footer):


    f = os.popen("./run %s %i" % (footer, s22_id))

    lines = f.readlines()

    e_monomer1 = float(lines[0].split()[2])
    e_monomer2 = float(lines[1].split()[2])
    e_complex  = float(lines[2].split()[2])

    e_interaction  = e_complex - e_monomer1 - e_monomer2

    return e_interaction

def replace_parameter(atom_type, parameter, multiplier revert=False):

    f = open(CPE_FILE, "r")
    lines = f.readlines()
    f.close()

    i = line_id[atom_type] + parameter

    tokens = lines[i].split()

    par = float(tokens[2]) * multiplier

    # print lines[i]
    # print lines[i][:24] + str(par) + "\n"

    lines[i] = lines[i][:24] + str(par) + "\n"

    if not revert:
        print "Changing %1s: %2s %10.7f -> %10.7f" % \
            (atom_type_names[atom_type], parameter_names[parameter], float(tokens[2]), par)
    else:
        print "Discarding %1s: %2s %10.7f -> %10.7f" % \
            (atom_type_names[atom_type], parameter_names[parameter], float(tokens[2]), par)


    f = open(CPE_FILE, "w")
    for line in lines:
        f.write(line)
    f.close()


output = mp.Queue()

def calc_energy(tid, output):

    irange = []

    if tid == 0: irange = range(81,87)
    if tid == 1: irange = range(87,93)
    if tid == 2: irange = range(93,98)
    if tid == 3: irange = range(98,103)

    de_list = []

    # print tid, irange, len(irange)
    for s22_id in irange:
        e_interaction = get_energy(s22_id, footer)
        de_interaction = e_interaction - energies[s22_id]

        print "%4i %-32s  %6.2f  %6.2f  %6.2f" % \
            (s22_id, names[s22_id][:-4], 
            energies[s22_id], e_interaction, de_interaction)

        de_list.append(de_interaction)

    output.put(de_list)


def compile_charmm():
    sys.stdout.write ('Compiling changes ... ')
    sys.stdout.flush()
    os.system("/home/andersx/projects/s22-charmm/structures/compile_charmm.sh")
    sys.stdout.write ('Done!\n')

if __name__ == "__main__":

    atom_types = [0, 1, 2, 3]
    parameters = [0, 1, 2, 3]

    iterations = 10

    # initial RMSD
    rmsd = 1.07883923066
    # rmsd = 0.0001

    pool = mp.Pool(processes=4)
    
    accept = 0

    for i in xrange(iterations):

        p = choice(parameters)
        a = choice(atom_types)

        m = normal(loc=1.0, scale=0.2)

        replace_parameter(a, p, m)
        compile_charmm()

        # de_list = [pool.map_async(calc_energy, args=(x,)) for x in  range(81, 102)]

        processes = [mp.Process(target=calc_energy, args=(x, output)) for x in range(4)]

        # Run processes
        for pr in processes:
            pr.start()
        
        # Exit the completed processes
        for pr in processes:
            pr.join()

        de_list = []

        for x in [output.get() for pr in processes]:
            de_list += x

        de_list = np.array(de_list)

        rmsd_new = np.sqrt(np.sum(de_list**2.0)/len(de_list))
        print "Max AE  [kcal/mol]  %6.2f" % np.amax(abs(de_list))
        print "RMSD    [kcal/mol]  %6.2f" % np.sqrt(np.sum(de_list**2.0)/len(de_list))
        print "Mean     [kcal/mol]  %6.2f" % np.average(de_list)
        print "MAD     [kcal/mol]  %6.2f" % np.mean(abs(de_list))

        print "Old RMSD: ", rmsd
        print "New RMSD: ", rmsd_new 
        print "At step number", i+1, "out of", iterations
        if rmsd_new < rmsd:
            print "ACCEPT"
            rmsd = rmsd_new
            accept += 1
        else:
            print "REJECT"
            replace_parameter(a, p, 1.0/m, revert=True)
        print "Acceptance ratio: %4.1f" % (accept/(i + 1.0)* 100.0)

compile_charmm() 
