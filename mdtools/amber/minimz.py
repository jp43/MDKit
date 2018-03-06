import os
import sys
import stat
import shutil
import subprocess
import argparse

import ambertools
from mdtools.utility import utils
from mdtools.utility import mol2
from mdtools.utility import reader

leap_default_settings = {'solvate': False, 'PBRadii': None, 'forcefield': 'leaprc.ff14SB'}

def do_minimization(file_r, files_l=None, keep_hydrogens=False):
    """ do_minimization(file_r, files_l=None, keep_hydrogens=False)

    Performs Amber minimization

    Parameters
    ----------
    file_r: filename for receptor (.pdb)
    files_l: list of filenames (.mol2) for ligand, when ligand-protein complex

    Steps
    -----
    antechamber, parmchk (ligand-protein complex)"""

    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = 'minimz'
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of ligand and receptor files
    file_r = os.path.abspath(file_r)
    if files_l:
        if isinstance(files_l, str):
            files_l = [files_l] # make it a list
        l_tmp = []
        for file_l in files_l:
            l_tmp.append(os.path.abspath(file_l))
        files_l = l_tmp

    # change working directory
    os.chdir(workdir)

    # only keep atom lines
    with open(file_r, 'r') as tmpf:
        with open('protein.pdb', 'w') as recf:
            for line in tmpf:
                # check if non-hydrogen atom line
                if line.startswith(('ATOM','TER')):
                    recf.write(line)
            # if last line not TER, write it
            if not line.startswith('TER'):
                recf.write('TER\n')

    # prepare receptor
    ambertools.prepare_receptor('protein.pdb', file_r, keep_hydrogens=keep_hydrogens)

    # amber minimization
    do_amber_minimization('protein.pdb', files_l)
    os.chdir(curdir)

def prepare_minimization_config_file(script_name, ligname):

    with open(script_name, 'w') as minf:
        script ="""In-Vacuo minimization with restraints
&cntrl
 imin=1,
 maxcyc=5000,
 ntb=0,
 ncyc=1000,
 ntmin=1,
 ntpr=5,
 cut=10.0,
 ibelly=1,
 bellymask=':%(ligname)s'
&end\n"""%locals()
        minf.write(script)

def prepare_and_minimize(output_file, ligname):

    # run tleap
    subprocess.check_output('tleap -f leap.in > /dev/null', shell=True, executable='/bin/bash')
    shutil.copyfile('start.inpcrd', 'start.rst')

    prepare_minimization_config_file('min.in', ligname)
    try:
        # run minimization
        utils.run_shell_command('sander -O -i min.in -o min.out -c start.inpcrd -p start.prmtop -ref start.rst -r end.inpcrd')
        shutil.move('end.inpcrd', 'start.inpcrd')
        status = 0 # minimization finished normaly
    except subprocess.CalledProcessError as e:
        # the minimization failed
        return e.returncode

    # get output configuration
    utils.run_shell_command('cpptraj -p start.prmtop -y start.inpcrd -x %s'%output_file)
    return status

def do_amber_minimization(file_r, files_l):

    tmpfile_l = 'lig.mol2'
    for idx, file_l in enumerate(files_l):
        shutil.copyfile(file_l, tmpfile_l)

        # do not regenerate the charges using antechamber when preparing the ligand
        ambertools.prepare_ligand(file_r, tmpfile_l, 'complex.pdb', charge_method=None)

        # prepare tleap config file
        ambertools.prepare_leap_config_file('leap.in', file_r, tmpfile_l, 'complex.pdb')

        ligname = ambertools.get_ligand_name('lig.mol2')
        status = prepare_and_minimize('complex-out.pdb', ligname)

        if status == 0:
            is_ligand = False
            # get ligand atom positions from complex file
            with open('complex-out.pdb', 'r') as cf:
                with open('protein-%s-out.pdb'%(idx+1), 'w') as recf:
                    with open('lig-out.pdb', 'w') as ligf:
                        for line in cf:
                            if line.startswith(('ATOM', 'HETATM')):
                                if line[17:20] == ligname:
                                    is_ligand = True
                            if is_ligand:
                                ligf.write(line)
                            else:
                                recf.write(line)
                            is_ligand = False

            mol2file = 'lig-%s-out.mol2'%(idx+1)
            mol2.pdb2mol2('lig-out.pdb', mol2file, file_l)

def create_arg_parser():

    parser = argparse.ArgumentParser(description="Run Amber Minimization")

    parser.add_argument('-l',
        type=str,
        dest='input_file_l',
        nargs='+',
        default=None,
        help = 'Ligand coordinate file(s): .mol2')

    parser.add_argument('-r',
        type=str,
        dest='input_file_r',
        required=True,
        help = 'Receptor coordinate file(s): .pdb')

    parser.add_argument('-removeh',
        dest='removeh',
        action='store_true',
        help = 'Remove Hydrogens')

    return parser

def run():

    parser = create_arg_parser()
    args = parser.parse_args()

    do_minimization(args.input_file_r, files_l=args.input_file_l, keep_hydrogens=not args.removeh)
