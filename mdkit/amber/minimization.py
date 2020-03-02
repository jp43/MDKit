import os
import sys
import stat
import shutil
import subprocess
import argparse

import ambertools
from mdkit.utility import utils
from mdkit.utility import mol2
from mdkit.utility import reader

leap_default_settings = {'solvate': False, 'PBRadii': None, 'forcefield': 'leaprc.ff14SB'}

def do_minimization_after_docking(file_r, files_l, keep_hydrogens=False, charge_method=None, ncyc=None, maxcyc=None, cut=None, amber_version='14'):
    """ do_minimization_after_docking(file_r, files_l, keep_hydrogens=False)

    Performs minimization with Amber after docking with MOBPredictor

    Parameters
    ----------
    file_r: file for receptor (.pdb)
    files_l: list of file (.mol2) containing ligand poses

    Steps
    -----
    antechamber, parmchk (ligand-protein complex)"""

    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = 'em'
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of ligand and receptor files
    file_r_abs = os.path.abspath(file_r)
    files_l_abs = []
    if isinstance(files_l, str):
        files_l_abs.append(os.path.abspath(files_l))
    else:
        for file_l in files_l:
            files_l_abs.append(os.path.abspath(file_l))

    # change working directory
    os.chdir(workdir)

    # prepare receptor
    ambertools.prepare_receptor('protein.pdb', file_r_abs, keep_hydrogens=keep_hydrogens)

    # amber minimization
    do_amber_minimization_after_docking('protein.pdb', files_l_abs, charge_method=charge_method, ncyc=ncyc, maxcyc=maxcyc, cut=cut, version=amber_version)
    os.chdir(curdir)

def prepare_minimization_config_file(script_name, ligname, ncyc=None, maxcyc=None, cut=None):
    """Prepare input file for constrained minimization"""
    with open(script_name, 'w') as minf:
        script ="""In-Vacuo minimization with constraints
&cntrl
imin=1,
ntmin=1,
maxcyc=%(maxcyc)s,
ncyc=%(ncyc)s,
ntb=0,
ntpr=100,
cut=%(cut)s,
ibelly=1,
bellymask=':%(ligname)s'
&end\n"""%locals()
        minf.write(script)

def prepare_and_minimize(output_file, ligname, charge_method=None, ncyc=None, maxcyc=None, cut=None):
    """Prepare and run in-vacuo minimization with constraints on all protein atoms"""

    # run tleap
    subprocess.check_output('tleap -f leap.in > /dev/null', shell=True, executable='/bin/bash')
    shutil.copyfile('start.inpcrd', 'start.rst')

    prepare_minimization_config_file('min.in', ligname, ncyc=ncyc, maxcyc=maxcyc, cut=cut)
    try:
        # run minimization
        utils.run_shell_command('sander -O -i min.in -o min.out -c start.inpcrd -p start.prmtop -ref start.rst -r end.inpcrd')
        shutil.move('end.inpcrd', 'start.inpcrd')
        status = 0 # minimization finished normally
    except subprocess.CalledProcessError as e:
        # the minimization failed
        return e.returncode

    # get output configuration
    utils.run_shell_command('cpptraj -p start.prmtop -y start.inpcrd -x %s'%output_file)
    return status

def do_amber_minimization_after_docking(file_r, files_l, charge_method=None, ncyc=None, maxcyc=None, cut=None, version='14'):

    for idx, file_l in enumerate(files_l):
        if idx == 0:
            shutil.copyfile(file_l, 'ligand.mol2')
            # re-generate charges
            ambertools.prepare_ligand(file_r, 'ligand.mol2', 'complex.pdb', charge_method=charge_method, version=version)
            shutil.copyfile('ligand.mol2', 'ligand-ref.mol2')
            ligname = ambertools.get_ligand_name('ligand-ref.mol2')
        else:
            # if not first one, do not regenerate the charges, copy charges generated the first time
            coords_l = mol2.get_coordinates(file_l)
            struct = mol2.Reader('ligand-ref.mol2').next()
            struct = mol2.replace_coordinates(struct, coords_l)
            mol2.Writer().write('ligand.mol2', struct)
            ambertools.prepare_ligand(file_r, 'ligand.mol2', 'complex.pdb', charge_method=None, version=version)

        # prepare tleap config file
        ambertools.prepare_leap_config_file('leap.in', file_r, 'ligand.mol2', 'complex.pdb', solvate=False, version=version)
        status = prepare_and_minimize('complex-out.pdb', ligname, ncyc=ncyc, maxcyc=maxcyc, cut=cut)

        if status == 0:
            is_ligand = False
            # get ligand atom positions from complex file
            with open('complex-out.pdb', 'r') as cf:
                with open('ligand-out.pdb', 'w') as ligf:
                    for line in cf:
                        if line.startswith(('ATOM', 'HETATM')):
                            if line[17:20] == ligname:
                                is_ligand = True
                        if is_ligand:
                            ligf.write(line)
                        is_ligand = False

            basename_l = os.path.basename(file_l)
            suffix, ext = os.path.splitext(basename_l)
            mol2file_out = suffix + '-out' + ext
            mol2.pdb2mol2('ligand-out.pdb', mol2file_out, file_l)
