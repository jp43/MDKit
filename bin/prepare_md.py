#!/usr/bin/env python
import os
import sys
import math
import shutil
import argparse
import subprocess
import numpy as np

from mdkit.utility import utils
from mdkit.namd import namdtools
from mdkit.amber import ambertools

# ----- SET CONSTANTS ------
# constants
kT = 0.5924849 # kT at 25C (298.15 K) in kcal/mol
na = 6.022140857e23 # Avogadro number

# units conversions
j2kcal = 2.3900573e-4
da2kg = 1.660539e-27
m2a = 1e10
ps2s = 1e-12

# steps included in MD preparation
known_steps = ['prep', 'min', 'nvt', 'npt', 'md']

# born radii for each implicit solvation scheme
known_pbradii = {1: 'mbondi', 2: 'mbondi2', 5: 'mbondi2', 7: 'bondi', 8: 'mbondi3'}

parser = argparse.ArgumentParser(description="Prepare MD simulations with Amber")

parser.add_argument('-r',
    dest='file_r',
    required=False,
    default=None,
    help="PDB input file for protein (1 structure)")

parser.add_argument('-l',
    dest='file_l',
    required=False,
    default=None,
    nargs='+',
    help="Input files for ligand (.pdb, .mol2)")

parser.add_argument('-addions',
    dest='addions',
    type=float,
    default=0.0,
    help="Add specific concentration of ions (Na+, Cl-) and neutralize the system")

parser.add_argument('-box',
    dest='box',
    type=str,
    default='parallelepiped',
    choices=['parallelepiped', 'octahedron'],
    help="Box type used for preparation")

parser.add_argument('-boxsize',
    dest='boxsize',
    type=float,
    default=10.0,
    help="Box size used for preparation")

parser.add_argument('-c',
    dest='charge_method',
    default='gas',
    help="Method used to assign partial charges when preparing the ligand (default: 'gas')")

parser.add_argument('-cut',
    dest='cut',
    default=None,
    help="Cutoff for non-bonded interactions in angstroms (default (explicit): 10.0, default (implicit, vacuo): 999.0)")

parser.add_argument('-hem',
    dest='hem',
    action='store_true',
    help="Skip generation of charges for unrecognized atoms with antechamber (similar to tutorial on HEME group: http://ambermd.org/tutorials/advanced/tutorial20/mcpbpy_heme.html)")

parser.add_argument('-keeph',
    dest='keeph',
    action='store_true',
    default=False,
    help="Keep hydrogens during amber preparation!")

parser.add_argument('-igb',
    dest='igb',
    type=int,
    default=5,
    help="GBSA index (used with implicit solvent, default: 5)")

parser.add_argument('-iwrap',
    dest='iwrap',
    type=int,
    default=0,
    help="Value of iwrap for production run (default: 1, i.e., the coordinates written to the restart and trajectory files are not wrapped into a primary box)")

parser.add_argument('-lb',
    dest='lower_bound',
    default=8.0,
    help="Lower bound for binding site")

parser.add_argument('-membrane',
    dest='membrane',
    action='store_true',
    help="Prepare simulation of membrane protein in lipids. Protein structure is assumed to be provided using CHARMM membrane builder (PDB output file from step 5).")

parser.add_argument('-namd',
    dest='namd',
    action='store_true',
    help="Prepare files for NAMD only.")

parser.add_argument('-np',
    dest='ncpus',
    type=int,
    default=None,
    help="Number of cpus used for the simulations (default (serial, largemem): 1, default (gpu): 16)")

parser.add_argument('-maxcyc_s',
    dest='maxcyc_s',
    default=2000,
    help="Minimization (solvent): maxcyc (default: 2000)")

parser.add_argument('-maxcyc',
    dest='maxcyc',
    default=5000,
    help="Minimization (full): maxcyc (default: 5000)")

parser.add_argument('-ncyc_s',
    dest='ncyc_s',
    default=500,
    help="Minimization (solvent): ncyc (default: 500)")

parser.add_argument('-ncyc',
    dest='ncyc',
    default=2000,
    help="Minimization (full): ncyc (default: 2000)")

parser.add_argument('-nstlim_h',
    dest='nstlim_h',
    default=10000,
    help="NVT: nstlim (default: 10000)")

parser.add_argument('-nstlim_e',
    dest='nstlim_e',
    default=20000,
    help="NPT: nstlim (default: 10000)")

parser.add_argument('-nstlim',
    dest='nstlim',
    default=100000,
    help="Production: nstlim (default: 100000)")

parser.add_argument('-ntpr',
    dest='ntpr',
    default=50,
    help="Production: ntpr and ntwr for restart files (default: 50)")

parser.add_argument('-ntwprt',
    dest='ntwprt',
    type=int,
    default=0,
    help="Atoms to include in trajectory file, 0: include solvent, 1: no solvent (default: 0)")

parser.add_argument('-ntwx',
    dest='ntwx',
    default=0,
    help="Production: ntwx (default: 0)")

parser.add_argument('-p',
    dest='partition',
    type=str,
    default="serial",
    help="Partition to run MD on")

parser.add_argument('-pbradii',
    dest='pbradii',
    action='store_true',
    default=False,
    help="Set PBRadii option consistent with igb value.")

parser.add_argument('-rst',
    dest='restraints',
    type=str,
    default=None,
    help="Restraints (None, backbone)")

parser.add_argument('-s',
    dest='script_name',
    required=False,
    default='run_md.sh',
    help="bash script filename")

parser.add_argument('-solvent',
    dest='solvent',
    type=str,
    default='explicit',
    choices=['explicit', 'implicit', 'vacuo'],
    help="solvent model (explicit, implicit, vacuo)")

parser.add_argument('-st',
    dest='step',
    nargs='+',
    default=known_steps,
    help="steps to be conducted (prep, min, nvt, npt, md)")

parser.add_argument('-temp',
    dest='temp',
    default=298.0,
    help="temperature (default: 298.0K)")

parser.add_argument('-w',
    dest='workdir',
    required=False,
    default=None,
    help="name of working directory")

parser.add_argument('-water',
    dest='water',
    type=str,
    default="TIP3P",
    help="water model")

args = parser.parse_args()

# set cutoff value if not set
if args.cut is None:
    if args.solvent == 'explicit':
        args.cut = 10.0
    else:
        args.cut = 999.0

if args.ncpus is None:
    if args.partition == 'gpu':
        args.ncpus = 16
    else:
        args.ncpus = 8

locals().update(args.__dict__)

# get amber version
amber_version = utils.check_amber_version()

# save current directory
pwd = os.getcwd()

if args.file_r:
    # get absolute paths of ligand and receptor files
    file_r_abspath = os.path.abspath(args.file_r)
else:
    raise ValueError('Receptor file should be provided!')

if args.file_l:
    files_l = [os.path.abspath(file_l) for file_l in args.file_l]

# ----- SET WORKING DIRECTORY ------
if args.workdir:
    workdir_r = os.path.abspath(args.workdir)
else:
    dir_r = os.path.dirname(file_r_abspath)
    file_r_prefix, ext = os.path.splitext(os.path.basename(file_r_abspath))

    workdir_r = dir_r + '/' + file_r_prefix

# ------- STEP 1: preparation -------
if 'prep' in args.step:
    shutil.rmtree(workdir_r, ignore_errors=True)
    os.mkdir(workdir_r)

    workdir = workdir_r + '/common'
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir)
    os.chdir(workdir)

    ambertools.prepare_receptor('protein.pdb', file_r_abspath, keep_hydrogens=args.keeph, membrane=args.membrane)

    # ligand preparation
    if args.file_l:
        files_l_new = []
        # preparing ligand(s)...
        for file_l in files_l:
            name = ambertools.get_ligand_name(file_l)

            prefix, ext = os.path.splitext(file_l)
            file_l_new = 'ligand_%s%s'%(name, ext)
            shutil.copyfile(file_l, file_l_new)
            files_l_new.append(file_l_new)

        files_l = files_l_new
        mol2files_l = ambertools.prepare_ligand('protein.pdb', files_l, 'complex.pdb', charge_method=args.charge_method, version=amber_version, skip_unrecognized_atoms=args.hem)
    else:
        mol2files_l = None

    if args.solvent == 'explicit' and not args.membrane:
        solvate = True
    else:
        solvate = False

    if args.pbradii:
        PBRadii = known_pbradii[args.igb]
    else:
        PBRadii = None

    ambertools.prepare_leap_config_file('leap.in', 'protein.pdb', mol2files_l, 'complex.pdb', solvate=solvate, box=args.box, distance=args.boxsize, model=args.water, version=amber_version, PBRadii=PBRadii, membrane=args.membrane)
    utils.run_shell_command('tleap -f leap.in')

    if args.addions != 0.0 and args.solvent == 'explicit':
        nna, ncl = ambertools.get_ions_number('leap.log', concentration=args.addions)
        ambertools.prepare_leap_config_file('leap.in', 'protein.pdb', mol2files_l, 'complex.pdb', solvate=solvate, box=args.box, distance=args.boxsize, nna=nna, ncl=ncl, model=args.water, version=amber_version, PBRadii=PBRadii)
        utils.run_shell_command('tleap -f leap.in')

    if args.namd:
        ligname = []
        if args.file_l:
            for file_l in mol2files_l:
                ligname.append(ambertools.get_ligand_name(file_l))
        namdtools.create_constrained_pdbfile('namd_equil_res.pdb', 'start.pdb', ligname)

    os.chdir(pwd)

if args.partition == 'gpu':
    script = """#!/bin/bash
#SBATCH --time=100-00:00
#SBATCH --partition=%(partition)s
#SBATCH --job-name="md"
#SBATCH --cpus-per-task=%(ncpus)s

set -e
\n"""%locals()
else:
    script = """#!/bin/bash
#SBATCH --time=100-00:00
#SBATCH --partition=%(partition)s
#SBATCH --job-name="md"
#SBATCH --ntasks=%(ncpus)s
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

set -e\n"""%locals()

if args.partition == 'gpu':
    if amber_version in ['14', '15']:
        script += """export AMBERHOME=/nfs/r510-2/pwinter/PharmaApps/build/amber/gpu/amber14
export CUDA_HOME=/usr/local/cuda-9.1
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
source $AMBERHOME/amber.sh
export CUDA_VISIBLE_DEVICES=`/usr/bin/select_gpu.py`
\n"""%locals()
    elif amber_version in ['16', '17']:
        script += """export AMBERHOME=/nfs/r510-2/pwinter/PharmaApps/build/amber/amber16/gpu/amber16
export CUDA_HOME=/usr/local/cuda-9.1
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
source $AMBERHOME/amber.sh
export CUDA_VISIBLE_DEVICES=`/usr/bin/select_gpu.py`
\n"""%locals()
    exe = 'pmemd.cuda'
else:
    if args.ncpus == 1:
      exe = 'sander'
    elif args.ncpus > 1:
      exe = 'mpirun -np %(ncpus)s sander.MPI'%locals()
    else:
        raise ValueError('Number of CPUS (-np) should be greater or equal to 1')

script += "cd %(workdir_r)s\n"""%locals()

# ------- RESTRAINTS ---------
if args.restraints:
    if args.restraints == "backbone":
        restraints_lines = """\nntr=1, restraint_wt=5.0,
restraintmask='@CA,C,N,O&!:WAT',"""
    else:
        restraints_lines = args.restraints
else:
    restraints_lines = ""

# ------- STEP 2: minimization -------
if 'min' in args.step:
    workdir = workdir_r + '/min'
    shutil.rmtree(workdir, ignore_errors=True) 
    os.mkdir(workdir)

    if args.solvent == 'explicit':
        solvent_mask = ambertools.get_solvent_mask(workdir_r+'/common/start.pdb', residues='WAT,Na+,Cl-')
 
        if args.restraints:
            ref_flag = ' -ref min1.rst'
        else:
            ref_flag = ''

        with open(workdir+'/min1.in', 'w') as min1f:
            min1f.write("""Minimization 1 solvent + ions
&cntrl
imin=1,
maxcyc=%(maxcyc_s)s,
ncyc=%(ncyc_s)s,
ntpr=100,
ntb=1,
cut=%(cut)s,
ntr=1,restraint_wt=100.0,
restraintmask='!%(solvent_mask)s&!@H='
/\n"""%locals())

        with open(workdir+'/min2.in', 'w') as min2f:
            min2f.write("""Minimization 2 (full structure)
&cntrl
imin=1,
maxcyc=%(maxcyc)s,
ncyc=%(ncyc)s,
ntpr=100,
ntb=1,
cut=%(cut)s,
/\n"""%locals())

        script += """\ncd min
# minimization with restrains
%(exe)s -O -i min1.in -p ../common/start.prmtop -c ../common/start.inpcrd -o min1.out -r min1.rst -ref ../common/start.inpcrd

# minimization (full structure)
%(exe)s -O -i min2.in -p ../common/start.prmtop -c min1.rst -o min2.out -r min2.rst%(ref_flag)s

# convert final structure to PDB
cpptraj -p ../common/start.prmtop -y min2.rst -x min2.pdb
cd ..\n"""%locals()

    elif args.solvent in ['vacuo', 'implicit']:

        # input files follow examples at: http://ambermd.org/tutorials/basic/tutorial1/section3.htm and
        # http://ambermd.org/tutorials/basic/tutorial1/section4.htm
        if args.solvent == 'vacuo':
            title = 'In-Vacuo minimization'
            igb_line = '\nigb=0,'
        elif args.solvent == 'implicit':
            title = 'Minimization (implicit solvent)'
            igb_line = '\nigb=%i,'%args.igb

        with open(workdir+'/min.in', 'w') as minf:
            minf.write("""%(title)s
&cntrl
imin=1,
maxcyc=%(maxcyc)s,
ncyc=%(ncyc)s,
ntb=0,%(igb_line)s
cut=%(cut)s,
/\n"""%locals())

        script += """\ncd min
# minimization
%(exe)s -O -i min.in -p ../common/start.prmtop -c ../common/start.inpcrd -o min.out -r min.rst -ref ../common/start.inpcrd

# convert final structure to PDB
cpptraj -p ../common/start.prmtop -y min.rst -x min.pdb
cd ..\n"""%locals()


# ------- STEP 3: heating -------
if 'nvt' in args.step and args.solvent == 'explicit':
    workdir = workdir_r + '/nvt'

    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    if args.restraints:
        ref_flag = ' -ref ../min/min2.rst'
    else:
        ref_flag = ''

    with open(workdir+'/nvt.mdin', 'w') as nvtf:
        nvtf.write("""Equilibration in NVT ensemble
&cntrl
nstlim=%(nstlim_h)s,
dt=0.002,
ntpr=1000,
ig=-1,
tempi=0.0, temp0=%(temp)s,
ntt=1,
tautp=2.0,
ntb=1, ntp=0,
cut=%(cut)s,
ntc=2, ntf=2,%(restraints_lines)s
/\n"""%locals())

    script += """\ncd nvt
# run heating
%(exe)s -O -i nvt.mdin -o nvt.mdout -c ../min/min2.rst -r nvt.rst -x nvt.mdcrd -inf nvt.mdinfo -p ../common/start.prmtop%(ref_flag)s
cd ..\n"""%locals()

# ------- STEP 4: equilibration -------
if 'npt' in args.step and args.solvent == 'explicit':
    workdir = workdir_r + '/npt'

    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    if args.restraints:
        ref_flag = ' -ref ../nvt/nvt.rst'
    else:
        ref_flag = ''

    with open(workdir+'/npt.mdin', 'w') as nptf:
        nptf.write("""Equilibration in NPT ensemble
&cntrl
nstlim=%(nstlim_e)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0, tautp=2.0,
ntb=2, pres0=1.0, ntp=1, taup=2.0,
cut=%(cut)s,
irest=1, ntx=5,
ntc=2, ntf=2,
ntpr=500,%(restraints_lines)s
/\n"""%locals())

    script += """\ncd npt
# run equilibration
%(exe)s -O -i npt.mdin -o npt.mdout -c ../nvt/nvt.rst -r npt.rst -x npt.mdcrd -inf npt.mdinfo -p ../common/start.prmtop%(ref_flag)s
cd ..\n"""%locals()


# ------- STEP 5: production -------
if 'md' in args.step:
    mddir = 'md'
    workdir = workdir_r + '/' + mddir

    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    if args.restraints and solvent == 'explicit':
        ref_flag = " -ref ../npt/npt.rst"
    else:
        ref_flag = ""

    if solvent == 'explicit':
        ntwprt_lines = ""
        if args.ntwprt == 1:
            with open(workdir_r+'/common/start.pdb') as pdbf:
                for line in pdbf:
                    if line.startswith(('ATOM', 'HETATM')):
                        resname = line[17:20].strip()
                        if resname == 'WAT':
                            last_protein_atom = int(line[6:11].strip())-1
                            break
            ntwprt_lines = " ntwprt=%s,"%last_protein_atom
 
        with open(workdir+'/md.mdin', 'w') as mdf:
            mdf.write("""MD production
&cntrl
nstlim=%(nstlim)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0, tautp=2.0,
ntb=2, pres0=1.0, ntp=1, taup=2.0,
cut=%(cut)s,
irest=1, ntx=5,
iwrap=%(iwrap)s,
ntc=2, ntf=2,
ntpr=%(ntpr)s, ntwr=%(ntpr)s, ntwx=%(ntwx)s,%(ntwprt_lines)s%(restraints_lines)s
/\n"""%locals())

        script += """\ncd %(mddir)s
# production run
%(exe)s -O -i md.mdin -o md.mdout -c ../npt/npt.rst -r md.rst -x md.mdcrd -inf md.mdinfo -p ../common/start.prmtop%(ref_flag)s
cd ..\n"""%locals()

    else: # implicit solvent
        if args.solvent == 'vacuo':
            title = 'MD Production (In-Vacuo)'
            igb_line = '\nigb=0,'
        elif args.solvent == 'implicit':
            title = 'MD Production (implicit solvent)'
            igb_line = '\nigb=%i,'%args.igb

        with open(workdir+'/md.mdin', 'w') as mdf:
            mdf.write("""%(title)s
&cntrl
nstlim=%(nstlim)s,
dt=0.002,
ig=-1,
tempi=%(temp)s, temp0=%(temp)s,
ntt=3, gamma_ln=1.0,
ntb=0,%(igb_line)s
cut=%(cut)s,
ntpr=%(ntpr)s, ntwr=%(ntpr)s, ntwx=%(ntwx)s
/\n"""%locals())

        script += """\ncd %(mddir)s
# production run
%(exe)s -O -i md.mdin -o md.mdout -c ../min/min.rst -r md.rst -x md.mdcrd -inf md.mdinfo -p ../common/start.prmtop%(ref_flag)s
cd ..\n"""%locals()

if args.step != ['prep']:
    with open(args.script_name, 'w') as ff:
        ff.write(script)
