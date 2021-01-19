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
    dest='files_l',
    required=False,
    default=None,
    nargs='+',
    help="Input files for ligand (.pdb, .mol2)")

parser.add_argument('-addions',
    dest='addions',
    type=float,
    default=0.0,
    help="Add specific concentration of ions (Na+, Cl-) and neutralize the system")

parser.add_argument('-amd',
    dest='amd',
    type=str,
    default=None,
    help="Accelerated MD options (default: accelerated MD disabled)")

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

parser.add_argument('-efn',
    dest='efn',
    default=0,
    type=int,
    help="Normalize electic field by box length in Z direction (voltage in kcal/(mol*e))")

parser.add_argument('-efz',
    dest='efz',
    default=None,
    help="Electric field in Z direction (in kcal/(mol*A*e))")

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

parser.add_argument('-irest',
    dest='irest',
    type=int,
    default=1,
    help="Production: irest (default: 1)")

parser.add_argument('-iwrap',
    dest='iwrap',
    type=int,
    default=0,
    help="Value of iwrap for production run (default: 1, i.e., the coordinates written to the restart and trajectory files are not wrapped into a primary box)")

parser.add_argument('-membrane',
    dest='membrane',
    action='store_true',
    default=False,
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
    default=10000,
    help="Minimization (solvent): maxcyc (default: 10000)")

parser.add_argument('-maxcyc',
    dest='maxcyc',
    default=10000,
    help="Minimization (full): maxcyc (default: 10000)")

parser.add_argument('-ncyc_s',
    dest='ncyc_s',
    default=5000,
    help="Minimization (solvent): ncyc (default: 5000)")

parser.add_argument('-ncyc',
    dest='ncyc',
    default=5000,
    help="Minimization (full): ncyc (default: 5000)")

parser.add_argument('-nstlim_nvt',
    dest='nstlim_nvt',
    default=250000,
    help="NVT: nstlim (default: 250000)")

parser.add_argument('-nstlim_npt',
    dest='nstlim_npt',
    default=250000,
    help="NPT: nstlim (default: 250000)")

parser.add_argument('-nstlim',
    dest='nstlim',
    default=500000,
    help="Production: nstlim (default: 500000)")

parser.add_argument('-ntpr',
    dest='ntpr',
    type=int,
    default=50,
    help="Production: ntpr and ntwr for restart files (default: 50)")

parser.add_argument('-ntwprt',
    dest='ntwprt',
    type=int,
    choices=[0, 1],
    default=0,
    help="Atoms to include in trajectory file, 0: all atoms, 1: no solvent (default: 0)")

parser.add_argument('-ntwx',
    dest='ntwx',
    type=int,
    default=0,
    help="Production: ntwx (default: 0)")

parser.add_argument('-nwaters',
    dest='nwaters',
    type=int,
    default=None,
    help="Number of waters (to be used with targeted MD)")

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

parser.add_argument('-ref',
    dest='reffile',
    type=str,
    default=None,
    help="Specify reference structure (to be used with targeted MD)")

parser.add_argument('-res',
    dest='resfile',
    type=str,
    default=None,
    help="Specify restraint file")

parser.add_argument('-rst',
    dest='rstfile',
    required=False,
    help=".rst file used when restarting MD (works only when -st md)")

parser.add_argument('-rsttop',
    dest='rsttop',
    required=False,
    help=".prmtop file used when restarting MD (works only when -st md)")

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

parser.add_argument('-tgtmd',
    dest='tgtmd',
    type=str,
    default=None,
    help="Targeted MD options (default: targeted MD disabled)")

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

# check electric field option
if not amber_version in ['16', '17'] and efz is not None:
    sys.exit('Electric field option only supported for Amber 16 or 17')
elif efz:
    ef_lines = "\nefz=%s,"%efz
    if efn == 1:
        ef_lines += " efn=1,"
else:
    ef_lines = "" 

## get accelerated MD lines
amd_lines = utils.extract_amd_lines(amd)

## get targeted MD lines
tgtmd_lines = utils.extract_tgtmd_lines(tgtmd)

if rstfile is not None and step != ['md']:
    sys.exit("rst flag can only be used with -st md")
elif step == ['md']:
    pass

## save current directory
pwd = os.getcwd()

if file_r is not None:
    # get absolute paths of ligand and receptor files
    file_r_abspath = os.path.abspath(file_r)
elif 'prep' in step:
    raise ValueError('Receptor file should be provided!')

if files_l is not None:
    files_l_abspath = [os.path.abspath(filename) for filename in files_l]

# ----- SET WORKING DIRECTORY ------
if workdir is not None:
    workdir_abspath = os.path.abspath(workdir)
else:
    dir_r = os.path.dirname(file_r_abspath)
    file_r_prefix, ext = os.path.splitext(os.path.basename(file_r_abspath))
    workdir_abspath = dir_r + '/' + file_r_prefix

# ------- STEP 1: preparation -------
if 'prep' in step:
    shutil.rmtree(workdir_abspath, ignore_errors=True)
    os.mkdir(workdir_abspath)

    workdir_curr = workdir_abspath + '/common'
    shutil.rmtree(workdir_curr, ignore_errors=True)
    os.makedirs(workdir_curr)
    os.chdir(workdir_curr)

    ambertools.prepare_receptor('protein.pdb', file_r_abspath, keep_hydrogens=keeph, membrane=membrane)
    # ligand preparation
    if files_l is not None:
        files_l_new = []
        # preparing ligand(s)...
        for filename in files_l_abspath:
            name = ambertools.get_ligand_name(filename)

            prefix, ext = os.path.splitext(filename)
            filename_new = 'ligand_%s%s'%(name, ext)
            shutil.copyfile(filename, filename_new)
            files_l_new.append(filename_new)

        mol2files_l = ambertools.prepare_ligand('protein.pdb', files_l_new, 'complex.pdb', charge_method=charge_method, \
        version=amber_version, skip_unrecognized_atoms=hem)
    else:
        mol2files_l = None

    # should solvate??
    if solvent == 'explicit' and not membrane:
        solvate = True
    else:
        solvate = False

    # get pdbradii
    if pbradii:
        PBRadii = known_pbradii[igb]
    else:
        PBRadii = None

    if not nwaters or solvent in ['implicit', 'vacuo']:
        ambertools.prepare_leap_config_file('leap.in', 'protein.pdb', mol2files_l, 'complex.pdb', solvate=solvate, \
        box=box, distance=boxsize, model=water, version=amber_version, PBRadii=PBRadii, membrane=membrane)
        utils.run_shell_command('tleap -f leap.in')
    else:
        removed_waters, dbest, cbest = ambertools.get_removed_waters('protein.pdb', mol2files_l, 'complex.pdb', nwaters, boxsize, step=0.01, ntries=5, version=amber_version)
        ambertools.prepare_leap_config_file('leap.in', 'protein.pdb', mol2files_l, 'complex.pdb', solvate=solvate, \
        box=box, distance=dbest, closeness=cbest, removed_waters=removed_waters, model=water, version=amber_version, PBRadii=PBRadii, membrane=membrane)
        utils.run_shell_command('tleap -f leap.in')

    if addions != 0.0 and solvate:
        nna, ncl = ambertools.get_ions_number('leap.log', concentration=addions, version=amber_version)
        ambertools.prepare_leap_config_file('leap.in', 'protein.pdb', mol2files_l, 'complex.pdb', solvate=solvate, \
        box=box, distance=boxsize, nna=nna, ncl=ncl, model=water, version=amber_version, PBRadii=PBRadii)
        utils.run_shell_command('tleap -f leap.in')

    if namd:
        ligname = []
        for file_l in mol2files_l:
            ligname.append(ambertools.get_ligand_name(file_l))
        namdtools.create_constrained_pdbfile('namd_equil_res.pdb', 'start.pdb', ligname)
    os.chdir(pwd)

script = """#!/bin/bash
#$ -N md_298K_IDP
#$ -q r730gpuRTX2080ti
#$ -S /bin/bash
#$ -V 
#$ -cwd 
set -e
\n"""%locals()

if partition == 'gpu':
    exe = 'pmemd.cuda'
    cpptrajexe = 'cpptraj.cuda'
else:
    if ncpus == 1:
      exe = 'sander'
    elif ncpus > 1:
      exe = 'mpirun -np %(ncpus)s sander.MPI'%locals()
    else:
        raise ValueError('Number of CPUS (-np) should be greater or equal to 1')
    cpptrajexe = 'cpptraj'
script += "cd %s\n"""%os.path.relpath(workdir_abspath)

# ------- STEP 2: minimization -------
if 'min' in step:
    workdir_curr = workdir_abspath + '/min'
    shutil.rmtree(workdir_curr, ignore_errors=True) 
    os.mkdir(workdir_curr)

    if solvent == 'explicit':
        solvent_mask = ambertools.get_solvent_mask(workdir_abspath+'/common/start.pdb', residues='WAT,Na+,K+,Cl-')

        with open(workdir_curr+'/min1.in', 'w') as min1f:
            min1f.write("""Minimization 1 solvent + ions
&cntrl
imin=1,
maxcyc=%(maxcyc_s)s,
ncyc=%(ncyc_s)s,
ntpr=100,
ntp=0, ntb=1,
cut=%(cut)s,
ntr=1,restraint_wt=10.0,
restraintmask='!%(solvent_mask)s&!@H='
/\n"""%locals())

        with open(workdir_curr+'/min2.in', 'w') as min2f:
            min2f.write("""Minimization 2 (full structure)
&cntrl
imin=1,
maxcyc=%(maxcyc)s,
ncyc=%(ncyc)s,
ntpr=100,
ntp=0, ntb=1,
cut=%(cut)s,
/\n"""%locals())

        script += """\ncd min
# minimization with restraints
%(exe)s -O -i min1.in -p ../common/start.prmtop -c ../common/start.inpcrd -o min1.out -r min1.rst -ref ../common/start.inpcrd

# minimization (full structure)
%(exe)s -O -i min2.in -p ../common/start.prmtop -c min1.rst -o min2.out -r min2.rst
cd ..\n"""%locals()

    elif solvent in ['vacuo', 'implicit']:
        # source: http://ambermd.org/tutorials/basic/tutorial1/section3.htm
        if solvent == 'vacuo':
            title = 'In-Vacuo minimization'
            igb_line = '\nigb=0,'
        elif solvent == 'implicit':
            title = 'Minimization (implicit solvent)'
            igb_line = '\nigb=%i,'%igb

        with open(workdir_curr+'/min.in', 'w') as minf:
            minf.write("""%(title)s
&cntrl
imin=1,
maxcyc=%(maxcyc)s,
ncyc=%(ncyc)s,
ntp=0, ntb=0,%(igb_line)s
cut=%(cut)s,
/\n"""%locals())

        script += """\ncd min
# minimization
%(exe)s -O -i min.in -p ../common/start.prmtop -c ../common/start.inpcrd -o min.out -r min.rst -ref ../common/start.inpcrd

# convert final structure to PDB
%(cpptrajexe)s -p ../common/start.prmtop -y min.rst -x min.pdb
cd ..\n"""%locals()


# ------- STEP 3: heating -------
if 'nvt' in step and solvent == 'explicit':
    workdir_curr = workdir_abspath + '/nvt'

    shutil.rmtree(workdir_curr, ignore_errors=True)
    os.mkdir(workdir_curr)

    if membrane:
        # source: http://ambermd.org/tutorials/advanced/tutorial16
        with open(workdir_curr+'/nvt1.mdin', 'w') as nvtf:
            nvtf.write("""Heating in NVT ensemble
&cntrl
nstlim=2500,
dt=0.002,
ntpr=1000,
ig=-1,
ntt=3, gamma_ln=1.0, 
ntp=0, ntb=1,
cut=%(cut)s,
ntc=2, ntf=2
ntr=1,restraint_wt=10.0,
nmropt=1,
restraintmask='!%(solvent_mask)s&!@H='
/
&wt type='TEMP0', istep1=0, istep2=2500, value1=0.0, value2=100.0
/
&wt type='END'
/\n"""%locals())

        with open(workdir_curr+'/nvt2.mdin', 'w') as nvtf:
            nvtf.write("""Heating in NPT ensemble (anisotropic barostat)
&cntrl
nstlim=%(nstlim_nvt)s,
dt=0.002,
ntpr=1000,
ig=-1,
ntt=3, gamma_ln=1.0, 
ntb=2, ntp=2,
pres0=1.0, taup=2.0,
cut=%(cut)s,
irest=1, ntx=5,
ntc=2, ntf=2
ntr=1,restraint_wt=10.0,
nmropt=1,
restraintmask='!%(solvent_mask)s&!@H='
/
&wt type='TEMP0', istep1=0, istep2=%(nstlim_nvt)s, value1=100.0, value2=%(temp)s
/
&wt type='END'
/\n"""%locals())

        script += """\ncd nvt
# run heating 1
%(exe)s -O -i nvt1.mdin -o nvt1.mdout -c ../min/min2.rst -r nvt1.rst -x nvt1.mdcrd -inf nvt1.mdinfo -p ../common/start.prmtop -ref ../min/min2.rst
# run heating 2
%(exe)s -O -i nvt2.mdin -o nvt2.mdout -c nvt1.rst -r nvt2.rst -x nvt2.mdcrd -inf nvt2.mdinfo -p ../common/start.prmtop -ref nvt1.rst
cd ..\n"""%locals()

    else:
        with open(workdir_curr+'/nvt.mdin', 'w') as nvtf:
            nvtf.write("""Equilibration in NVT ensemble
&cntrl
nstlim=%(nstlim_nvt)s,
dt=0.002,
ntpr=1000,
ig=-1,
ntt=1, tautp=2.0, 
ntp=0, ntb=1,
cut=%(cut)s,
ntc=2, ntf=2
nmropt=1,
restraintmask='!%(solvent_mask)s&!@H='
/
&wt type='TEMP0', istep1=0, istep2=%(nstlim_nvt)s, value1=0.0, value2=%(temp)s
/
&wt type='END'
/\n"""%locals())

        script += """\ncd nvt
# run heating
%(exe)s -O -i nvt.mdin -o nvt.mdout -c ../min/min2.rst -r nvt.rst -x nvt.mdcrd -inf nvt.mdinfo -p ../common/start.prmtop -ref ../min/min2.rst
cd ..\n"""%locals()

# ------- STEP 4: equilibration -------
if 'npt' in step and solvent == 'explicit':
    workdir_curr = workdir_abspath + '/npt'

    shutil.rmtree(workdir_curr, ignore_errors=True)
    os.mkdir(workdir_curr)

    if membrane:
        # source: http://ambermd.org/tutorials/advanced/tutorial16
        with open(workdir_curr+'/npt.mdin', 'w') as nptf:
            nptf.write("""Equilibration in NPT ensemble
&cntrl
nstlim=%(nstlim_npt)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0,
ntp=2, ntb=2,
pres0=1.0, taup=2.0,
cut=%(cut)s,
irest=1, ntx=5,
ntc=2, ntf=2,
ntpr=500
/
&ewald
skinnb=5, ! Increase skinnb to avoid skinnb errors
/\n"""%locals())

        script += "\ncd npt\n# run equilibration"
        for idx in range(10):
            if idx == 0:
                script += """\n%s -O -i npt.mdin -o npt%i.mdout -c ../nvt/nvt2.rst -r npt%i.rst -x npt%i.mdcrd -inf npt%i.mdinfo \
-p ../common/start.prmtop"""%(exe,idx+1,idx+1,idx+1,idx+1)
            else:
                script += """\n%s -O -i npt.mdin -o npt%i.mdout -c npt%i.rst -r npt%i.rst -x npt%i.mdcrd -inf npt%i.mdinfo \
-p ../common/start.prmtop"""%(exe,idx+1,idx,idx+1,idx+1,idx+1)
        script += "\ncd ..\n"

    else:
        with open(workdir_curr+'/npt.mdin', 'w') as nptf:
            nptf.write("""Equilibration in NPT ensemble
&cntrl
nstlim=%(nstlim_npt)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0,
ntp=1, ntb=2,
pres0=1.0, taup=2.0,
cut=%(cut)s,
irest=1, ntx=5,
ntc=2, ntf=2,
ntpr=500
/
&ewald
skinnb=5, ! Increase skinnb to avoid skinnb errors
/\n"""%locals())

        script += """\ncd npt
# run equilibration
%(exe)s -O -i npt.mdin -o npt.mdout -c ../nvt/nvt.rst -r npt.rst -x npt.mdcrd -inf npt.mdinfo -p ../common/start.prmtop
cd ..\n"""%locals()


# ------- STEP 5: production -------
if 'md' in step:
    if rstfile is None:
        mddir = 'md'
        workdir_curr = workdir_abspath + '/' + mddir

        shutil.rmtree(workdir_curr, ignore_errors=True)
        os.mkdir(workdir_curr)

        startpdb = workdir_abspath + '/common/start.pdb'
        prmtop = os.path.relpath(workdir_abspath+'/common/start.prmtop', workdir_curr) 
    else:
        rstdir = os.path.dirname(rstfile)
        if rstdir:
            rstdir = os.path.relpath(os.path.dirname(rstfile))
        rstdir_s = rstdir.split('/')
        if len(rstdir_s) > 1:
            workdir_rst = '/'.join(rstdir_s[:-1])
        else:
            workdir_rst = '.'
        if workdir_rst == os.path.relpath(workdir_abspath):
            idx = 1
            mddir = 'md'
            while os.path.isdir(workdir_abspath+'/'+mddir):
                if idx == 1:
                    mddir = 'md_ext'
                else:
                    mddir = 'md_ext_' + str(idx)
                idx += 1
        else:
            mddir = 'md'
        workdir_curr = workdir_abspath + '/' + mddir
        shutil.rmtree(workdir_curr, ignore_errors=True)
        os.makedirs(workdir_curr)

        if rsttop is None:
            startpdb = workdir_rst + '/common/start.pdb'
            prmtop = workdir_rst+'/common/start.prmtop'
        else:
            topdir = os.path.relpath(os.path.dirname(rsttop))
            topdir_s = topdir.split('/')
            if len(topdir_s) > 1:
                topdir_rst = '/'.join(topdir_s[:-1])
            else:
                topdir_rst = '.'
            startpdb = topdir_rst + '/common/start.pdb'
            prmtop = topdir_rst + '/common/start.prmtop'
        prmtop = os.path.relpath(prmtop, workdir_curr)

    if resfile is not None:
        ref_flag = ' -ref %s'%(os.path.relpath(resfile, workdir_curr))
    elif reffile is not None:
        ref_flag = ' -ref %s'%(os.path.relpath(reffile, workdir_curr))
    else:
        ref_flag = ""

    if solvent == 'explicit':
        # get ntwprt string from option value
        if ntwprt == 0:
            ntwprt_str = ""
        else:
            last_protein_atom = ambertools.get_last_solute_atom_num(startpdb)
            ntwprt_str = "\nntwprt=%s,"%last_protein_atom

        # fix barostat options if electric field is enabled
        if efz is None and membrane:
            barostat = "\nntp=2, ntb=2,\npres0=1.0, taup=2.0,"
        elif efz is None and not membrane:
            barostat = "\nntp=1, ntb=2,\npres0=1.0, taup=2.0,"
        else:
            barostat = "\nntp=0, ntb=1,"

        if irest == 0:
            rstline = "irest=0, ntx=1, tempi=0.0,"
        else:
            rstline = "irest=1, ntx=5"

        if resfile is not None:
            resline = "\nntr=1, restraint_wt=10.0,\nrestraintmask='!:842-13396&!@H='"
        else:
            resline = ""

        if membrane:
            if rstfile is not None:
                rstfile = os.path.relpath(rstfile, workdir_curr)
            else:
                rstfile = "../npt/npt10.rst"

            # source: http://ambermd.org/tutorials/advanced/tutorial16
            with open(workdir_curr+'/md.mdin', 'w') as mdf:
                mdf.write("""MD production
&cntrl
nstlim=%(nstlim)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0,%(barostat)s
cut=%(cut)s,
%(rstline)s,
iwrap=%(iwrap)s,
ntc=2, ntf=2,
ntpr=%(ntpr)s, ntwr=%(ntpr)s, ntwx=%(ntwx)s,%(ntwprt_str)s%(ef_lines)s%(amd_lines)s%(tgtmd_lines)s%(resline)s
/\n"""%locals())

            script += """\ncd %(mddir)s
# run equilibration
%(exe)s -O -i md.mdin -o md.mdout -c %(rstfile)s -r md.rst -x md.mdcrd -inf md.mdinfo -p %(prmtop)s%(ref_flag)s
cd ..\n"""%locals()

        else:
            if rstfile is not None:
                rstfile = os.path.relpath(rstfile, workdir_curr)
            else:
                rstfile = "../npt/npt.rst"
            with open(workdir_curr+'/md.mdin', 'w') as mdf:
                mdf.write("""MD production
&cntrl
nstlim=%(nstlim)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0,%(barostat)s
cut=%(cut)s,
%(rstline)s,
iwrap=%(iwrap)s,
ntc=2, ntf=2,
ntpr=%(ntpr)s, ntwr=%(ntpr)s, ntwx=%(ntwx)s,%(ntwprt_str)s%(ef_lines)s%(amd_lines)s%(tgtmd_lines)s
/\n"""%locals())

            script += """\ncd %(mddir)s
# production run
%(exe)s -O -i md.mdin -o md.mdout -c %(rstfile)s -r md.rst -x md.mdcrd -inf md.mdinfo -p %(prmtop)s%(ref_flag)s
cd ..\n"""%locals()

    else: # implicit solvent
        if solvent == 'vacuo':
            title = 'MD Production (In-Vacuo)'
            igb_line = '\nigb=0,'
        elif solvent == 'implicit':
            title = 'MD Production (implicit solvent)'
            igb_line = '\nigb=%i,'%igb

        if rstfile is not None:
            rstfile = os.path.relpath(rstfile, workdir_curr)
        else:
            rstfile = "../min/min.rst"
        with open(workdir_curr+'/md.mdin', 'w') as mdf:
            mdf.write("""%(title)s
&cntrl
nstlim=%(nstlim)s,
dt=0.002,
ig=-1,
temp0=%(temp)s,
ntt=3, gamma_ln=1.0,
ntb=0,%(igb_line)s
cut=%(cut)s,
ntpr=%(ntpr)s, ntwr=%(ntpr)s, ntwx=%(ntwx)s
/\n"""%locals())

        script += """\ncd %(mddir)s
# production run
%(exe)s -O -i md.mdin -o md.mdout -c %(rstfile)s -r md.rst -x md.mdcrd -inf md.mdinfo -p ../common/start.prmtop%(ref_flag)s
cd ..\n"""%locals()

if step != ['prep']:
    with open(script_name, 'w') as ff:
        ff.write(script)
