#!/usr/bin/env python
import os
import sys
import argparse
import shutil
from glob import glob

parser = argparse.ArgumentParser(description="Restart MD simulations with Amber")

parser.add_argument('-gpu',
    action='store_true',
    default=False,
    help="Run simulations on gpu")

parser.add_argument('-np',
    dest='ncpus',
    type=int,
    default=1,
    help="Number of cpus used for the simulations (default: 1)")

parser.add_argument('-w',
    dest='workdir',
    required=True,
    help="name of working directory")

args = parser.parse_args()

if args.gpu:
    exe = 'pmemd.cuda'
elif args.ncpus == 1:
    exe = 'sander'
elif args.ncpus > 1:
    exe = 'mpirun -np %(ncpus)s sander.MPI'%locals()
else:
    raise ValueError('Number of CPUS (-np) should be greater or equal to 1')

# script used to restart AMBER md upon failure 
# initial MD simulations are supposed to be prepared with the prepare_md.py script

# assume the correct MD directory is given
if os.path.isdir(args.workdir):
    commondir = args.workdir + '/common'
    if os.path.isdir(commondir):
        # check md directories
        if glob(args.workdir+'/md_*'):
            index_md = [int(os.path.basename(dir)[3:]) for dir in glob(args.workdir+'/md_*')]
            index_md_max = max(index_md)
            last_mddir = args.workdir + '/md_%i'%index_md_max
            new_mddir = args.workdir + '/md_%i'%(index_md_max+1)
            new_index = index_md_max + 1
        elif os.path.isdir(args.workdir+'/md'):
            last_mddir = args.workdir + '/md'
            new_mddir = args.workdir + '/md_1'
            new_index = 1
        else:
            raise ValueError("No MD folder found!")
    else:
        raise ValueError("Directory %s does not exist!"%commondir)
else:
    raise ValueError("Directory %s does not exist!"%args.workdir)

ntpr = 0
nstlim = 0
# get number of initial steps
with open(last_mddir+'/md.mdin', 'r') as mdin:
    for line in mdin:
        line_no_comment = line.split('#')[0]
        line_s = line_no_comment.split(',')
        for option in line_s:
            if option.strip().startswith('ntpr'):
                ntpr = int(option.split('=')[1])
            elif option.strip().startswith('nstlim'):
                nstlim = int(option.split('=')[1])

restart_last_md = True
mdinfo = last_mddir+'/md.mdinfo'
if os.path.isfile(mdinfo):
    with open(mdinfo, 'r') as infof:
        for line in infof:
            # check number of steps run in the last simulation
            if 'Total steps' in line and 'Completed' in line and 'Remaining' in line:
                restart_last_md = False
                nstlim_info = int(line[15:25])
                ncompleted_info = int(line[39:49])
                nremaining_info = int(line[63:73])
                assert nstlim == nstlim_info
                restart_last_md = False
                if nremaining_info <= ntpr:
                    sys.exit('Dynamics already done! No need to extend it!')

if restart_last_md:
    if new_index <= 1:
        raise ValueError("MD info not found (1st run)!")
    else:
        print "No MD info found, restarting dynamics from last run"
    if new_index == 2:
        new_index = new_index - 1
        new_mddir = last_mddir
        last_mddir = args.workdir + '/md'
    else:
        new_index = new_index - 1
        new_mddir = last_mddir
        last_mddir = args.workdir + '/md_%i'%(new_index-1)

    # purging unused files
    for ff in glob(last_mddir+'/*'):
        if os.path.basename(ff) != 'md.mdin':
            shutil.rmtree(ff, ignore_errors=True)
else:
    # create new md dir
    os.mkdir(new_mddir)
    # write new mdin file (use the read me of the old)
    shutil.copyfile(last_mddir+'/md.mdin', new_mddir+'/md.mdin')

    with open(new_mddir+'/md.mdin.tmp', 'w') as tmpf:
        # get number of initial steps
        with open(new_mddir+'/md.mdin', 'r') as mdin:
            for line in mdin:
                line_no_comment = line.split('#')[0]
                if 'nstlim' in line_no_comment:
                    new_line = ''
                    line_s = line_no_comment.split(',')
                    for option in line_s:
                        if option.strip().startswith('nstlim'):
                            new_option = 'nstlim=%i'%nremaining_info
                        else:
                            new_option = option
                        if new_line:
                            new_line += ','
                        new_line += new_option
                else:
                    new_line = line
                tmpf.write(new_line)
    shutil.move(new_mddir+'/md.mdin.tmp', new_mddir+'/md.mdin')

last_rstfile = last_mddir+'/md.rst'
if os.path.isfile(last_rstfile):
     last_rstfile_rel = os.path.relpath(last_mddir+'/md.rst', new_mddir)
new_mddir_abs = os.path.abspath(new_mddir)

with open('restart_md_%s.sh'%new_index, 'w') as shf:
   contents ="""#!/bin/bash
cd %(new_mddir_abs)s

# production run
%(exe)s -O -i md.mdin -o md.mdout -c %(last_rstfile_rel)s -r md.rst -x md.mdcrd -inf md.mdinfo -p ../common/start.prmtop
cd ..\n"""%locals()
   shf.write(contents)
