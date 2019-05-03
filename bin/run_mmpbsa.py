#!/usr/bin/env python
from __future__ import with_statement

import sys
import os
import subprocess
import argparse

class MMPBSAConfigError(Exception):
    pass

class MMPBSAWorker(object):

    def prepare_pbsa_input_file(self, startframe=1, endframe=1000):

        # write mmpbsa input file
        with open('mm.in', 'w') as file:
            script ="""Input file for running PB
&general
  startframe=%(startframe)s, endframe=%(endframe)s, keep_files=0, verbose=2,
/
&pb 
  istrng=0.0, inp=1, cavity_surften=0.0072, cavity_offset=0, 
  prbrad=1.6,
/   
&gb
  igb=5, saltcon=0.0, surften=0.005, surfoff=0,
/"""% locals()
            file.write(script)

    def prepare_entropy_input_file(self, startframe=1, endframe=1000):

        # write mmpbsa with entropy input file
        with open('mm_entropy.in', 'w') as file:
            script ="""Input file for running PB with entropy
&general
  startframe=%(startframe)s, endframe=%(endframe)s, keep_files=0, verbose=2,
/
&nmode
  nmstartframe=1, nminterval=1, nmode_igb=1, nmode_istrng=0.0,
/"""% locals()
            file.write(script)

    def create_arg_parser(self):

        parser = argparse.ArgumentParser(description='Run MMPBSA simulations..')

        parser.add_argument('-s',
            dest='solvent_mask',
            type=str,
            default=':WAT,Na+,Cl-',
            help='Mask used for receptor!')

        parser.add_argument('-n',
            dest='ligand_mask',
            type=str,
            default=':LIG',
            help='Mask used for ligand')

        parser.add_argument('-p',
            dest='complex_prmtop',
            type=str,
            required=False,
            default='../common/start.prmtop')

        parser.add_argument('-r',
            dest='receptor_prmtop',
            type=str,
            required=False,
            default='rec.prmtop')

        parser.add_argument('-l',
            dest='ligand_prmtop',
            type=str,
            required=False,
            default='ligand.prmtop')

        parser.add_argument('-d',
            dest='dcdfile',
            required=False,
            default='../md.dcd',
            type=str)

        parser.add_argument('-nt',
            dest='ncpus',
            required=True, 
            type=int)

        parser.add_argument('-en',
            dest='is_entropy',
            action='store_true',
            default=False)

        return parser

    def prepare_mmpbsa(self, complex_prmtop, receptor_prmtop, ligand_prmtop, solvent_mask, ligand_mask, complex_dry_prmtop='complex_dry.prmtop'):

        cmdline = """ante-MMPBSA.py -p %(complex_prmtop)s -c %(complex_dry_prmtop)s -r %(receptor_prmtop)s -l %(ligand_prmtop)s \
-s \'%(solvent_mask)s\' -n \'%(ligand_mask)s\'"""%locals()
        print cmdline
        subprocess.check_output(cmdline, shell=True)

    def do_mmpbsa(self, args):

        ncpus = args.ncpus
        complex_prmtop = args.complex_prmtop
        receptor_prmtop = args.receptor_prmtop
        ligand_prmtop = args.ligand_prmtop
        solvent_mask = args.solvent_mask
        ligand_mask = args.ligand_mask
        dcdfile = args.dcdfile

        self.prepare_mmpbsa(complex_prmtop, receptor_prmtop, ligand_prmtop, solvent_mask, ligand_mask)

        if args.is_entropy:
            self.prepare_entropy_input_file(startframe=1, endframe=1000)
            mmpbsa_cmdline = """mpiexec -n %(ncpus)s MMPBSA.py.MPI -O -i mm_entropy.in -o mm_entropy.out -sp %(complex_prmtop)s -cp complex_dry.prmtop \
-rp %(receptor_prmtop)s -lp %(ligand_prmtop)s -y %(dcdfile)s"""%locals()
            print mmpbsa_cmdline
        else:
            self.prepare_pbsa_input_file(startframe=1, endframe=1000)
            mmpbsa_cmdline = """mpiexec -n %(ncpus)s MMPBSA.py.MPI -O -i mm.in -o mm.out -sp %(complex_prmtop)s -cp complex_dry.prmtop \
-rp %(receptor_prmtop)s -lp %(ligand_prmtop)s -y %(dcdfile)s"""%locals()
            print mmpbsa_cmdline

        subprocess.check_output(mmpbsa_cmdline, shell=True)

    def run(self):

        parser = self.create_arg_parser()
        args = parser.parse_args()
        self.do_mmpbsa(args)

if __name__ == '__main__':
    MMPBSAWorker().run()
