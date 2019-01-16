import os
import sys
import argparse
import shutil
from glob import glob
from mdtools.utility import mol2
from mdtools.amber import ambertools

parser = argparse.ArgumentParser(description="Prepare rescoring of MD frames")

parser.add_argument('-l',
    type=str,
    dest='files_l',
    nargs='+',
    required=True,
    help='Files ligands')

parser.add_argument('-r',
    type=str,
    dest='files_r',
    nargs='+',
    required=True,
    help='Files receptor')

parser.add_argument('-rref',
    type=str,
    dest='file_r_ref',
    required=True,
    help='File receptor ref')

parser.add_argument('-lref',
    type=str,
    dest='file_l_ref',
    required=True,
    help='File ligand ref')

parser.add_argument('-ligname',
    type=str,
    dest='ligname',
    default='LIG',
    help='Ligand name (default: LIG)')

args = parser.parse_args()

values = ambertools.compute_rmsd(args.files_r, args.files_l, args.file_r_ref, args.file_l_ref, rmsddir='rmsd', ligname=args.ligname, cleanup=False)
shutil.copyfile('rmsd/rmsd.dat', 'rmsd.dat')
shutil.rmtree('rmsd', ignore_errors=True)
