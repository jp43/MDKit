import os
from glob import glob
import subprocess
import numpy as np

# required programs
required_amber_programs = ['tleap', 'antechamber', 'parmchk']

def run_shell_command(cmd):
    try:
        subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print e.output
        raise

def check_amber_version():
    if os.environ.get('AMBERHOME'):
        for exe in required_amber_programs:
            try:
                subprocess.check_call('which %s > /dev/null'%exe, shell=True)
            except subprocess.CalledProcessError:
                raise ValueError('Executable %s is not found in your PATH!'%exe)
        docfile = glob(os.environ.get('AMBERHOME')+'/doc/Amber*.pdf')
        amber_version = os.path.basename(docfile[0])[5:-4]
        try:
            int(amber_version)
            return amber_version
        except ValueError:
            raise ValueError("Amber version not detected")
    else:
        raise ValueError("AMBERHOME is not set! Check Amber manual") 

def center_of_geometry(coords):
    cog = np.array([0.0, 0.0, 0.0])
    coords = np.array(coords)
    num_atom = coords.shape[0]
    for idx in xrange(num_atom):
        cog = cog + coords[idx]
    cog = cog / num_atom
    return cog


