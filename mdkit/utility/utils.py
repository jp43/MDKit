import os
from glob import glob
import subprocess
import reader
import numpy as np

required_amber_programs=['tleap', 'antechamber', 'parmchk']

def run_shell_command(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
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
        raise ValueError("AMBERHOME is not set! Check Amber manual!") 

def center_of_geometry(coords):
    cog = np.array([0.0, 0.0, 0.0])
    coords = np.array(coords)
    num_atom = coords.shape[0]
    for idx in xrange(num_atom):
        cog = cog + coords[idx]
    cog = cog / num_atom
    return cog

def get_box_dimensions(pdbfile, mask=[]):
    # used when prepare pdbfile from CHARMM Membrane Builder
    pdbf = reader.open(pdbfile)
    coords = []
    for line in pdbf.next()['ATOM']:
        if (mask and line[2].strip() in mask) or not mask:
            coords.append(map(float,line[4:]))
    coords = np.array(coords)
    dist = np.amax(coords, axis=0) - np.amin(coords, axis=0)
    return dist
