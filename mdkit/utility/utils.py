import os
from glob import glob
import subprocess
import reader
import shlex
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

def extract_amd_lines(parsed_string):
    """Set options for enhanced sampling methods"""
    if parsed_string is not None:
        options = parsed_string.split(',')
        noptions = len(options)/2
 
        if len(options) != 2*noptions:
            raise IOError('Number of options provided for aMD should be even!')
        if 'alphap' in options and 'alphad' in options:
            iamd = 3
        elif 'alphad' in options:
            iamd = 2
        elif 'alphap' in options:
            iamd = 1
        else:
            sys.exit("Option alphad or alphap should be provided with -amd flag")
        lines = "\niamd=%i"%iamd

        for idx in range(noptions):
            key = str(options[2*idx])
            value = str(options[2*idx+1])
            if idx%2 == 0:
                lines += "\n"
            else:
                lines += ", "
            lines += "%s=%s"%(key, value)
        return lines
    else:
        return ""

def extract_tgtmd_lines(parsed_string):
    if parsed_string is None:
        return ""

    args_tgtmd = {}
    tgtmd_options = ['tgtmdfrc', 'tgtrmsmask', 'tgtfitmask', 'tgtrmsd']

    key = ""; value = ""
    for item in parsed_string.split(","):
        if item in tgtmd_options:
            if key and value:
                args_tgtmd[key] = value
                value = ""
            key = item
        elif key and not value:
            value = item
        elif key and value:
            value += ',' + item
    if key and value:
        args_tgtmd[key] = value

    tgtmd_lines = ''
    if args_tgtmd:
        tgtmd_lines = "\nitgtmd=1"
        for key, value in args_tgtmd.iteritems():
            if key in ['tgtrmsmask', 'tgtfitmask']:
                value = '\'' + value + '\''
            tgtmd_lines += ", %s=%s"%(key, value)
    return tgtmd_lines
