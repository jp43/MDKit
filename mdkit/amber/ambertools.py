import os
import sys
import csv
import re
import shutil
import subprocess
import numpy as np

from mdkit.utility import reader
from mdkit.utility import utils

import clustering

def get_nwaters(logfile):

    with open(logfile, 'r') as logf:
        for line in logf:
            line_s = line.split()
            if len(line_s) == 3 and line_s[0] == 'Added' and line_s[-1] == 'residues.':
                return int(line_s[1])

def get_removed_waters(file_r, files_l, file_c, nwaters_tgt, boxsize, step=0.01, ntries=5):

    if not files_l:
        file_c = file_r

    # determine the number of solute residues
    nsolutes = int(subprocess.check_output("echo `cpptraj -p %(file_c)s -mr '*' | awk '{print NF - 1;}'`"%locals(), shell=True))

    print "Targeted number of water residues:", nwaters_tgt
    print "Number of residues found for solute:", nsolutes

    distance = [boxsize]
    for idx in range(ntries):
        distance.append(boxsize + (idx+1)*step)
        distance.append(boxsize - (idx+1)*step)

    diff_best = 1e10
    for d in distance:
        c = 1.0
        lastdiff = 1e10
        while True:
            prepare_leap_config_file('leap.in', file_r, files_l, file_c, solvate=True, distance=d, closeness=c)
            utils.run_shell_command('tleap -f leap.in > leap.log')
            nwaters = get_nwaters('leap.log')
            diff = nwaters - nwaters_tgt
            #print d, c, diff, nwaters
            if diff > 0:
                if lastdiff < 0:
                    if diff < diff_best:
                        diff_best = diff
                        dbest = d
                        cbest = c
                        nwaters_best = nwaters
                    break
                else:
                   lastdiff = diff
                c += 0.01
            elif diff < 0:
                c -= 0.01
                lastdiff = diff
            elif diff == 0:
                diff_best = diff
                dbest = d
                cbest = c
                nwaters_best = nwaters
                break
        if diff == 0:
            break

    print "Closest number of water residues found:", nwaters_best
    print "Removing %i water residues..."%diff_best

    removed_waters = [nsolutes + nwaters_best - diff_best + idx + 1 for idx in range(diff_best)]
    return removed_waters, dbest, cbest

def compute_rmsd(files_r, files_l, file_r_ref, file_l_ref, rmsddir='rmsd', ligname='LIG', cleanup=True):

    files_r_abs = [os.path.abspath(file_r) for file_r in files_r]
    files_l_abs = [os.path.abspath(file_l) for file_l in files_l]

    file_r_ref_abs = os.path.abspath(file_r_ref)
    file_l_ref_abs = os.path.abspath(file_l_ref)

    nfiles_r = len(files_r_abs)
    nfiles_l = len(files_l_abs)

    if nfiles_r != nfiles_l:
        raise ValueError('Not the same number of files for receptors and ligands')

    pwd = os.getcwd()

    shutil.rmtree(rmsddir, ignore_errors=True)
    os.mkdir(rmsddir)
    os.chdir(rmsddir)

    # prepare ref receptor 
    prepare_receptor('protein.pdb', file_r_ref_abs)
    shutil.copyfile(file_l_ref_abs, 'ligand.mol2')
    prepare_ligand('protein.pdb', 'ligand.mol2', 'protein-ligand.pdb')
    
    os.mkdir('PDB')
    files_l_new = ['PDB/ligand-%i.mol2'%(idx+1) for idx in range(nfiles_r)]
    files_rl = ['PDB/protein-ligand-%i.pdb'%(idx+1) for idx in range(nfiles_r)]

    for idx in range(nfiles_r):
        prepare_receptor('PDB/protein-%i.pdb'%(idx+1), files_r_abs[idx])
        shutil.copyfile(files_l_abs[idx], files_l_new[idx])
        prepare_ligand('PDB/protein-%i.pdb'%(idx+1), files_l_new[idx], files_rl[idx])
    
    clustering.prepare_leap_config_file('leap.in', files_r_abs, ['ligand.mol2']+files_l_new, ['protein-ligand.pdb']+files_rl)
    subprocess.check_output('tleap -f leap.in', shell=True)

    lines_trajin = ""
    for idx in range(nfiles_r):
        lines_trajin += "trajin %s\n"%files_rl[idx]
    lines_trajin = lines_trajin[:-1]

    # write cpptraj config file to cluster frames
    with open('cpptraj.in', 'w') as file:
        contents ="""parm protein-ligand.prmtop
trajin protein-ligand.pdb
%(lines_trajin)s
rms first @CA,C,N&!:%(ligname)s
rms nofit :%(ligname)s&!@H= out rmsd.dat"""% locals()
        file.write(contents)

    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True)
    os.chdir(pwd)

    values = np.loadtxt(rmsddir+'/rmsd.dat') 
    values = values[1:,1]
    if cleanup:
        shutil.rmtree(rmsddir)

    return values

def get_solvent_mask(pdbfile, residues='WAT'):

    solvent_residues = map(str.strip, residues.split(','))

    # set residue numbers
    resnum = None
    with open(pdbfile) as pdbf:
        for line in pdbf:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                if resname in solvent_residues:
                    if not resnum:
                        resnum_in = line[22:27].strip()
                    resnum = line[22:27].strip()

    resnum_fin = resnum
    solvent_mask = ':%s-%s'%(resnum_in, resnum_fin)

    return solvent_mask

def get_ions_number(logfile, concentration=0.15, version='14'):

    # Nions = Cions * Nwater * 1/55.5 where 55.5 M is the concentration of pure water
    with open(logfile, 'r') as lf:
        for line in lf:
            line_s = line.strip().split()
            if len(line_s) > 2:
                if line_s[0] == 'Added' and line_s[2] == 'residues.':
                    nresidues = int(line_s[1])
                    ncl = int(round(nresidues * concentration / 55.5))
                    nna = ncl
            if line.startswith("WARNING: The unperturbed charge"):
                net_charge = int(round(float(line_s[7])))
                if net_charge > 0:
                    ncl += abs(net_charge)
                elif net_charge < 0:
                    nna += abs(net_charge)
    return nna, ncl

def load_PROTON_INFO():

    filename = os.path.dirname(os.path.abspath(__file__)) + '/PROTON_INFO'
    info = {}

    with open(filename) as ff:
        ff.next() # skip first line
        for line in ff:
            line_s = line.split()
            is_residue_line = len(line_s) == 2 and line_s[1].isdigit()
            is_hydrogen_line = len(line_s) >= 4 and \
                all([c.isdigit() for c in line_s[:4]])
            is_heavy_atom_line = not is_residue_line and \
                not line_s[0].isdigit()

            if is_residue_line:
                resname = line_s[0]
                info[resname] = []
            elif is_hydrogen_line:
                info[resname].extend(line[15:].split())
            elif is_heavy_atom_line:
                info[resname].extend(line_s)

    no_h_residues = ['PRO']
    for resname in info:
        if resname not in no_h_residues:
            info[resname].append('H')

    info['ACE'] = ['CH3', 'C', 'O', 'HT1', 'HT2', 'HT3']
    info['NME'] = ['CH3', 'C', 'N', 'HT1', 'HT2', 'HT3', 'H']
    return info

def load_atomic_ions():
    """Load formal charge libraries of monoatomic ions"""

    filename = os.path.dirname(os.path.abspath(__file__)) + '/atomic_ions.cmd'
    info = {}
    with open(filename) as ff:
        for line in ff:
            if line.startswith('i = createAtom'):
                charge = float(line.split()[5])
                is_new_atom = True
            elif line.startswith('r = createResidue') and is_new_atom:
                resname = line.split()[3]
                info[resname] = charge
            elif not line.strip():
                is_new_atom = False
    return info

def get_ligand_name(filename):

    # finding ligand name...
    prefix, ext = os.path.splitext(filename)

    r_ligand = reader.open(filename)
    struct = r_ligand.next()
    if ext == '.pdb':
        idxname = 2
    elif ext == '.mol2':
        idxname = 7
    else:
        raise ValueError('File format for file_l should be .pdb or .mol2!')
    names = [row[idxname] for row in struct['ATOM']]
    name = list(set(names))

    if len(name) != 1:
        raise IOError('More than one ligand found in file %s'%filename)
    else:
        name = name[0]
    return name

def correct_hydrogen_names(file_r, keep_hydrogens=False):

    chainIDs = []
    atoms_info = load_PROTON_INFO()

    nremoved = 0
    removed_lines = []

    chainID = None
    is_first_residue = True
    first_residues = []

    # determine which residues are first residues
    with open(file_r, 'r') as rf:
        for line in rf:
            if line.startswith('ATOM'): # atom line
                resnum = line[22:27].strip()

                if is_first_residue:
                    first_residues.append(resnum)
                    is_first_residue = False

            elif line.startswith('TER'):
                is_first_residue = True

    resnum = ''
    with open(file_r, 'r') as rf:
        with open('tmp.pdb', 'w') as wf:
            for line in rf:
                remove_line = False
                if line.startswith('ATOM'): # atom line
                    resname = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    chainID = line[21:22].strip()
                    resnum = line[22:27].strip()

                    if resname in atoms_info:
                       # atom (if atom name starts with a digit, correct it)
                        if atom_name[0].isdigit():
                            atom_name = atom_name[1:] + atom_name[0]

                        # check if hydrogen should be removed
                        if atom_name[0] == 'H':
                            is_hydrogen_from_nterminal = resnum in first_residues and atom_name == 'H'
                            is_hydrogen_known = atom_name in atoms_info[resname] and not is_hydrogen_from_nterminal
                            if keep_hydrogens and not is_hydrogen_known:
                                #print line
                                remove_line = True
                                removed_lines.append(line)
                                #print hydrogens_info[resname], atom_name
                                nremoved += 1
                            elif not keep_hydrogens:
                                remove_line = True
                                nremoved += 1
                        # check if non-hydrogen atom should be removed
                        else:
                            is_atom_known = atom_name in atoms_info[resname]
                            if not is_atom_known:
                                remove_line = True
                                removed_lines.append(line)
                                nremoved += 1

                if not remove_line:
                    wf.write(line)
    #print '\n'.join(removed_lines)
    shutil.move('tmp.pdb', file_r)
    #print "Number of atom lines removed: %s" %nremoved

def run_antechamber(infile, outfile, at='gaff', c='gas', logfile='antechamber.log', version='14', skip_unrecognized_atoms=False):
    """ use H++ idea of running antechamber multiple times with bcc's charge method to estimate the appropriate net charge!!"""
    suffix, ext = os.path.splitext(infile)
    ext = ext[1:]

    supported_charge_methods = ['gas', 'bcc', 'mul']
    max_net_charge = 30
    net_charge = [0]

    unrecognized_atom = False
    if c and c.lower() != 'none':
        c_lower = c.lower()
        if c_lower not in supported_charge_methods:
            raise ValueError("Charge method %s not recognized! Charge method should be one of %s"%(c_lower,', '.join(supported_charge_methods)))

        for nc in range(max_net_charge):
            net_charge.extend([nc+1,-(nc+1)])

        for nc in net_charge:
            is_error = False
            is_wrong_charge = False
            if version in ['14', '15']:
                command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -c %(c)s -nc %(nc)s -du y -pf y &>> %(logfile)s'%locals()
            elif version in ['16', '17']:
                command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -c %(c)s -nc %(nc)s -du y -pf y -dr no &>> %(logfile)s'%locals()  
            utils.run_shell_command('echo "# command used: %(command)s" > %(logfile)s'%locals()) # print command in logfile
            try:
                utils.run_shell_command(command)
            except subprocess.CalledProcessError:
                with open(logfile, 'r') as lf:
                    for line in lf:
                        if line.startswith("No Gasteiger parameter for atom"):
                            unrecognized_atom = True
                if unrecognized_atom and skip_unrecognized_atoms:
                    break 
            with open(logfile, 'r') as lf:
                for line in lf:
                    line_s = line.split()
                    if c_lower == 'gas':
                        if version in ['14', '15'] and 'does not equal to the total formal charge' in line:
                            nc_suggested = int(float(line_s[8][1:5]))
                            is_wrong_charge = True
                        elif version in ['16', '17'] and 'Warning: The net charge of the molecule' in line and 'does not equal the' in line:
                            nc_suggested = int(float(line_s[16][1:5]))
                            is_wrong_charge = True
                        if is_wrong_charge:
                            if nc_suggested == nc:
                                return
                            break
                    elif c_lower in ['bcc', 'mul']:
                        if version in ['14', '15'] and 'Error: cannot run' in line and 'sqm' in line:
                            is_wrong_charge = True
                        elif version in ['16', '17'] and 'Cannot properly run' in line and 'sqm' in line:
                            is_wrong_charge = True
                    if 'Error' in line:
                        is_error = True
                if not is_error and not is_wrong_charge:
                    break

        if is_wrong_charge:
            raise ValueError("No appropriate net charge was found to run antechamber's %s charge method"%c)
        elif is_error:
            raise ValueError("Error when running antechamber!")

    if c is None or c.lower() == 'none' or (unrecognized_atom and skip_unrecognized_atoms): # do not regenerate charges
       if version in ['14', '15']:
           command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -du y -pf y > %(logfile)s'%locals()
       elif version in ['16', '17']:
           command = 'antechamber -i %(infile)s -fi %(ext)s -o %(outfile)s -fo mol2 -at %(at)s -du y -pf y -dr no > %(logfile)s'%locals()
       utils.run_shell_command(command)

def prepare_leap_config_file(script_name, file_r, files_l, file_rl, solvate=False, PBRadii=None, forcefield='ff14SB', nna=0, ncl=0, box='parallelepiped', distance=10.0, closeness=1.0, model='TIP3P', version='14', membrane=False):
 
    solvation_line = ""
    pbradii_lines = ""
    add_ions_lines = ""
    ligand_lines = ""
    ions_libraries_lines = ""
    loadpdb_line = ""
    set_box_line = "" # used when membrane is True

    tip3p_models = ['TIP3P', 'TIP3PF', 'POL3', 'QSPCFW']
    tip4p_models = ['TIP4P', 'TIP4PEW']
    spc_models = ['QSPCFW', 'SPC', 'SPCFW']

    if model.upper() in tip3p_models:
        solvent_model = 'tip3p'
    elif model.upper() in tip4p_models:
        solvent_model = 'tip4pew'
    elif model.upper() in spc_models:
        solvent_model = 'spce'
    else:
        raise ValueError('Solvent model %s unknown, should be one of !'% \
(model,', '.joined(tip3p_models + tip4p_models + spc_models)))


    if version in ['14', '15']:
        forcefield_lines = 'source leaprc.' + forcefield 
    elif version in ['16', '17']:
        forcefield_lines = 'source leaprc.protein.' + forcefield
        forcefield_lines += '\n' + 'source leaprc.water.' + solvent_model
    else:
        raise ValueError("Amber version %s not supported!"%version)

    if solvate:
        boxtype = model.upper() + 'BOX'
        if box.lower() == 'parallelepiped':
            solvation_line = "\nSolvateBox complex %s %.2f %.2f"%(boxtype,distance,closeness)
        elif box.lower() == 'octahedron':
            solvation_line = "\nSolvateOct complex %s %.2f %.2f"%(boxtype,distance,closeness)
        else:
            raise ValueError('Box type %s not recognized'%box) 

        if nna > 0:
            add_ions_lines += "\naddions complex Na+ %i"%nna
        if ncl > 0:
            add_ions_lines += "\naddions complex Cl- %i"%ncl
        suffix_ions_libraries = solvent_model
        if nna > 0 or ncl > 0: 
            ions_libraries_lines = """\nloadamberparams frcmod.ionsjc_%(suffix_ions_libraries)s
loadamberparams frcmod.ionslm_1264_%(suffix_ions_libraries)s"""%locals()
    elif membrane:
        forcefield_lines += '\nsource leaprc.lipid14'
        box_dx, box_dy, box_dz = utils.get_box_dimensions(file_r, mask=['WAT'])
        set_box_line = "\nset complex box { %.3f %.3f %.3f }"%(box_dx, box_dy, box_dz)

    if PBRadii:
        pbradii_lines = "\nset default PBRadii %s"%PBRadii

    if files_l:
        forcefield_lines += '\nsource leaprc.gaff'
        if isinstance(files_l, basestring):
            files_l = [files_l]
        for idx, file_l in enumerate(files_l):
            file_l_prefix, ext = os.path.splitext(file_l)
            file_l_prefix = os.path.basename(file_l_prefix)
            name = get_ligand_name(file_l)
            ligand_lines += "\n%(name)s = loadmol2 %(file_l)s\nloadamberparams %(file_l_prefix)s.frcmod"%locals()
        loadpdb_line = "complex = loadPdb %s"%file_rl
    else:
        loadpdb_line = "complex = loadPdb %s"%file_r

    with open(script_name, 'w') as leapf:
        script ="""%(forcefield_lines)s%(ions_libraries_lines)s%(ligand_lines)s%(pbradii_lines)s
%(loadpdb_line)s%(solvation_line)s%(add_ions_lines)s%(set_box_line)s
saveAmberParm complex start.prmtop start.inpcrd
savePdb complex start.pdb
quit\n"""%locals()
        leapf.write(script)

def prepare_receptor(file_r_out, file_r, keep_hydrogens=False, membrane=False):

    # only keep atom lines
    with open(file_r, 'r') as tmpf:
        with open(file_r_out, 'w') as recf:
            for line in tmpf:
                # check if non-hydrogen atom line
                if line.startswith(('ATOM', 'HETATM', 'TER')):
                    recf.write(line)
            # if last line not TER, write it
            if not line.startswith('TER'):
                recf.write('TER\n')

    if membrane:
        charmmlipid2amber(file_r_out)

    # remove atoms and hydrogen with no name recognized by AMBER
    correct_hydrogen_names(file_r_out, keep_hydrogens=keep_hydrogens)

def prepare_ligand(file_r, files_l, file_rl, charge_method='gas', version='14', skip_unrecognized_atoms=False):

    if isinstance(files_l, basestring):
        files_l = [files_l]

    mol2files_l = []

    shutil.copyfile(file_r, file_rl)
    for file_l in files_l:
        file_l_prefix, ext = os.path.splitext(file_l)
        file_l_prefix = os.path.basename(file_l_prefix)

        mol2file = file_l_prefix + '.mol2'
        run_antechamber(file_l, 'tmp.mol2', at='gaff', c=charge_method, version=version, skip_unrecognized_atoms=skip_unrecognized_atoms)

        shutil.move('tmp.mol2', mol2file)
        utils.run_shell_command('parmchk -i %s -f mol2 -o %s.frcmod'%(mol2file, file_l_prefix))
        if version in ['14', '15']:
            utils.run_shell_command('antechamber -i %s -fi mol2 -o %s.pdb -fo pdb &> /dev/null'%(mol2file, file_l_prefix))
        elif version in ['16', '17']:
            utils.run_shell_command('antechamber -i %s -fi mol2 -o %s.pdb -fo pdb -dr no &> /dev/null'%(mol2file, file_l_prefix))

        mol2files_l.append(mol2file)
        with open(file_rl, 'a') as ffrl:
            with open(file_l_prefix + '.pdb', 'r') as ffl:
                for line in ffl:
                    if line.startswith(('ATOM', 'HETATM')):
                        ffrl.write(line)
            ffrl.write('TER\n')

    return mol2files_l

def charmmlipid2amber(filename):

    amberhome = os.getenv('AMBERHOME')
    if amberhome == None:
        print "Error: AMBERHOME not set. Set the AMBERHOME environment variable."
        sys.exit(2)
    else:
        # Hard coded into here. This has to be changed if the data file is moved.
        convert_filename = "%s/dat/charmmlipid2amber/charmmlipid2amber.csv" % (amberhome)

    input_file  = open(filename, 'r') # File to be processed
    input_file_list = [] # File to be processed as a list
    for line in input_file:
        input_file_list.append(line)
    input_file.close()

    # Process residues
    residue_list = [] # List of residue numbers
    residue_start = [] # List of lines where residues start. Line numbers start at 1.
    residue_end = [] # List of lines where residues end, including TER card if present. Line numbers start at 1.
    it0 = 1
    previous_residue = "" # Residue number of the previous line (including TER card). Set to "" if it is not an ATOM record or TER card.
    current_residue = "" # Residue number of the current line (including TER card). Set to "" if it is not an ATOM record or TER card.
    # Split file into residues by checking columns 23-27 (<99999 residues):
    # First line:
    if (input_file_list[0][0:6] == "ATOM  " or
        input_file_list[0][0:6] == "HETATM"):
        residue_list.append(input_file_list[0][22:27])
        residue_start.append(it0)
        previous_residue = input_file_list[0][22:27] 
    elif line[0:3] == "TER":
        previous_residue = ""
    else:
        previous_residue = ""
    it0+=1
    # Rest of lines:
    for line in input_file_list[1:]:
        if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
            current_residue = line[22:27]
        elif line[0:3] == "TER":
            current_residue = previous_residue
        else:
            current_residue = ""
        if previous_residue != current_residue:
            # Previous line was not an ATOM or TER:
            if previous_residue == "":
                residue_list.append(current_residue)
                residue_start.append(it0)
                previous_residue = current_residue
            # Current line is not an ATOM or TER:
            elif current_residue == "":
                residue_end.append(it0-1)
                previous_residue = current_residue
            # Previous and current line are ATOM or TER:
            else:
                residue_list.append(current_residue)
                residue_start.append(it0)
                residue_end.append(it0-1)
                previous_residue = current_residue
        it0+=1
    # If the last residue is not closed, define the end:
    if current_residue != "":
        residue_end.append(it0-1)

    # Process substition dictionaries
    try:
        csv_file = open(convert_filename, 'r') # csv file with all substitutions
    except IOError, err:
        print "Error: ", str(err)
        sys.exit(1)
    # Skip header line of csv file. Line 2 contains dictionary keys:
    csv_file.readline()
    csv_file_reader = csv.DictReader(csv_file) # Dictionary csv reader
    replace_dict = {} # Dictionary of atom name and residue name search and replace
    order_dict = {} # Dictionary of atom name and residue name order
    ter_dict = {} # Dictionary of whether residue should have a TER card based on atom name and residue name. All atom name and residue name in a residue with a TER card will return True.
    num_atom_dict = {} # Dictionary of number of atoms in current residue for the search string
    for line in csv_file_reader:
        replace_dict[line["search"]] = line["replace"]
        order_dict[line["search"]] = int(line["order"])
        ter_dict[line["search"]] = (line["TER"] == "True")
        num_atom_dict[line["search"]] = int(line["num_atom"])
    csv_file.close()

    # Do substitions
    # The search and replace is columns 13-21 in the PDB file:
    # 13-16: atom name
    # 17:      alternate location indicator
    # 18-20: residue name
    # 21:      sometimes used for the residue name
    output_file_list = [] # File to be written in list form (after processing)
    residue_substituted = False # For error checking. True if a substitution occurs.
    for it1 in range(0, len(residue_list)):
        # residue_start and residue_end indices start at 1. 
        # input_file_list indices start at 0.
        input_residue = input_file_list[residue_start[it1]-1:
            residue_end[it1]]
        output_residue = []
        # Process residue only if first atom is in the replacement dictionary:
        if input_residue[0][12:21] in replace_dict:
            residue_substituted = True
            # Check if length of input_residue is correct:
            # Count TER cards in residue for residue length arithmetic
            n_TER_cards = 0
            for line in input_residue:
                if line[0:3] == "TER":
                    n_TER_cards += 1
            if len(input_residue)-n_TER_cards != num_atom_dict[input_residue[0][12:21]]:
                print "Error: Number of atoms in residue does not match number of atoms in residue in replacement data file"
                sys.exit(1)
            output_residue = len(input_residue)*[0]
            for it2 in range(0, len(input_residue)):
                line = input_residue[it2]
                if line[0:3] != "TER":
                     search = line[12:21]
                     output_residue[order_dict[search]] = re.sub(search, 
                        replace_dict[search], line)
                else:
                     output_residue[it2] = "TER   \n"
            if ter_dict[input_residue[0][12:21]] == True:
                if input_residue[-1][0:3] != "TER":
                    output_residue.append("TER   \n")
        else:
            output_residue = input_residue
        output_file_list.extend(output_residue)
    # Check if any residues were substituted:
    if residue_substituted == False:
        print "Error: No residues substituted"
        sys.exit(1)

    resname_last = ''
    resnum_last = 0
    resum_new = 0
    new_output_file_list = []
    for line in output_file_list:
        if line.startswith(('ATOM', 'HETATM')):
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            resnum = int(line[22:27])
            if resname != resname_last or resnum != resnum_last:
                resum_new = (resum_new+1)%100000
            # (1) Add TER lines if needed
            if (resnum_last > resnum) or (resname == 'WAT' and atomname == 'O') or atomname in ['Na+', 'Cl-', 'K+']:
                new_output_file_list.append('TER\n')
            # (2) fix atom names
            if resname == 'ILE' and atomname == 'CD':
                atomname_new = ' CD1'
            else:
                atomname_new = line[12:16]
            # (3) fix residue names
            if resname == 'HSD':
                resname_new = 'HIS'
            else:
                resname_new = line[17:20]
            resname_last = resname
            resnum_last = resnum
            resum_new_str = str(resum_new)
            if len(resum_new_str) <= 4:
                resum_new_str = (4-len(str(resum_new)))*' ' + resum_new_str + ' '
            elif len(resum_new_str) == 5:
                resum_new_str = (5-len(str(resum_new)))*' ' + resum_new_str
            else:
                raise Exception("Something went wrong!")
            line_new = line[:12] + atomname_new + line[16:17] + resname_new + line[20:22] + resum_new_str + line[27:]
        else:
            line_new = line
        new_output_file_list.append(line_new)

    last_line = ''
    output_file_list = []
    for line in new_output_file_list:
        if not last_line.startswith('TER') or not line.startswith('TER'):
            output_file_list.append(line)
        last_line = line

    # Write output
    try:
        output_file = open(filename, 'w') # File to be written
    except IOError, err:
        print "Error: ", str(err)
        sys.exit(1)
    for line in output_file_list:
        output_file.write(line)
    output_file.close()
