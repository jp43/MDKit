import os

def create_constrained_pdbfile(pdbfile, startpdb, ligname):

    with open(startpdb, 'r') as startfile:
        with open(pdbfile, 'w') as posresfile:
            for line in startfile:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    if res_name == 'WAT': # water molecules
                        newline = line[0:30] + '%8.3f'%0.0 + line[38:]
                    elif res_name in ligname: # atoms of the ligand
                        if atom_name.startswith(('C', 'N', 'O')):
                            newline = line[0:30] + '%8.3f'%50.0 + line[38:]
                        else:
                            newline = line[0:30] + '%8.3f'%0.0 + line[38:]
                    else: # atoms of the protein
                        if atom_name in ['C', 'CA', 'N', 'O']:
                            newline = line[0:30] + '%8.3f'%50.0 + line[38:]
                        else:
                            newline = line[0:30] + '%8.3f'%0.0 + line[38:]
                else:
                    newline = line
                print >> posresfile, newline.replace('\n','')

def create_steered_constrained_pdbfile(pdbfile, startpdb, ligname):

    has_found_first_atom = False 
    has_found_last_atom = False 
    # check pdb once to find first and last C alpha atoms
    with open(startpdb, 'r') as startfile:
        for line in startfile:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                atom_num = line[6:11].strip()
                if not has_found_first_atom and atom_name == 'CA':
                    atom_num_first = atom_num
                    has_found_first_atom = True 
                elif atom_name == 'CA':
                    atom_num_last = atom_num
                    has_found_last_atom = True

    if not has_found_first_atom:
        raise ValueError("First CA atom not found in %s"%startpdb)

    if not has_found_last_atom:
        raise ValueError("Last CA atom not found in %s"%startpdb)

    with open(startpdb, 'r') as startfile:
        with open(pdbfile, 'w') as posresfile:
            for line in startfile:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_num = line[6:11].strip()
                    if atom_num == atom_num_first:
                        newline = line[0:56] + "0.00  1.00\n"
                    elif atom_num == atom_num_last:
                        newline = line[0:56] + "1.00  0.00\n"
                    else:
                        newline = line[0:56] + "0.00  0.00\n"
                else:
                    newline = line
                posresfile.write(newline)
