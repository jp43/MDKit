import os

def create_constrained_pdbfile(pdbfile, ligname):

    with open('start.pdb', 'r') as startfile:
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

