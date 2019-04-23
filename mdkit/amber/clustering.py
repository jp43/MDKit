import os
import stat
import shutil
import subprocess
import argparse
import ambertools

def cluster_trajectory(crdfile, prmtop, mode='clustering', cutoff=None, nclusters=None, cleanup=True, trajin='', ligname='LIG'):
    # get current directory
    curdir = os.getcwd()

    # create directory where clustering (or fit, 2rmsd...) will be performed
    workdir = mode
    shutil.rmtree(workdir, ignore_errors=True) # overwrite by default
    os.mkdir(workdir)

    # get abolute paths of trajectory and prmtop files
    if os.path.isfile(crdfile):
        crdfile_abs = os.path.abspath(crdfile)
    else:
        raise ValueError("No crd file found: %s"%crdfile)

    if os.path.isfile(prmtop):
        prmtop_abs = os.path.abspath(prmtop)
    else:
        raise ValueError("No prmtop file found: %s"%prmtop)

    # change working directory
    os.chdir(workdir)

    if mode == 'clustering':
        print "Clustering trajectory..."
    elif mode == 'fit':
        print "RMSD fit trajectory..."
    elif mode == 'rmsd2d':
        print "Computing RMSD matrix from trajectory..."

    # run clustering (or fit, 2rmsd...) using cpptraj
    prepare_cpptraj_config_file_trajectory('cpptraj.in', crdfile_abs, prmtop_abs, cutoff=cutoff, nclusters=nclusters, mode=mode, trajin=trajin, ligname=ligname)
    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True)

    os.chdir(curdir)
    print "done."

def prepare_rmsd2d(crdfile, prmtop, trajin=''):

    mask_solvent = ':WAT,Na+,Cl-'

    # read number of frames in trajectory using cpptraj
    with open('cpptraj.in', 'w') as inf:
        contents = """parm %(prmtop)s
trajin %(crdfile)s %(trajin)s
strip %(mask_solvent)s outprefix noW
trajout traj.mdcrd\n"""% locals()
        inf.write(contents)

    # use debug option to guarantee that the number of frames will show up
    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True)
    with open('cpptraj.log') as logf:
        for line in logf:
            if 'frames and processed' in line:
                # get total number of frames in trajectory
                nframes = int(line.split()[-2])

    # from trajin options provided, deduce indices of frames that will actually be loaded
    if trajin:
        trajin_s = map(int, trajin.split())
        first = trajin_s[0]
        if len(trajin_s) > 1:
            last = min(trajin_s[1], nframes)
            if len(trajin_s) > 2:
                off = trajin_s[2]
            else:
                off = 1
        else:
            last = nframes
            off = 1
    else:
        first = 1
        last = nframes
        off = 1

    nframes = (last-first)/off
    indices = [first + idx*off for idx in range(nframes+1)]
    return indices


def cluster_poses(files_r, files_l, mode='clustering', cutoff=None, nclusters=None, cleanup=True, ligname='LIG'):

    # get current directory
    curdir = os.getcwd()

    # create directory where minimization will be performed
    workdir = mode
    shutil.rmtree(workdir, ignore_errors=True)
    os.mkdir(workdir)

    # get full path of ligand files
    if len(files_l) >= 2:
        files_l = [os.path.abspath(file_l) for file_l in files_l]
    else:
        raise ValueError('At least 2 ligand files are required for clustering')

    # get full path of receptor files
    if len(files_r) == 1:
        nfiles_l = len(files_l)
        files_r = [os.path.abspath(files_r[0]) for idx in range(nfiles_l)]
    else:
        files_r = [os.path.abspath(file_r) for file_r in files_r]

    # Check if same number of receptors and ligands 
    if len(files_r) != len(files_l):
        raise ValueError('Number of receptors and ligands should be the same!')

    # change working directory
    os.chdir(workdir)

    new_files_r = []
    # prepare receptors
    for idx, file_r in enumerate(files_r):
        new_file_r = 'protein-%s.pdb'%idx
        # prepare receptor
        ambertools.prepare_receptor(new_file_r, file_r)
        new_files_r.append(new_file_r)

    if mode == 'clustering':
        print "Clustering poses..."
    elif mode == 'fit':
        print "RMSD fit poses..."
    elif mode == '2rmsd':
        print "Computing RMSD matrix of poses..."

    # run clustering (or fit, 2rmsd...) using cpptraj
    do_amber_clustering(new_files_r, files_l, mode, cutoff=cutoff, nclusters=nclusters, cleanup=cleanup, ligname=ligname)
    os.chdir(curdir)

def prepare_leap_config_file(filename, files_r, files_l, files_rl, forcefield='leaprc.ff14SB'):

    linespdb = ""
    for idx, file_rl in enumerate(files_rl):
        if idx == 0:
            linespdb += """p = loadPdb %s
saveAmberParm p protein-ligand.prmtop protein-ligand.inpcrd
savepdb p %s\n"""%(file_rl,file_rl)
        else:
            linespdb += """p = loadPdb %s
savepdb p %s\n"""%(file_rl,file_rl)

    name = ambertools.get_ligand_name(files_l[0])

    linespdb = linespdb[:-1]
    with open(filename, 'w') as ff:
        contents ="""source %(forcefield)s
source leaprc.gaff
%(name)s = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
%(linespdb)s
quit"""% locals()
        ff.write(contents)

def prepare_cpptraj_config_file_trajectory(filename, crdfile, prmtop, cutoff=None, nclusters=None, mode='clustering', trajin='', ligname='LIG'):
    mask_solvent = ':WAT,Na+,Cl-'

    mask = ':%s&!@H='%ligname
    maskfit = '@CA,C,N&!:%s'%ligname

    # write cpptraj config file to cluster frames
    with open(filename, 'w') as file:
        if mode == 'clustering':
            if cutoff and nclusters:
                ValueError('Both cutoff value and nclusters provided. Only one of those parameters should be given!')
            elif cutoff:
                option = " epsilon %s "%cutoff
            elif nclusters:
                option = " clusters %s"%nclusters
            else:
                option = ""
            contents = """parm %(prmtop)s
trajin %(crdfile)s %(trajin)s
strip %(mask_solvent)s
rms first %(maskfit)s
cluster %(mask)s nofit summary summary.dat info info.dat repout frame repfmt pdb%(option)s\n"""% locals()
            file.write(contents)
        elif mode == 'fit':
            contents = """parm %(prmtop)s
trajin %(crdfile)s %(trajin)s
strip %(mask_solvent)s
rms first %(maskfit)s
strip !:%(ligname)s
trajout struct.pdb multi\n"""% locals()
            file.write(contents)
        elif mode == 'rmsd2d':
            # Needs to know the frames indices to compute 2D RMSD matrix
            frames_idxs = prepare_rmsd2d(crdfile, prmtop, trajin=trajin)
            new_prmtop = 'noW.' + os.path.basename(prmtop)
            lines_rms = ""
            for idx, frame_idx in enumerate(frames_idxs):
                jdx = idx + 1
                lines_rms += """reference traj.mdcrd %(frame_idx)s [ref%(jdx)s]
rms ref [ref%(jdx)s] %(maskfit)s 
rms ref [ref%(jdx)s] %(mask)s nofit out rmsd_%(jdx)s.txt\n""" %locals()
            contents = """parm %(new_prmtop)s
trajin traj.mdcrd %(trajin)s
%(lines_rms)s"""% locals()

            file.write(contents)

def prepare_cpptraj_config_file_poses(filename, files_rl, prmtop, cutoff=None, nclusters=None, mode='clustering', ligname='LIG'):

    mask = ':%s&!@H='%ligname
    maskfit = '@CA,C,N&!:%s'%ligname

    lines_trajin = ""
    for file_rl in files_rl:
        lines_trajin += "trajin %s\n"%(file_rl)

    # remove last \n
    lines_trajin = lines_trajin[:-1]

    # write cpptraj config file to cluster frames
    with open(filename, 'w') as file:
        if mode == 'clustering':
            if cutoff and nclusters:
                ValueError('Both cutoff value and nclusters provided. Only one of those parameters should be given!')
            elif cutoff:
                option = " epsilon %s "%cutoff
            elif nclusters:
                option = " clusters %s"%nclusters
            else:
                option = ""

            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
cluster %(mask)s nofit%(option)s summary summary.dat info info.dat\n"""% locals()
            file.write(contents)
        elif mode == 'pca':
            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
createcrd md-trajectories
run
crdaction md-trajectories matrix covar name covar %(mask)s
runanalysis diagmatrix covar out evecs.dat vecs 2 name myEvecs
crdaction md-trajectories projection md-pca modes myEvecs %(mask)s out pca.out\n"""% locals()
            file.write(contents)
        elif mode == 'fit':
            contents = """parm %(prmtop)s
%(lines_trajin)s
rms first %(maskfit)s
trajout ref.rst restart onlyframes 1
trajout struct.pdb multi\n"""% locals()
            file.write(contents)
        elif mode == 'rmsd2d':
            lines_rms = ""
            for idx, file_rl in enumerate(files_rl):
                jdx = idx + 1
                lines_rms += """reference %(file_rl)s [ref%(jdx)s]
rms ref [ref%(jdx)s] %(maskfit)s 
rms ref [ref%(jdx)s] %(mask)s nofit out rmsd_%(jdx)s.txt\n""" %locals()
            contents = """parm protein-ligand.prmtop
%(lines_trajin)s
%(lines_rms)s"""% locals()
            file.write(contents)

def do_amber_clustering(files_r, files_l, mode, cutoff=None, nclusters=None, cleanup=False, ligname='LIG'):

    # (A) Prepare ligand and PDB files
    os.mkdir('PDB')
    files_rl = []
    for idx, file_l in enumerate(files_l):
        file_r = files_r[idx]
        file_rl = 'PDB/protein-ligand-%s.pdb'%(idx+1)
        shutil.copyfile(file_l, 'ligand.mol2')
        ambertools.prepare_ligand(file_r, 'ligand.mol2', file_rl)
        files_rl.append(file_rl)
        os.remove(file_r)
        if idx == 0:
            shutil.copyfile('ligand.mol2', 'ligand_ref.mol2')

    # (B) Run tleap
    prepare_leap_config_file('leap.in', files_r, files_l, files_rl)
    subprocess.check_output('tleap -f leap.in > leap.log', shell=True)

    # (C) Run cpptraj
    prepare_cpptraj_config_file_poses('cpptraj.in', files_rl, 'protein-ligand.prmtop', cutoff=cutoff, nclusters=nclusters, mode=mode, ligname=ligname)
    subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True)

    if cleanup:
        # (D) remove PDB folder and other large files
        shutil.rmtree('PDB', ignore_errors=True)
        os.remove('leap.log')
        os.remove('protein-ligand.prmtop')

def create_arg_parser():

    parser = argparse.ArgumentParser(description="Run Amber Clustering")

    parser.add_argument('-l',
        type=str,
        dest='files_l',
        nargs='+',
        default=None,
        help = 'Ligand coordinate file(s): .mol2')

    parser.add_argument('-r',
        type=str,
        dest='files_r',
        nargs='+',
        default=None,
        help = 'Receptor coordinate file(s): .pdb')

    parser.add_argument('-cleanup',
        dest='cleanup',
        action='store_true',
        default=False,
        help = 'Remove intermediate files')

    parser.add_argument('-ligname',
        dest='ligname',
        default='LIG',
        help = 'Name of the ligand to build the masks')

    parser.add_argument('-mode',
        type=str,
        dest='mode',
        default='clustering',
        help = 'Cpptraj mode (clustering, pca, fit, rmsd2d)')

    parser.add_argument('-n',
        type=str,
        dest='nclusters',
        default=None,
        help = 'Number of clusters for clustering analysis')

    parser.add_argument('-rmsd',
        type=str,
        dest='cutoff',
        default=None,
        help = 'RMSD cutoff for clustering analysis')

    parser.add_argument('-traj',
        action='store_true',
        dest='trajectory',
        default=False,
        help = 'Cluster trajectory (requires -y and -p options)')

    parser.add_argument('-trajin',
        type=str,
        dest='trajin',
        default='',
        help = 'First frame, last frame and offset for trajin when loading trajectory')

    parser.add_argument('-p',
        type=str,
        dest='topfile',
        default=None,
        help = 'Topology file: .prmtop (used with -traj flag)')

    parser.add_argument('-y',
        type=str,
        dest='crdfile',
        default=None,
        help = 'Trajectory file: .mdcrd (used with -traj flag)')

    return parser

def check_args(args):

    if args.trajectory:
        required_args = {'.prmtop file': (args.topfile, '-p'), '.mdcrd file': (args.crdfile, '-y')}
    else:
        required_args = {'ligand file': (args.files_l, '-l'), 'protein file': (args.files_l, '-r')}

    for arg, value in required_args.iteritems():
        if not value[0]:
            raise IOError("No %s provided (%s option)"%(arg, value[1]))

def run():

    parser = create_arg_parser()
    args = parser.parse_args()
    check_args(args)

    if args.trajectory:
        cluster_trajectory(args.crdfile, args.topfile, mode=args.mode, cutoff=args.cutoff, nclusters=args.nclusters, cleanup=args.cleanup, \
ligname=args.ligname, trajin=args.trajin)
    else:
        cluster_poses(args.files_r, args.files_l, mode=args.mode, cutoff=args.cutoff, nclusters=args.nclusters, cleanup=args.cleanup, ligname=args.ligname)
