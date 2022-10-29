import sys
import setuptools

def exit_with_error(head, body=''):
    _print_admonition('error', head, body)
    sys.exit(1)

# check Python version
if not (sys.version_info[0] == 2 and sys.version_info[1] >= 6):
    exit_with_error("You need Python 2.6.x or Python 2.7.x to install the MDKit package!")

setuptools.setup(name='mdkit',
    version="0.1.4",
    packages=['mdkit', 'mdkit.amber', 'mdkit.namd', 'mdkit.utility'],
    package_data = {'mdkit.amber': ['PROTON_INFO', 'atomic_ions.cmd']},
    scripts = ['bin/prepare_md.py', 'bin/pdb2mol2'],
    install_requires = ['networkx<=2.2', 'numpy<=1.16', 'decorator<=4.4.2', 'cython<=0.29'],
    license='LICENSE.txt',
    description='Tools to prepare and analyze MD simulations conducted with popular MD packages')
