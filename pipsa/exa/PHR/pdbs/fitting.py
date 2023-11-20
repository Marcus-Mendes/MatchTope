#!/usr/bin/env python
import glob,os,sys,shutil, subprocess, warnings
#import __main__
#__main__.pymol_argv = [ 'pymol', '-qei' ]
#import pymol
#pymol.finish_launching()
#from pymol.cgo import *
#from pymol import cmd

# pymol environment
moddir='/usr/share/pymol'
sys.path.insert(0, moddir)
os.environ['PYMOL_PATH'] = os.path.join(moddir, '')
warnings.filterwarnings("ignore", category=FutureWarning)

# pymol launching
import pymol
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
#Fitting Process:
pwd = subprocess.check_output(['pwd'])
pwd = pwd[:-1]
os.chdir(pwd)
cmd = pymol.cmd
for file in glob.glob("*.pdb"):
	open(file)
	os.rename(file, "target.pdb")
	os.rename("modelo.txt", "modelo.pdb")
	cmd.load ("target.pdb")
	cmd.load ("modelo.pdb")	
	cmd.align('target', 'modelo', cutoff=10.0, cycles=5, gap=-10.0, extend=-0.5, max_gap=50, object=None, matrix='/miniconda/envs/matchtope_env/share/pymol/data/pymol/matrices/BLOSUM62', mobile_state=0, target_state=0, quiet=1, max_skip=0, transform=1, reset=0)
	cmd.save ("target_fit.pdb", "target", -1)
	os.rename("target_fit.pdb", file)
	os.system("rm target.pdb")
	os.rename("modelo.pdb", "modelo.txt")
	cmd.reinitialize()
	
for pdbs in glob.glob("*.pdb"):
	updated_data = ''
	# opening the file
	with open(pdbs, 'r+') as file:
	    # read the file content
	    file_content = file.readlines()

	    # iterate over the content
	    for line in file_content:

	        # removing last word
	        updated_line = ' '.join(line.split(' ')[:-5])

	        # appending data to the variable
	        updated_data += f'{updated_line}\n'

	    # removing the old data
	    file.seek(0)
	    file.truncate()

	    # writing the new data
	    file.write(updated_data)