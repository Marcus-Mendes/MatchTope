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
	os.rename(file, "target.dump")
	cmd.load ("target.dump")
	cmd.load ("modelo.txt")	
	cmd.align('target', 'modelo', cutoff=10.0, cycles=5, gap=-10.0, extend=-0.5, max_gap=50, object=None, matrix='BLOSUM62', mobile_state=0, target_state=0, quiet=1, max_skip=0, transform=1, reset=0)
	cmd.save ("target_fit", "target", -1)
	os.rename("target_fit", file)
	os.system("rm target.dump")
	cmd.reinitialize()
	
