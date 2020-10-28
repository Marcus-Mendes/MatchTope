#!/usr/bin/env python
import glob,os,sys,shutil
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

# pymol launching
import pymol
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
#def test():
os.chdir("/home/marcus/Downloads/sw/pipsa/exa/PHR/pdbs/")
cmd = pymol.cmd
for file in glob.glob("*.pdb"):
	#print(file)
	#if file.endswith(".pdb"):
	os.rename(file, "alvo.pdb") #ok
	cmd.load ("alvo.pdb")
	cmd.load ("modelo.txt")	
	cmd.do ("align alvo, modelo")
	#cmd.save ("alvo_fit", "alvo", 0)
	cmd.save ("alvo_fit.pdb", "alvo", -1)
	os.rename("alvo_fit.pdb", file)	
	os.system("rm alvo.pdb")
	#os.system("cd /home/marcus/Dropbox/Doutorado/NetTope/")
	cmd.reinitialize()
	
