#!/usr/bin/env python3

import os
import sys
import File_operat_mi as FO

def mode_select(ARG):
	'''
	Select the way that the inputs are processed.
	The function should be called as mode_select(sys.argv)
	'''
	if '-g' in ARG:
		return 'gaussian'
	elif '-o' in ARG:
		return 'orca'
	elif '-xyz' in ARG:
		return 'xyz'
	else:
		return 'Please specify the requested input format using either \'-g\' (gaussian), \'-o\' (orca) or \'-xyz\' (xyz format)!'

def main():
	'''
	A program to generate input files for Gaussian (.gjf) or ORCA (.inp) calculations.
	Usage examples:
		a) ./makeinput -i input.gjf -s sample.gjf (or -s sample.inp)
		b) ./makeinput -i input.xyz -s sample.gjf (or -s sample.inp)
		c) ./makeinput -i output.log -s sample.gjf (or -s sample.inp)
		d) ./makeinput -i output.log -o (or -g or -xyz)
		e) ./makeinput -s sample.gjf (or -s sample.inp)
		f) ./makeinput -o (or -g or -xyz)
	- a) covers the case where an arbitrary Gaussian input is
	  transformed into a Gaussian or ORCA input with the desired command line.
	- The "arbitrary input" can also be a .xyz or Gaussian output (.log) file.
	- c-d) takes the last geometry from a Gaussian output and 
	  makes an ORCA (-o), Gaussian (-g) or .xyz (-xyz) "input" file.
	- e) does the same as a) but for all .gjf files in the folder
	- f) does the same as d) but for all .log files in the folder
	- Processing ORCA outputs (.out) is a work in progress.
	'''
	# set up some defaults
	mode = 'orca' # default mode is to write inputs in orca format
	samplefiles = False
	redox="no"
	mark1='CH MU'
	mark2='COORD'
	mark3='SLVT'
	mark4='GENBASIS'
	mark5='XYZFILENAME'
	infile = False
	logfile = False
	basisfile = False
	# parse the arguments
	for n in range(len(sys.argv)):
		if '-g' in sys.argv[n]:
			mode = 'gaussian'
		elif '-o' in sys.argv[n]:
			mode = 'orca'
		elif '-xyz' in sys.argv[n]:
			mode = 'xyz'
		if '-s' in sys.argv[n]:
			if '.gjf' in sys.argv[n+1] or '.inp' in sys.argv[n+1] or '.dat' in sys.argv[n+1]: # a sample file must follow the -s option
				samplefiles = [sys.argv[n+1]]
			else: # if a folder name is given after -s
				samplefiles = os.listdir(sys.argv[n+1])
				# use the absolute paths for the sample files:
				samplefiles = [os.path.join(os.path.abspath(sys.argv[n+1]), i) for i in samplefiles]
		if '-i' in sys.argv[n]: # an input file from where the xyz data is taken must follow the -i option
			if '.gjf' in sys.argv[n+1] or '.xyz' in sys.argv[n+1]:
				infile = sys.argv[n+1]
			elif '.log' in sys.argv[n+1]:
				logfile = sys.argv[n+1]
		if '-b' in sys.argv[n]:
			basisfile = sys.argv[n+1]
		if '-r' in sys.argv[n]:
			redox = "red"
		if '-x' in sys.argv[n]: # o is for ORCA so ox should not be searched for safety
			redox = "ox"
	# do the work
	if samplefiles != False and logfile == False and infile != False:
		for i in samplefiles:
			FO.modify_input(infile, i, mark1, mark2, mark3, mark4, mark5)
	elif samplefiles != False and logfile != False and infile == False:
		for i in samplefiles:
			FO.makeinput1(logfile,True,'gaussian',redox)
			infile = str(logfile)[:-4] + '_coords.gjf'
			FO.modify_input(infile, i, mark1, mark2, mark3, mark4, mark5)
			os.remove(infile)
	elif samplefiles == False and logfile != False and infile == False:
		FO.makeinput1(logfile,True,mode,redox)
	elif samplefiles != False and logfile == False and infile == False: # it will make inputs from .gjf files (!= samplefile) in the folder
		files = []
		for i in os.listdir():
			if '.gjf' in i and '_oxidized' not in i and '_reduced' not in i:
				files.append(i)
			if '.gjf' in i and samplefiles != False:
				files.append(i)
		for i in files:
			for j in samplefiles:
				if i != j and '.gjf' in i and os.path.basename(j)[:-4] in i:
					if basisfile == False: # for OPCs (no genbasis)
						FO.modify_input(i, j, mark1, mark2, mark3, mark4, mark5)
					else: # for TMPCs, need to provide a basis file via the -b option
						FO.modify_input(i, j, mark1, mark2, mark3, mark4, mark5, basisfile)
	elif samplefiles == False and logfile == False and infile == False: # it will make inputs from .log files in the folder
		files = os.listdir()
		for i in files:
			if '.log' in i:
				FO.makeinput1(i,True,mode,redox)

		
if __name__ == "__main__":
    main()


