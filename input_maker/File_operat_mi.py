import os
import sys
#import Energia_4 as energ

ATOMLISTA = {11: 'Na', 35: 'Br', 47: 'Ag', 77: 'Ir', 79: 'Au', 29: 'Cu', 58: 'Ce', 26: 'Fe', 28: 'Ni', 44: 'Ru', 45: 'Rh', 17: 'Cl', 16: 'S', 14: 'Si', 15: 'P', 6: 'C', 8: 'O', 1: 'H', 7: 'N', 5: 'B', 9: 'F'}


def makeinput1(logfile,printtofile=True,mode='gaussian',redox="no"):
	'''
	Extract the final coordinates from a Gaussian output.
	Use the printtofile argument to set behavior:
		-- printtofile == True --> returns a basic input file based on the mode argument
		-- printtofile == Flase --> returns three strings: coordinates, charge, multiplicity
	Options for mode:
		-- 'gaussian'
		-- 'orca'
		-- 'xyz'
	'''
	parse = False
	ch = 0
	mu = 1
	ch_mu_done = False
	iname = str(logfile)[:-4] + '_coords'
	with open(logfile) as inp:
		xyz_lines = '' # initialize xyz_lines for safety reasons
		N = 0 # xyz line counter (= number of atoms counter)
		for line in inp:
			if "Charge" in line and "Multiplicity" in line and ch_mu_done == False:
				q = line.split()
				ch = int(q[2])
				mu = int(q[5])
				if mu%2 == 0: # even multiplicity (e.g., radicals)
					if redox == "ox" or redox == "red":
						mu = mu - 1
				else: # even multiplicity (e.g., singlet ground state)
					if redox == "ox" or redox == "red":
						mu = mu + 1
				if redox == "ox":
					ch = ch + 1
					iname += '_oxidized'
				elif redox == "red":
					ch = ch - 1
					iname += '_reduced'
				ch_mu_done = True
			elif "Standard orientation:" in line or "Input orientation:" in line: # there is no "Standard orientation:" in "nosymm" runs
				parse = True
				xyz_lines = '' # re-initialize xyz_lines so that it will always contain the last set of coordinates
				N = 0 # re-initialize the line counter as well
			elif "Rotational constants (GHZ):" in line or 'Distance matrix (angstroms):' in line:
				parse = False
			if parse == True:
				line.strip()
				dat = line.split()
				if dat != [] and dat[0].isdigit():
					N += 1
					atomnumber = dat[1]
					atomtype = ATOMLISTA[int(atomnumber)]
					x = dat[3]
					y = dat[4]
					z = dat[5]
					xyz_lines += '%-2s' % str(atomtype)+str('%28.8f' % float(x))+str('%14.8f' % float(y))+str('%14.8f' % float(z))+"\n"
	if xyz_lines == '':
		print('Error: no coordinates were found!')
		sys.exit(3)
	if printtofile == True:
		if mode == 'gaussian':
			iname += '.gjf'
			ifile = open(iname, "w")
			ifile.write('%chk=cp.chk\n')
			ifile.write('#p hf/3-21g\n')
			ifile.write('\n')
			ifile.write('Title\n')
			ifile.write('\n')
			ifile.write(str(ch) + ' ' + str(mu) + '\n')
			ifile.write(xyz_lines)
			ifile.write('\n')
			ifile.close()
			print(iname + ' has been writen!')
			#sys.exit(1)
		elif mode == 'orca':
			iname += '.inp'
			ifile = open(iname, "w")
			ifile.write('! B3LYP def2-SVP\n')
			ifile.write('\n')
			ifile.write('* xyz ' + str(ch) + ' ' + str(mu) + '\n')
			ifile.write(xyz_lines)
			ifile.write('*\n')
			ifile.write('\n')
			ifile.close()
			print(iname + ' has been writen!')
			#sys.exit(1)
		elif mode == 'xyz':
			iname += '.xyz'
			ifile = open(iname, "w")
			ifile.write(str(N) + '\n')
			ifile.write('\n')
			ifile.write(xyz_lines)
			ifile.write('\n')
			ifile.close()
			print(iname + ' has been writen!')
			#sys.exit(1)
		else:
			print('The selected file write mode is not valid. Please use mode=\'gaussian\' (-g), \'orca\' (-o) or \'xyz\' (-xyz)')
			sys.exit(2)
	else:
		return xyz_lines, ch, mu

def makeinput2(xyzfile):
	'''
	Extract coordinates from an input file.
	The input file is given through the 'xyzfile' argument and is tested to work with the following:
		-- gaussian input (.gjf)
		-- xyz file (.xyz)
	Returns three strings (coordinates, charge, multiplicity) to be used in other functions.
	'''
	parse = False
	ch = 0
	mu = 1
	with open(xyzfile) as inp:
		xyz_lines = '' # initialize xyz_lines for safety reasons
		N = 0 # xyz line counter (= number of atoms counter)
		got_ch_mu = False
		for line in inp:
			line.strip()
			dat = line.split()
			if dat != []:
				if len(dat) == 2 and dat[0].lstrip('-').isdigit() and dat[1].lstrip('-').isdigit() and got_ch_mu == False:
					ch = dat[0]
					mu = dat[1]
					got_ch_mu = True
					print('The following charge and multiplicity has been read from ' + xyzfile + ': ' + str(ch) + ' ' + str(mu))
				elif len(dat) == 4:
					if dat[0] in ATOMLISTA.values():
						N += 1
						atomtype = dat[0]
						x = dat[1]
						y = dat[2]
						z = dat[3]
						xyz_lines += '%-2s' % str(atomtype)+str('%28.8f' % float(x))+str('%14.8f' % float(y))+str('%14.8f' % float(z))+"\n"
					else:
						print('Warning: \'' + dat[0] + '\' is not recognized as a valid atom and will be skipped!')
	#print(xyz_lines)
	return xyz_lines, ch, mu

def get_solvent(name,program="gaussian"):
	if program == "gaussian":
		G = {
	      'DAAA': 'Dichloromethane', 'POZ': 'n,n-DiMethylAcetamide', 
	      'PTZ': 'n,n-DiMethylAcetamide', 'PA': 'n,n-DiMethylFormamide', 
	      'PDI': 'Dichloromethane', 'NCE': 'Acetonitrile', 'Me2-Acr': 'n,n-DiMethylFormamide', 
	      'Mes-Acr': 'Acetonitrile', 'BOH-Acr': 'Acetonitrile', 'BF3-Acr': 'Acetonitrile', 
	      'Ph-Acr': 'Acetonitrile', 'FeTP': 'Acetonitrile', 
	      'Ir-': 'n,n-DiMethylFormamide', 'Ni_co': 'n,n-DiMethylFormamide', 
	      'Ru-': 'Acetonitrile', 'Au-alkyne': 'Acetonitrile', 'Cu-phen': 'Dichloromethane', 
	      'Cu-xant': 'Acetonitrile', 'Rh_dimer': 'TetraHydroFuran',
	      'Eos_': 'Acetonitrile', 'CA_': 'Acetonitrile', 'Rh_6G':'Methanol', 'Rh_B':'Ethanol',
	      }
	elif program == "orca":
		# using DMF instead of n,n-DiMethylAcetamide
		G = {
	      'DAAA': 'CH2Cl2', 'POZ': 'DMF', 
	      'PTZ': 'DMF', 'PA': 'DMF', 
	      'PDI': 'CH2Cl2', 'NCE': 'Acetonitrile', 'Me2-Acr': 'DMF', 
	      'Mes-Acr': 'Acetonitrile', 'BOH-Acr': 'Acetonitrile', 'BF3-Acr': 'Acetonitrile', 
	      'Ph-Acr': 'Acetonitrile', 'FeTP': 'Acetonitrile', 
	      'Ir-': 'DMF', 'Ni_co': 'DMF', 
	      'Ru-': 'Acetonitrile', 'Au-alkyne': 'Acetonitrile', 'Cu-phen': 'CH2Cl2', 
	      'Cu-xant': 'Acetonitrile', 'Rh_dimer': 'THF',
	      'Eos_': 'Acetonitrile', 'CA_': 'Acetonitrile', 'Rh_6G':'Methanol', 'Rh_B':'Ethanol',
	      }
	else:
		print("Please provide either gaussian or orca to get solvent.")
	for i in G.keys():
		if i in name:
			return G[i]

def get_genbasis(basisfile):
	with open(basisfile, 'r') as file:
		data = file.read()
	return data

def modify_input(xyzfile, smpl_fullpath, mark1='CH MU', mark2='COORD', mark3='SLVT', mark4='GENBASIS', mark5='XYZFILENAME', basisfile='metal_basis_def2-TZVP.txt'):
	'''
	- Pastes xyz lines from an input into a provided sample input.
	- A new file is created and the sample is kept.
	- The sample file is expected to contain placeholder strings 
	  specified with the 'mark1' and 'mark2' arguments.
	'''
	smpl = os.path.basename(smpl_fullpath)
	if xyzfile == smpl:
		return 0
	if smpl[-4:] == ".inp": #  or smpl[-4:] == ".dat"
		program = "orca"
	else:
		program = "gaussian"
	if '_coords' in xyzfile: # files with such names are obtained with makeinput1()
		iname = xyzfile.replace('_coords','')
		if smpl[-4:] not in iname:
			iname = str(iname)[:-4] + '_' + str(smpl)
	else:
		iname = str(xyzfile)[:-4] + '_' + str(smpl)
	xyz_lines, ch, mu = makeinput2(xyzfile)
	with open(smpl_fullpath) as inp: # parse the sample file
		ifile = open(iname, "w") # create new input file after opening the smpl
		for line in inp:
			if mark1 in line:
				ch_mu = str(ch)+' '+str(mu)
				line = line.replace(mark1,ch_mu)
				if mark5 in line:
					xyz_name = iname[:-4] + ".xyz"
					line = line.replace(mark5,xyz_name)
				ifile.write(line)
# 				ifile.write(str(ch)+' '+str(mu))
# 				ifile.write('\n')
			elif mark2 in line:
				ifile.write(xyz_lines)
				#ifile.write('\n')
			elif mark3 in line:
				solvent = get_solvent(xyzfile,program)
				line = line.replace(mark3,solvent)
				ifile.write(line)
			elif mark4 in line:
				ecp = get_genbasis(basisfile)
				line = line.replace(mark4,ecp)
				ifile.write(line)
			elif mark5 in line:
				xyz_name = iname[:-4] + ".xyz"
				line = line.replace(mark5,xyz_name)
				ifile.write(line)
			else:
				ifile.write(line)
		print(str(iname) + ' has been written!')
		ifile.close()





