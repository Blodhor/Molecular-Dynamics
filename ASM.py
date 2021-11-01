# -*- coding: utf-8 -*-
'''
AMBER molecular dynamics Simulation Manager (works on AMBER18+ and AMBERtools18+). 

Automates the process of performing a molecular dynamics simulation using the Amber/Ambertools computational packages.

Needs AMBER/AMBERtools installed and with its binaries on the computer's PATH. Works on Linux operating systems.

Author: Braga, B. C.
E-mail: bruno.braga@ufms.br
'''

import re
import random
from os import system as cmd
from os import stat
from copy import deepcopy as cp

Titratable_Residue_Names = ['AS4', 'Gl4', 'HIP', 'CYS', 'LYS', 'TYR', 'DAP', 'DCP', 'DG', 'DT', 'AP', 'CP', 'G', 'U', 'HEH', 'PRN']
random.seed()

def IntTransNumb(a):
	'''Translate string number into integer'''

	if '*' not in a and 'E' in a.upper():
		# 10e6
		a_e = a.upper().split('E')
		res = int(a_e[0])**int(a_e[1])
	elif '*' in a and 'E' in a.upper() and '**' not in a:
		# 5*10e6
		temp = a.split('*')
		a_e  = temp[1].upper().split('E')
		res  = int(temp[0])*int(a_e[0])**int(a_e[1])
	elif 'E' not in a.upper() and '**' in a:
		temp0 = a.split('**')
		f = False
		for i in temp0:
			f = f or '*' in i
		if not f:
			# 10**6
			res = int(temp0[0])**int(temp0[1])
		else:
			# 6*10**6
			temp1 = temp0[0].split('*')
			res = int(temp1[0])*int(temp1[1])**int(temp0[1])
	elif 'E' not in a.upper() and '**' not in a and '*' in a:
		# 10*6
		temp = a.split('*')
		res = int(temp[0])*int(temp[1])
	else:
		# 52
		res = int(a)

	return res

class Amber_par:
	'''Creation method of the Amber (or Ambertools) input parameters.'''
	mode  = {'L': (10**-1,1), 'M': (1,10), 'GM': (10,50), 'GH': (50,100), 'GU': (100,500), 'CT': ''}
	steps = 10**5
	AnnealEqFactor = 0.1

	def __init__(self, system='WillThisWork.pdb', ligand='already_docked_positions.pdb', with_docking= False, simulation='CpHMD', pH=7.0, pHstep=1.0, simulation_length='low',
	custom_mode = False, custom_min = 10**4, custom_a = 10**6, custom_e = 10**7, custom_p = 10**8,
	exp_solv= False, solvent_in='water', box_cuttoff= 12.0, information_cycles=200, gpu=0, prot_res='AS4 SER HIP',
	mpi_use= False, mpicomp='mpiexec', mpicores=2, prep_stop= False, chosen_mut_flag = False,
	chosen_mut=['SER_160-MET'], active_site={'SER':[160],'HIS':[237],'ASP':[206]}):
		'''

Parameters
----------
system: pdb file.

-----
Docking options are not complete on version 1.1 and will work properly only on future versions!

ligand: The present code doesn't have a docking module (optmization of the ligand's atoms' positions in the active site of the system[if it is an enzyme]).
		At the present version, it is needed a third party programm to get this PDB. 

with_docking: True/False choice for the creation of a docked ligand-system topology and coordinates file with LEAP.
-----

gpu: Gpu which would be used (if there is two gpu's the choice is between number 0 and 1).

simulation: Goal of the simulation; between a usual Molecular Dynamis (MD) with only a thermal stability analysis or a constant ph dynamics (CpHMD)

pH: Equilibration and Production stages' solvent pH. (Default=7.0)

pHstep: If goal is CpHMD, this variable defines the pH-unit intervals between two productions to be made.

simulation_length: Duration of an average stage's simulation, between 0.1 and 500 (where this factor is multiplied by a default quantity [self.steps] of time steps of 2fs or 4fs with HMR method). Based on the class' dict, self.mode, keys:
	'L': (10**-1,1);
	'M': (1,10);
	'GM': (10,50);
	'GH': (50,100);
	'GU': (100,500); 
	where the second position of the tuple is only for production/CpHMD and the first for annealing and equilibration.
	'CT': custom option.

custom_mode: True/False flag for 'simulation_length' custom.

custom_min: (For simulation_length = CT only) Max number of Minimization steps. 

custom_a: (For simulation_length = CT only) Duration of the Annealing stage.

custom_e: (For simulation_length = CT only) Duration of the Equilibration stage.

custom_p: (For simulation_length = CT only) Duration of the Production stage.

information_cycles: Number of times information regarding energy, rmsd, rmsf, coordinates, etc, are saved (too high of a number requires too much memory, but too few cycles can't help in the system properties' analysis). Directly proportional to the number of frames.

exp_solv: True/False choice for explicit solvent.

solvent_in: Choice of explicit solvent used. Solvents available: water; methanol; chloroform; N-methyacetamide; and urea.

box_cuttoff: If you are explicitly solvating, this flag indicates the cuttoff in angstrom for the solvate box (nonbonded cutoff in minimization, annealing...).

prot_res: Residues to tritate, for CpHMD.

mpi_use: Like the name implies, it's a flag that decides if we are using MPI or not.

mpicomp: MPI compiler like mpiexec/mpirun. [OpenRTE 2.1.1/Open MPI 2.1.1]

mpicores: (Integer) number of cores to use.

prep_stop: If 'True', the manager will stop after preparing all inputs and the shell script file which starts the minimization and all the rest of the simulation. If 'False', the simulation runs completly. 

--------------------------------------------------------------------------------------
| The following is only for Amber_mutation class (a "grandchild" class of Amber_par) |
--------------------------------------------------------------------------------------

active_site: Dict. of aminoacids from the active site (if the system is an enzyme). Ex:{residueA: residueA_id, residueB: residueB_id,...}. Option needed for the random active site mutation. 

chosen_mut_flag: (Boolean) States if a specific mutation will be made (one or more residue changes).

chosen_mut: Specific mutation chosen.'''

		cmd('rm hmr-in.sh simulation.sh rmsdf.cpptraj BeforeManager.out AfterManager.out')
		cmd('ls > BeforeManager.out')

		# Creating object's attributes
		bool_class           = 'Amber_mutation' in str(type(self)) 
		if bool_class:
			self.exp_solv    = True
		else:
			self.exp_solv    = exp_solv
		self.pdb             = system
		self.goal            = simulation
		self.pH              = pH
		self.hmr             = True
		self.gpunumber       = gpu
		self.cph             = False
		self.simulation_mode = simulation_length
		self.mode_custom     = custom_mode
		self.mode_custom_min = custom_min
		self.mode_custom_a   = custom_a
		self.mode_custom_e   = custom_e
		self.mode_custom_p   = custom_p
		if self.simulation_mode == 'L':
			self.AnnealEqFactor  = 1
		self.prmtop          = 'system.prmtop'
		self.rms             = 'system.rms'
		self.rmsf            = 'system_rmsf.agr'
		self.solvent         = solvent_in.lower()
		self.cuttoff         = float(box_cuttoff)
		self.mpi_use         = mpi_use
		self.mpi             = mpicomp
		self.cores           = int(mpicores)
		self.prep_stop       = prep_stop
		self.pH_step         = pHstep
		
		# Sorting protonation and active site information for later verification
		prot_ = prot_res.split()
		for i in range(len(prot_)):
			if prot_[i] == 'AS4':
				prot_[i] = 'ASP'
			elif prot_[i] == 'HIP':
				prot_[i] = 'HIS'
			elif prot_[i] == 'GL4':
				prot_[i] = 'GLU'
		site_ = []
		for i in active_site.keys():
			site_.append(i)
		prot_.sort()
		site_.sort()
		self.ok = True
		
		# Docking-only options
		self.docking    = with_docking
		self.ligand     = ligand
		self.lig_mol    = ligand[:-3]+"mol2"
		self.lig_frcmod = ligand[:-3]+"frcmod"
		self.lig_name   = "DOK"

		if self.docking:
			f = open(self.ligand,'r+')
			temp = f.readlines()
		
			f.seek(0)
			for i in range(len(temp)):
				j = temp[i]
				data = j.split()
				if 'ATOM' in data[0] or 'ANISOU' in data[0] or 'HETATM' in data[0]:
					if len(data[2]) == 1 and data[2] != data[len(data)-1]:
						#HETATM    1  C   UNL     1     -27.682  -2.375   1.924  1.00  0.00      .068 A\n
						bit = j.find(data[len(data)-1],len(data[0])+1,len(j))
						j = j[:bit]
						temp[i] = j + data[2]
					if len(data[3]) == 3:
						bit2 = j.find(data[3],len(data[0])+1,len(j))
						jj = j[:bit2]+self.lig_name+j[bit2 +len(self.lig_name):]
						temp[i] = jj
						#HETATM    1  C   DOK     1     -27.682  -2.375   1.924  1.00  0.00      .068 A\n
					else:
						self.lig_name = data[3]

					if i == len(temp) -1:
						f.write(temp[i])
					else:
						f.write(temp[i]+'\n')
			f.close()
			cmd('reduce %s > ligand-reduced.pdb'%self.ligand)
			cmd('antechamber -fi pdb -fo mol2 -i ligand-reduced.pdb -o %s -c bcc -pf y -nc 0'%self.lig_mol)
			#
			# Docking related only:
			#  There is a problem with "DUN" and "hn" atom types on the mol2 file
			#  Needs to be fixed on later versions! 
			#
			cmd('parmchk2 -i %s -o %s -f mol2'%(self.lig_mol,self.lig_frcmod))

		# Mutation-only parameters
		self.pdb_string          = ''
		self.pdb_str_2Mut        = ''
		self.pdb_str_not2Mut     = ''
		self.mut_pos             = {}
		self.active_site         = active_site
		for i in self.active_site:
			self.active_site[i].sort()
		self.beforeAmber_res2mut = {}
		self.not_rand_mut        = chosen_mut_flag
		self.oldid_mut           = {}
		self.chosen_mutation     = {}
		self.chosen_residues     = {}
		self.mut_count           = 0
		self.newid_info          = []
		res_newid                = ''
		resmut_newid             = {}
		mutid                    = {}

		# Adjusting (mutation) object's attributes
		if bool_class:
			if chosen_mut_flag:
				# chosen_mut := ['SER_160-MET']
				for ji in chosen_mut:
					pos_ji = ji.find('-',0,len(ji))
					# ji[:pos_ji] := 'SER_160'
					pos_resid = ji.find('_',0,pos_ji)
					# ji[:pos_resid] := 'SER'; ji[pos_resid+1:pos_ji] := '160'; ji[pos_ji+1:] := 'MET'
					self.chosen_mutation[(ji[:pos_resid], ji[pos_resid+1:pos_ji])] = ji[pos_ji+1:]
					# self.chosen_mutation := {('SER', '160'): 'MET', ('SER', '125'): 'MET'}
					if ji[:pos_resid] not in mutid:
						mutid[ji[:pos_resid]] = [ int(ji[pos_resid+1:pos_ji]) ]
					else:
						mutid[ji[:pos_resid]].append( int(ji[pos_resid+1:pos_ji]) )
					# mutid := {'SER': [160,125]}
				for ji in mutid:
					mutid[ji].sort()
				# mutid := {'SER': [125,160]}
				self.beforeAmber_res2mut = cp(mutid)
				# resmut_newid := {'ARG': [ [34, ['coords1','coords2',...]], [53, ['coords1','coords2',...]], ... ], ... }
				resmut_newid = self.adjusting_firstep(chosen_res=mutid)
			else:
				self.beforeAmber_res2mut = cp(self.active_site)
				for i in prot_:
					if i not in site_ and input("Residue %s not in active site. Do you wanna procced? (Please answer only with [Yes/No])\n"%i).lower() == "no":
						self.ok = False

			# The method adjusting_firstep is written on the Amber_mutation class
			res_newid = self.adjusting_firstep(chosen_res=self.active_site)

		# res_newid now have all coords of the active_site atoms, and with this, it's trivial to look for the new active site atoms' id (after pdb4amber)!
		if self.ok:
			self.protonation_res = prot_
			bool_pdb = self.pdb != 'WillThisWork.pdb'
			if bool_pdb and self.goal.lower() == 'cphmd':
				self.cph= True
				cmd('pdb4amber -i %s -o system-amber.pdb --constantph --dry'%self.pdb)
				self.pdb ='system-amber.pdb'
			elif bool_pdb:
				cmd('pdb4amber -i %s -o system-amber.pdb --dry'%self.pdb)
				self.pdb ='system-amber.pdb'

			if bool_class:
				if chosen_mut_flag:
					self.oldid_mut = cp(self.chosen_mutation)
					new_ambermutid = {}
					for i in resmut_newid:
						if i == 'ASP':
							new_ambermutid['AS4'] = resmut_newid['ASP']
						elif i == 'HIS':
							new_ambermutid['HIP'] = resmut_newid['HIS']
						elif i == 'GLU':
							new_ambermutid['GL4'] = resmut_newid['GLU']
						else:
							new_ambermutid[i] = resmut_newid[i]
					resmut_newid = new_ambermutid
					print("Old residue id:\n", mutid)
					# self.chosen_mutation := {('SER', '160'): 'MET', ('SER', '125'): 'MET'}
					self.chosen_residues = self.adjusting_finalstep(coords=resmut_newid)
					# self.chosen_mutation := {('SER', 'newid'): 'MET', ('SER', 'O1newid'): 'MET'}
					if self.chosen_residues != -1:
						print("New residue id:\n", self.chosen_residues)
				else:
					new_amber_site = {}
					for i in res_newid:
						if i == 'ASP':
							new_amber_site['AS4'] = res_newid['ASP']
						elif i == 'HIS':
							new_amber_site['HIP'] = res_newid['HIS']
						elif i == 'GLU':
							new_amber_site['GL4'] = res_newid['GLU']
						else:
							new_amber_site[i] = res_newid[i]
					res_newid = new_amber_site
					print("Old active site's id:\n",self.active_site)
					self.active_site = self.adjusting_finalstep(coords=res_newid)
					# self.active_site := {'SER':[newid, ...], ...}
					print("New active site's id:\n",self.active_site)

			if information_cycles <50:
				print('The number of information cycles is too short...ajusting to default value.') 
				self.info_factor = int(100)
			elif information_cycles > 1000:
				print('The number of information cycles is too high...ajusting to the upper limit.')
				self.info_factor = int(1000)
			else:
				self.info_factor = int(information_cycles)

	def hmr_transform_file(self, input_name='hmr-prmtop.in', prmtop_new='systemHMR.prmtop'):
		'''Creates the Parmed input file for a HMR transform in the topology file.

Parameters
----------

input_name: Name for the Parmed input.

prmtop_new: Name chosen for the new topoly file.'''

		cmd('rm %s'%input_name)
		f = open(input_name,'w')
		f.write('HMassRepartition\n')
		f.write('outparm %s\nquit'%prmtop_new)
		f.close()
		self.prmtop = prmtop_new

	def leap_in(self, checking=False, checking_cmd='charge', charge_fix=False, add_Na=0, add_Cl=0,
	antechamber=False, mol2='system.mol2', frcmod='system.frcmod'):
		'''Based on what kind of simulation, creates a leap.in file in order to use implicit (CpHMD) or explicit (MD) solvent and generates the topology.

Parameters
----------

checking: Decides if the code will check for abnormalities in the system before it tries to create the topoly file.

checking_cmd: Leap commands: 'charge' or 'check'.

charge_fix: Decides if the arguments 'add_Na' and 'add_Cl' will be used. OBS: Only possible for explicit solvent.

add_Na: Adds a positive integer charge in the system.

add_Cl: Adds a negative integer charge in the system.

antechamber: If True means leap library can't deal with your system and you chose to try creating the mol2 file with antechamber.

mol2: File created with antechamber.

frcmod: File created with parmck2.'''

		cmd('rm tleap.in')
		# CpHMD or not...it doesn't make a difference in leap if you just use all libraries in both
		# The flag 'w' on the open() function -> if there is a file 'tleap.in', it'll be overwritten
		f = open('tleap.in','w')
		f.write('source leaprc.gaff2\n')
		f.write('source leaprc.DNA.OL15\n')
		f.write('source leaprc.lipid17\n')
		f.write('source leaprc.water.tip3p\n')
		f.write('source leaprc.protein.ff14SB\n')
		f.write('source leaprc.constph\n')
		f.write('source leaprc.conste\n')
		if antechamber:
			f.write('SYS = loadmol2 %s\n'%mol2)
			f.write('loadamberparams %s\n'%frcmod)
		else:
			f.write('SYS = loadPDB %s\n'%(self.pdb))

		if self.docking:
			f.write('%s = loadmol2 %s\n'%(self.lig_name, self.lig_mol))
			f.write('NEW = combine {SYS %s}\n'%self.lig_name)
			f.write('loadamberparams %s\n'%self.lig_frcmod)
			sys_temp = 'NEW'
		else:
			sys_temp = 'SYS'
		
		# The two commands create periodic solvent boxes around the solute, which should be a UNIT. 
		#  solvateBox creates a cuboid box, while solvateOct creates a truncated octahedron.
		if self.exp_solv: 
			f.write('loadoff solvents.lib\n')
			if self.solvent == 'water':
				f.write('solvateOct %s TIP3PBOX %.2f\n'%(sys_temp, self.cuttoff))
			elif self.solvent == 'methanol':
				f.write('loadamberparams frcmod.meoh\n')
				f.write('solvateOct %s MEOHBOX %.2f\n'%(sys_temp, self.cuttoff))
			elif self.solvent == 'chloroform':
				f.write('loadamberparams frcmod.chcl3\n')
				f.write('solvateOct %s CHCL3BOX %.2f\n'%(sys_temp, self.cuttoff))
			elif self.solvent == 'n-methyacetamide':
				f.write('loadamberparams frcmod.nma\n')
				f.write('solvateOct %s NMABOX %.2f\n'%(sys_temp, self.cuttoff))
			elif self.solvent == 'urea':
				f.write('loadoff 8Mureabox.off\n')
				f.write('loadamberparams frcmod.urea\n')
				f.write('solvateOct %s UREABOX %.2f\n'%(sys_temp, self.cuttoff))

			if charge_fix:
				f.write('addIonsRand %s Cl- %d Na+ %d\n'%(sys_temp, add_Cl, add_Na))

		if not checking:
			f.write('savepdb %s systemLeap.pdb\n'%sys_temp)
			f.write('saveAmberParm %s %s systemLeap.rst7\n'%(sys_temp, self.prmtop))
		else:
			f.write('%s %s\n'%(checking_cmd, sys_temp))
		f.write('quit')
		f.close()

	def input_min(self, min_name='minimization'):
		'''Creates the input file for energy minimization.

Parameter
---------
min_name: Name for the minimization input. 
	OBS: It's useful to set the names like this or using object/class attributes to name every related file the same (only with different extensions), which makes it easier for the user later. 
'''
		cmd('rm %s.in'%min_name)
		f = open('%s.in'%min_name,'w')
		f.write('Energy Minimization Stage\n')
		f.write('&cntrl\n')
		# imin=1 -> Flag to run minimization.
		f.write(' imin=1,\n')
		# The maximum number of cycles of minimization.
		if not self.mode_custom:
			f.write(' maxcyc=10000,\n')
		else:
			f.write(' maxcyc=%d,\n'%self.mode_custom_min)
			change = 2000
			if self.mode_custom_min <= 2000:
				change = int(0.8*self.mode_custom_min)
		# If ntmin=1 then the method of minimization will be switched from steepest descent to conjugate gradient after ncyc cycles.
		f.write(' ntmin=1,\n')
		f.write(' ncyc=%d,\n'%change)
		# Every ntpr steps, energy information will be printed in human-readable form to files "mdout"
		f.write(' ntpr=250,\n')
		if not self.exp_solv:
			# Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models. For more info into igb, go to page 61 of Amber18's manual.
			f.write(' igb=2,\n')
			# ntb=0 no periodicity is applied and PME is off (default when igb > 0)
			f.write(' ntb=0,\n')
			# No cutoff. When igb > 0, the default is 9999.0 (effectively infinite)
			f.write(' cut=9999.0,\n')
		else:
			# No generalized Born term is used. (Default).
			f.write(' igb=0,\n')
			# This variable controls whether or not periodic boundaries are imposed on the system during the calculation of non-bonded interactions. ntb=1 constant volume;
			f.write(' ntb=1,\n')
			# This is used to specify the nonbonded cutoff, in Angstroms.
			f.write(' cut=%.2f,\n'%self.cuttoff)
			
		if self.cph:
			# Turn on constant pH 
			if not self.exp_solv:
				# You must set icnstph=1 to turn on constant pH in implicit solvent.
				f.write(' icnstph=1,\n')
				# Use the salt conc. CpHMD was parametrized for
				f.write(' saltcon=0.1,\n')
				# Never attempt to change prot. states (way bigger number of steps than the minimization itself)
				f.write(' ntcnstph=100000,\n')
				# Turn on positional restraints. Flag for restraining specified atoms in Cartesian space using a harmonic potential, if ntr > 0. The restrained atoms are determined by the restraintmask string. The force constant is given by restraint_wt. The coordinates are read in "restrt" format from the "refc" file.
				f.write(' ntr=1,\n')
			f.write(' restraint_wt=10,\n')
			# String that specifies the restrained atoms when ntr=1. '@CA,C,O,N'Restraints on the backbone atoms only.
			f.write(' restraintmask=\'!:WAT&@CA,C,O,N\',\n')
			
		f.write('/\n')
		f.close()

	def input_heat(self, annealing_name='annealing', mdsteps_factor=1, custom=False):
		'''Creates the input file for annealing.

Parameters
----------

annealing_name: Name for the annealing stage input.

mdsteps_factor: The factor by which the attibute self.steps will be multiplied (defining the simulation's length).
'''

		cmd('rm %s.in'%annealing_name)
		f = open('%s.in'%annealing_name,'w')
		f.write('Annealing\n&cntrl\n')
		# Switch for temperature scaling. Setting ntt=0 corresponds to the microcanonical (NVE) ensemble. ntt=1 Constant temperature, using the weak-coupling algorithm. ntt=2 Andersen-like temperature coupling scheme, in which imaginary "collisions" randomize the velocities to a distribution corresponding to temp0 every vrand steps (in between these "massive collisions", the dynamics is Newtonian). ntt=3 Use Langevin dynamics with the collision frequency given by gamma_ln
		f.write(' ntt=3,\n')
		# Option to read the initial coordinates, velocities and box size from the inpcrd file. = 1 Coordinates, but no velocities, will be read;
		f.write(' ntx=1,\n')
		# imin=0  -> Run molecular dynamics without any minimization.
		f.write(' imin=0,\n')
		# Flag to restart a simulation. = 0 run as a new simulation.
		f.write(' irest=0,\n')
		# Number of MD-steps to be performed.
		if not self.mode_custom:
			f.write(' nstlim=%d,\n'%(mdsteps_factor*self.steps))
		else:
			f.write(' nstlim=%d,\n'%self.mode_custom_a)
		# The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if it isn't. The use of Hydrogen Mass Repartitioning (HMR) together with SHAKE, allows the time step to be increased in a stable fashion by about a factor of two (up to .004) by slowing down the high frequency hydrogen motion in the system. To use HMR, the masses in the topology file need to be altered before starting the simulation. ParmEd can do this automatically with the HMassRepartition option.
		if self.hmr:
			f.write(' dt=0.004,\n')
		else:
			f.write(' dt=0.002,\n')

		# Every ntpr steps, energy information will be printed in human-readable form to files "mdout"
		f.write(' ntpr=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		# Every ntwx steps, the coordinates will be written to the mdcrd file.
		f.write(' ntwx=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		# Every ntwr steps during dynamics, the “restrt” file will be written, ensuring that recovery from a crash will not be so painful.
		f.write(' ntwr=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		# Flag for constant pressure dynamics. = 0 No pressure scaling.
		f.write(' ntp=0,\n')
		# Flag used to control which barostat to use in order to control the pressure. = 1 Berendsen barostat; = 2 Monte Carlo barostat.
		f.write(' barostat=2,\n')
		# Reference pressure (in units of bars, where 1 bar ≈ 0.987 atm) at which the system is maintained
		f.write(' pres0=1.0,\n')
		# Pressure relaxation time (in ps), when ntp>0. The recommended value is between 1.0 and 5.0 psec.
		f.write(' taup=2.0,\n')
		# Flag for SHAKE to perform bond length constraints. = 1 SHAKE is not performed; = 2 bonds involving hydrogen are constrained; = 3 all bonds are constrained.
		f.write(' ntc=2,\n')
		# Typically ntf=ntc. For water models, a special "three-point" algorithm is used. Consequently, to employ TIP3P set ntf=ntc=2.
		f.write(' ntf=2,\n')
		# Initial temperature.
		f.write(' tempi=10.0,\n')
		# Reference temperature at which the system is to be kept, if ntt > 0.
		f.write(' temp0=300,\n')
		# The seed for the pseudo-random number generator. If ig=-1, the random seed will be based on the current date and time, and hence will be different for every run.
		f.write(' ig=-1,\n')
		if not self.exp_solv:
			# This is used to specify the nonbonded cutoff, in Angstroms. For PME, the cutoff is used to limit direct space sum. When igb > 0, the default is 9999.0 (effectively infinite).
			f.write(' cut=9999.0,\n')
			# Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models.
			f.write(' igb=2,\n')
			# ntb=0 no periodicity is applied and PME is off (default when igb > 0)
			f.write(' ntb=0,\n')
		else:
			# No generalized Born term is used. (Default).
			f.write(' igb=0,\n')
			# This variable controls whether or not periodic boundaries are imposed on the system during the calculation of non-bonded interactions. ntb=1 constant volume;
			f.write(' ntb=1,\n')
			# This is used to specify the nonbonded cutoff, in Angstroms.
			f.write(' cut=%.2f,\n'%self.cuttoff)

		if not self.cph:
			# The collision frequency, in ps**-1 , when ntt = 3.
			f.write(' gamma_ln=3,\n')
			# nmropt=0, No nmr-type analysis will be done.
			f.write(' nmropt=0,\n')
		else:
			# The collision frequency, in ps**-1 , when ntt = 3.
			f.write(' gamma_ln=5.0,\n')
			# Turn on constant pH
			if not self.exp_solv:
				# You must set icnstph=1 to turn on constant pH in implicit solvent.
				f.write(' icnstph=1,\n')
				# Never attempt to change prot. states
				f.write(' ntcnstph=100000000,\n')
				# tol: Relative geometrical tolerance for coordinate resetting in shake. Recommended maximum: <0.00005. Angstrom Default 0.00001.
				f.write(' tol=0.000001,\n')
				# NRESPA: This variable allows the user to evaluate slowly-varying terms in the force field less frequently. For PME, "slowly-varying" (now) means the reciprocal sum. For generalized Born runs, the "slowly-varying" forces are those involving derivatives with respect to the effective radii, and pair interactions whose distances are greater than the "inner" cutoff, currently hard-wired at 8 Å. If NRESPA>1 these slowly-varying forces are evaluated every nrespa steps. The forces are adjusted appropriately, leading to an impulse at that step. If nrespa*dt is less than or equal to 4 fs the energy conservation is not seriously compromised. However if nrespa*dt > 4 fs the simulation becomes less stable. Note that energies and related quantities are only accessible every nrespa steps, since the values at other times are meaningless.
				f.write(' nrespa=1,\n')
				# Use the salt conc. CpHMD was parametrized for
				f.write(' saltcon=0.1,\n')
				# if ntr > 0 The restrained atoms are determined by the restraintmask string.
				f.write(' ntr=1,\n')
				# The restraint weight (2 kcal/mol/Ang^2) is not a weak restraint (nor is it particularly strong). The structure should not change conformation very much with this setting.
				f.write(' restraint_wt=2.0,\n')
				f.write(' restraintmask=\'@CA,C,O,N\',\n')

			f.write(' ioutfm=1,\n')
			# nmropt= 1 NMR restraints and weight changes will be read.
			f.write(' nmropt=1,\n')

		f.write('/\n\n&wt')
		if self.cph and self.exp_solv:
			f.write(' TYPE=\'TEMP0\', ISTEP1=0, ISTEP2=150000,\n')
			f.write(' VALUE1=10.0, VALUE2=300.0\n')
		else:
			if mdsteps_factor*self.steps >= 1000000:
				f.write(' TYPE=\'TEMP0\', ISTEP1=1, ISTEP2=500000,\n')
			else:
				f.write(' TYPE=\'TEMP0\', ISTEP1=1, ISTEP2=%d,\n'%(mdsteps_factor*self.steps/2))
			f.write(' VALUE1=10.0, VALUE2=300.0\n')

		f.write('/\n\n')
		f.write('&wt TYPE=\'END\'\n')
		f.write('\n/\n')
		f.close()

	def input_equil(self, step='Equilibration', mdsteps_factor=1, ph=7.0, trescnt=11):
		'''Creates the input files for equilibration or the production in CpHMD.
OBS: The equilibration stage is the part of the setup where we allow the structure to stabilize to its surroundings and thermodynamic constraints (e.g., constant temperature and/or pressure). In this stage, we will finally allow our protonation states to change via Monte Carlo moves (if we are going for CpHMD).

Parameters
----------

step: Name for the current stage input (either equilibration or production for CpHMD).

mdsteps_factor: The factor by which the attibute self.steps will be multiplied (defining the simulation's length.

ph: Chosen ph. It'll be used only integer values in this code (keep in mind it was just a small design choice, Amber allows it to be a double/float).

trescnt: The number of residues to titrate (the methods on this class defined, will look for this info automatically if needed).'''

		cmd('rm %s.in'%step)
		f = open('%s.in'%step,'w')
		f.write('%s\n'%step)
		f.write('&cntrl\n')
		f.write(' ntt=3,\n')
		# ntx=5 Coordinates and velocities will be read from either a NetCDF or a formatted (ASCII) coordinate file. The velocity information will only be used if irest=1
		f.write(' ntx=5,\n')
		f.write(' imin=0,\n')
		f.write(' irest=1,\n') 
		if not self.mode_custom:
			f.write(' nstlim=%d,\n'%(mdsteps_factor*self.steps))
		else:
			if 'CpHMD_prod' in step:
				f.write(' nstlim=%d,\n'%self.mode_custom_p)
			else:
				f.write(' nstlim=%d,\n'%self.mode_custom_e)
		if self.hmr:
			f.write(' dt=0.004,\n')
			# ntcnstph times the number of residues (TRESCNT in the cpin file) times our 'dt' (in ps) should be ~100 fs = 0.1 ps if possible.# ntcnstph = 100/ (dt*1000*TRESCNT)
			ntph = int(100/(4*trescnt))+1
		else:
			f.write(' dt=0.002,\n')
			ntph = int(100/(2*trescnt))+1
		f.write(' ntpr=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		f.write(' ntwx=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		f.write(' ntwr=%d,\n'%(int(mdsteps_factor*self.steps/self.info_factor)))
		f.write(' ntc=2,\n')
		f.write(' ntf=2,\n')
		if not (self.exp_solv and self.cph):
			f.write(' pres0=1.0,\n')
			f.write(' barostat=2,\n')
			# Pressure relaxation time (in ps), when NTP > 0. The recommended value is between 1.0 and 5.0 psec. Default value is 1.0
			f.write(' taup=1.0,\n')
			
		if not self.exp_solv:
			# This is used to specify the nonbonded cutoff, in Angstroms. For PME, the cutoff is used to limit direct space sum. When igb > 0, the default is 9999.0 (effectively infinite).
			f.write(' cut=9999.0,\n')
			# Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models.
			f.write(' igb=2,\n')
			# ntb=0 no periodicity is applied and PME is off (default when igb > 0)
			f.write(' ntb=0,\n')
			# Flag for constant pressure dynamics. This option should be set to 1 or 2 when Constant Pressure periodic boundary conditions are used.
			f.write(' ntp=0,\n')
		else:
			# No generalized Born term is used. (Default).
			f.write(' igb=0,\n')
			# This is used to specify the nonbonded cutoff, in Angstroms.
			f.write(' cut=%.2f,\n'%self.cuttoff)
			if 'CpHMD_prod' in step:
				# Cphmd explicit solvent (PRODUCTION)
				f.write(' ntb=1,\n') # default value when igb=ntp=0
				f.write(' ntp=0,\n') # default value
				# Format of the final 'restrt' file
				f.write(' ntxo=2,\n')# default value (NetCDF)
				# icnstph=2, indicates that CpHMD should be run in explicit solvent.
				f.write(' icnstph=2,\n')
				# Use the salt conc. CpHMD was parametrized for
				f.write(' saltcon=0.1,\n')
				f.write(' solvph=%.2f,\n'%ph)
				# Good results are expected with ntcnstph set between 100 and 500 steps (dt=0.002).
				f.write(' ntcnstph=100,\n') ## %d,\n'%ntph
				# ntrelax: This is the number of steps of relaxation dynamics you wish to run. 200 fs is sufficient to account for the bulk of the solvent relaxation.(This is unique to running in explicit solvent)
				f.write(' ntrelax=100,\n')
			else:
				# ntb=2 constant pressure
				f.write(' ntb=2,\n')
				# ntp=1 MD with isotropic position scaling
				f.write(' ntp=1,\n')
			
		if not self.cph:
			# The collision frequency can be doubled in the equilibration and production stages
			f.write(' gamma_ln=6,\n')
		else:
			if self.exp_solv and 'CpHMD_prod' not in step:
				f.write(' gamma_ln=1.0,\n')
			else:
				# The collision frequency can be doubled in the equilibration and production stages
				f.write(' gamma_ln=5,\n')
				
		# The format of coordinate and velocity trajectory files (mdcrd, mdvel and inptraj). ioutfm=1 -> Binary NetCDF trajectory
		f.write(' ioutfm=1,\n')
		f.write(' temp0=300.0,\n')
		f.write(' tempi=300.0,\n')
		f.write(' nmropt=0,\n') # default
		f.write(' ig=-1,\n')
		f.write('/\n')
		f.close()

	def prod_cph(self, mdsteps_factor=10, ph=[7.0,7.0], trescnt=11):
		'''Creates the input files for production or cphmd. For cphmd we will run simulations at pH values of ph[0] through ph[1] with 1 pH-unit intervals (the template input file for the cphmd production dynamics is the same as the equilibration stage).

Parameters
----------

mdsteps_factor: The factor by which the attibute self.steps will be multiplied (defining the simulation's length.

ph: Range of ph chosen.

trescnt: The number of residues to titrate (the methods on this class defined, will look for this info automatically if needed).'''
		ph_i  = ph[0]
		nm_ph = []
		while ph_i <= ph[1]:
			nm_ph.append( ('CpHMD_prod_%.2f'%ph_i,ph_i) )
			ph_i += self.pH_step

		if not self.cph:
			for pp in nm_ph:
				cmd('rm Production_%.2f.in'%pp[1])
				f = open('Production_%.2f.in'%pp[1],'w')
				f.write('Production_pH%.2f\n&cntrl\n'%pp[1])
				f.write(' ntt=3,\n')
				f.write(' gamma_ln=6,\n')
				f.write(' ntx=5,\n')
				f.write(' imin=0,\n')
				f.write(' irest=1,\n')
				if not self.mode_custom:
					f.write(' nstlim=%d,\n'%(mdsteps_factor*self.steps))
				else:
					f.write(' nstlim=%d,\n'%self.mode_custom_p)
				if self.hmr:
					f.write(' dt=0.004,\n')
				else:
					f.write(' dt=0.002,\n')

				f.write(' ntpr=%d,\n'%(int(mdsteps_factor*self.steps/(2*self.info_factor))))
				f.write(' ntwx=%d,\n'%(int(mdsteps_factor*self.steps/(2*self.info_factor))))
				f.write(' ntwr=%d,\n'%(int(mdsteps_factor*self.steps/(2*self.info_factor))))
				f.write(' ntp=0,\n')
				f.write(' barostat=2,\n')
				f.write(' pres0=1.0,\n')
				f.write(' taup=1.0,\n')
				if not self.exp_solv:
					# This is used to specify the nonbonded cutoff, in Angstroms. For PME, the cutoff is used to limit direct space sum. When igb > 0, the default is 9999.0 (effectively infinite).
					f.write(' cut=9999.0,\n')
					# Flag for using the generalized Born or Poisson-Boltzmann implicit solvent models.
					f.write(' igb=2,\n')
					# ntb=0 no periodicity is applied and PME is off (default when igb > 0)
					f.write(' ntb=0,\n')
				else:
					# No generalized Born term is used. (Default).
					f.write(' igb=0,\n')
					# MD-explicit solvent pH
					f.write(' solvph=%.2f,\n'%pp[1])
					# This variable controls whether or not periodic boundaries are imposed on the system during the calculation of non-bonded interactions. ntb=1 constant volume;
					f.write(' ntb=1,\n')
					# This is used to specify the nonbonded cutoff, in Angstroms.
					f.write(' cut=%.2f,\n'%self.cuttoff)

				f.write(' ntc=2,\n')
				f.write(' ntf=2,\n')
				f.write(' ioutfm=1,\n')
				f.write(' temp0=300.0,\n')
				f.write(' nmropt=0,\n')
				f.write(' ig=-1,\n/\n')
				f.close()
		else:
			for pointer in nm_ph:
				self.input_equil(step=pointer[0], mdsteps_factor=mdsteps_factor, ph=pointer[1], trescnt=trescnt)	

	def cpptraj_in(self, input_name='rmsdf.cpptraj', rms_file='system.rms', rmsf_file='system_rmsf.agr',
	min_name='minimization', annealing_name='annealing', equil_name='equilibration', phr=[0.0,7.0]):
		'''Creates the input file for trajectory analysis.

Parameters
----------

input_name: Cpptraj input file's name.

rms_file: RMS file's name.

rmsf_file: RMSF data file's name.

min_name: Minimization's files name.

annealing_name: Annealing's files name.

equil_name: Equilibration's files name.

phr: Ph range used to titrate.'''
		ph_i = phr[0]
		pH   = []
		while ph_i <= phr[1]:
			pH.append( ph_i )
			ph_i += self.pH_step

		if not self.cph:
			for pp in pH:
				cmd('rm MD_%.2f_%s'%(pp,input_name))
				f = open('MD_%.2f_%s'%(pp,input_name),'w')
				f.write('parm %s\n'%self.prmtop)
				f.write('trajin %s.nc\n'%annealing_name)
				f.write('trajin %s.nc\n'%equil_name)
				f.write('trajin Production_%.2f.nc\n'%pp)
				f.write('reference %s.rst7\n'%min_name)
				f.write('autoimage\n')
				# See pg. 684 Amer18 guide for rms option, pg. 682 for radgyr and pg. 638 for atomicfluct
				f.write('rms reference @CA,C,O,N,H&!(:WAT) mass out MD_%.2f_%s\n'%(pp,rms_file))
				# @ => atom  # : => residue
				# @CA,C,O,N,H&!(:WAT) := means all these atoms of backbone and not in water
				f.write('atomicfluct out MD_%.2f_%s @CA,C,O,N,H&!(:WAT)\n'%(pp,rmsf_file))
				f.close()
		else:
			for pp in pH:
				cmd('rm rmsdf_CpHMD_%.2f.cpptraj'%pp)
				f = open('rmsdf_CpHMD_%.2f.cpptraj'%pp,'w')
				f.write('parm %s\n'%self.prmtop)
				f.write('trajin CpHMD_prod_%.2f.nc\n'%pp)
				f.write('reference %s.rst7\n'%equil_name)
				# For future reference:
				#  distance dPET-sitio :PET :132,178,209 out dPET-sitio.dat #pg 649-650
				f.write('rmsd reference @CA,C,O,N,H&!(:WAT) first out ph%.2f_rmsd.dat mass\n'%(pp))
				f.write('atomicfluct out ph%.2f_rmsf.dat @CA,C,O,N,H&!(:WAT) byres\n'%(pp))
				f.write('radgyr ph%.2f-radgyr out ph%.2f_radgyr.dat @CA,C,O,N,H&!(:WAT) mass nomax\n'%(pp, pp))
				f.close()

class Amber_run(Amber_par):
	'''Shell script manager (creator of the simulation's executable)'''
	autor       = 'Braga, B. C.'
	email       = 'bruno.braga@ufms.br'
	linebreak   = '#'*60
	report_file = 'report_temp.py' 
	
	def leap_exec(self, mol2_check=False):
		'''Run tleap to prepare the system for Amber simulations.

Returns -2 if something is strange with the system charge (for explicit solvent only).
Returns -3 if there is any FATAL error in trying to build the topology.
Returns 0 if the topology was succesfully built.

Parameter
---------

mol2_check: If False, load system pdb on leap; if True, it'll run antechamber to get mol2 and then try leap.'''

		if not mol2_check:
			self.leap_in()
		else:
			# Leap library doesn't know your system so we'll try to create the mol2 file and load it on Leap
			cmd('reduce %s > system-amber-reduced.pdb'%self.pdb)
			cmd('antechamber -fi pdb -fo mol2 -i system-amber-reduced.pdb -o system-amber.mol2 -c bcc -pf y -nc 0')
			cmd('parmchk2 -i system-amber.mol2 -o system-amber.frcmod -f mol2')
			self.leap_in(antechamber=True, mol2='system-amber.mol2', frcmod='system-amber.frcmod')

		cmd('rm leap.log')
		cmd('tleap -f tleap.in')
		if self.exp_solv:
			tt = open('leap.log','r')
			temp = tt.readlines()
			tt.close()
			charge = 0.0
			for i in range(len(temp)):
				if 'unperturbed charge' in temp[-i]:
					data = temp[-i][len('The unperturbed charge of unit'):].split()
					charge = float(data[1][1:-1])
					break
			cha = int(charge)
			if charge - cha > 0.05:
				print("\nError with the charge\n")
				return -2
			elif cha > 0:
				self.leap_in(charge_fix=True, add_Na=3, add_Cl=cha+3)
			elif cha < 0:
				self.leap_in(charge_fix=True, add_Na=3-cha, add_Cl=3)

			# This will check for any fatal errors, proximity warnings don't matter (minimization can fix it).
			#  In order to check for true errors we need to delete all leap log before attempting the leapfix
			cmd('rm leap.log')
			cmd('tleap -f tleap.in')

		# In the check system, we are looking for something troubling like this:
		# FATAL:  Atom .R<ALA 100>.A<HG2 11> does not have a type.
		# With this kind of error it's impossible to build the topology 		
		f = open('leap.log','r')
		temp = f.readlines()
		error = []
		for i in temp:
			if 'FATAL:' in i:
				error.append(i)
		if len(error) > 0 and not mol2_check:
			self.leap_exec(mol2_check=True)
		elif len(error) > 0:
			for i in error:
				print(i[:-1])
			return -3
		else:
			return 0

	def simulation(self,arq='gpu', min_files='Minimization', heating_files='Annealing', equil_files='Equilibration',
	ph_range=[7.0,7.0]):
		'''Creation of the shell script that runs the Amber (or Ambertools) simulation packages. 

Returns -3 if there is any FATAL error in trying to build the topology.
Returns -2 if something is strange with the system charge (for explicit solvent only).
Returns -1 if in a CpHMD your chosen titratable residues are not in cpinutil.py database.
Returns 0 if chosen to stop the simulatio after prepararing the system and succesfully creating the shell file.

Parameters
----------
arq:
	'gpu' -> Means to run with pmemd.cuda;
	'sander' -> Means to run with the open-source package (Ambertools);
	'pmemd' -> Means to run only with pmemd (in case there is no gpu available).

ph_range: Initial and final ph-value for CpHMD, with self.pH_step pH-unit intervals.

min_files: Minimization files name (without .extension).

heating_files: Annealing files name (without .extension).

equil_files: Equilibration files name (without .extension).'''

		if self.pdb_string == '':
			# No mutation was done
			leap_checking = self.leap_exec()
			if leap_checking == -2:
				return -2
			elif leap_checking == -3:
				return -3

		# Cphmd variable. It doesn't matter its value if not self.cph
		res_cnt = 1 
		if self.cph:
			cmd('rm -r Cph')
			cmd('mkdir Cph')

			# The following will adjust the residues to names AMBER understands and test if any of them are in the AMBER's database.
			#  It's no walk in the park adding a residue in the database, feel free to try it yourself!
			ss = ''
			for i in range(len(self.protonation_res)):
				if self.protonation_res[i] == 'ASP':
					if ss != '':
						ss += ' '
					ss += 'AS4'
				elif self.protonation_res[i] == 'GLU':
					if ss != '':
						ss += ' '
					ss += 'GL4'
				elif self.protonation_res[i] == 'HIS':
					if ss != '':
						ss += ' '
					ss += 'HIP'
				elif self.protonation_res[i] in Titratable_Residue_Names:
					if ss != '':
						ss += ' '
					ss += self.protonation_res[i]

			if ss == '':
				print('Chosen titratable residues not in cpinutil.py database!\n')
				return -1

			if self.exp_solv:
				cmd('cpinutil.py -p %s -resnames %s -igb 2 -o system.cpin -op system_op.prmtop'%(self.prmtop,ss))
				self.prmtop = 'system_op.prmtop'
			else:
				cmd('cpinutil.py -p %s -resnames %s -igb 2 -o system.cpin'%(self.prmtop,ss))

			# TRESCNT is the number of chosen titratable residues for the simulation
			g = open('system.cpin','r')
			tg = g.readlines()
			g.close()
			for jj in tg:
				if 'TRESCNT' in jj:
					dat_1 = jj.find('TRESCNT',0,len(jj))
					dat_2 = jj.find(',',dat_1,len(jj))
					res_cnt = jj[dat_1 + 8:dat_2]
					break

		print('Preparation step completed.')
		# Creating input for minimization, annealing, equilibration and production
		if not self.mode_custom:
			self.input_min(min_name=min_files)
			self.input_heat(annealing_name=heating_files, mdsteps_factor=self.AnnealEqFactor*self.mode[self.simulation_mode][0])
			self.input_equil(step=equil_files, mdsteps_factor=self.mode[self.simulation_mode][0], trescnt=int(res_cnt))
			self.prod_cph(mdsteps_factor=self.mode[self.simulation_mode][1], ph=ph_range, trescnt=int(res_cnt))
		else:
			self.input_min(min_name=min_files)
			self.input_heat(annealing_name=heating_files)
			self.input_equil(step=equil_files, trescnt=int(res_cnt))
			self.prod_cph(ph=ph_range, trescnt=int(res_cnt))
			
		# Minimization-annealing-equilibration-prodMD/CpHMD
		cmd('rm simulation.sh')
		f = open('simulation.sh','w')
		description = 'Script - Molecular dynamics in AMBER'
		f.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		f.write('rm -r %s %s Analysis %s\n'%(min_files,heating_files,equil_files))
		f.write('mkdir %s %s Analysis %s\n'%(min_files,heating_files,equil_files))
		if self.hmr:
			self.hmr_run()
			f.write('sh hmr-in.sh\n')
		f.write('cp %s Analysis/ \n'%self.prmtop)
		# Report program will be created on our "home" directory
		self.report_program(min_name=min_files, heat_name=heating_files, eq_name=equil_files, phs=ph_range)

		f.write('sh simMin.sh\n')
		f.close()
		if arq.lower() == 'gpu':
			arq='pmemd.cuda'
		self.minimization_to_analysis(arq_type=arq, min_name=min_files, heat_name=heating_files, eq_name=equil_files, phs=ph_range)
		self.analysis(mini_name=min_files, heat_name=heating_files, eq_name=equil_files, phs=ph_range)
		cmd('ls > AfterManager.out')
		cmd('diff BeforeManager.out AfterManager.out > diff.out')
		print('\nThe Following files were created:')
		dii = open('diff.out','r')
		diit = dii.readlines()
		dii.close()
		for i in diit:
			if '>' in i and 'AfterManager' not in i:
				print(i[1:-1])
		print('')
		cmd('rm diff.out BeforeManager.out AfterManager.out')
		if not self.prep_stop:
			cmd('sh simulation.sh &')
		else:
			print('Stopping simulation. You can continue anytime with:\n $ sh simulation.sh &')
			return 0

	def report_program(self, reportFile="Simulation_Report.README", min_name='min_files', heat_name='heating_files', eq_name='equil_files', phs='ph_range'):
		f = open(self.report_file,'w')
		f.write('from os import system as cmd\n')
		f.write('from os import stat\n\n')
		f.write('def report(file="Stage.mdout", reportFile=\"%s\", inputFile = "Stage.in", flag_file = "ContinueMD.in", rst_rename= True, mdStageStep = 5*10**5):\n'%reportFile)
		f.write('\tf = open(file,"r")\n')
		f.write('\ttemp = f.readlines()\n')
		f.write('\tf.close()\n')
		f.write('\tnstep    = -1\n')
		f.write('\tresults  = False\n')
		f.write('\tall_done = False\n')
		f.write('\tnotes   = \'\'\n')
		f.write('\tcontinue_sim = True\n\n')
		## Error verification (reading mdout file)
		f.write('\tfor i in range(len(temp)):\n')
		f.write('\t\tif not results:\n')
		# Did the stage began?
		f.write('\t\t\tif "4" in temp[i] and "RESULT" in temp[i]:\n')
		f.write('\t\t\t\tresults = True\n')
		f.write('\t\telse:\n')
		# Since the stage began, "read" the last saved snapshot's step number
		f.write('\t\t\tif "NSTEP" in temp[i]:\n')
		f.write('\t\t\t\tnstep = i\n')
		# If it was displayed the timings, it means the stage finished without problems
		f.write('\t\t\tif "5" in temp[i] and "TIMING" in temp[i]:\n')
		f.write('\t\t\t\tall_done = True\n')
		f.write('\t\t\t\tbreak\n\n')
		## Error report 
		f.write('\tif not all_done:\n')
		# If there was some kind of problem
		f.write('\t\tif not results:\n')
		# Fatal error -> stop the simulation 
		f.write('\t\t\tnotes = \'%s stage didn\\\'t started due to unknown reasons, please recheck your system, the input and the mdout files!\\n\'%file[:-6]\n')
		f.write('\t\t\tcontinue_sim = False\n')
		f.write('\t\telse:\n')
		# Non-fatal error -> create input for restart
		f.write('\t\t\tline_step  = temp[nstep].split()\n')
		f.write('\t\t\tif "NSTEP=" in line_step[0] and len(line_step[0]) > 6:\n')
		f.write('\t\t\t\t#NSTEP=500\n')
		f.write('\t\t\t\tfinal_step = int(line_step[0][6:])\n')
		f.write('\t\t\telif "NSTEP=" == line_step[0]:\n')
		f.write('\t\t\t\t#NSTEP= 500\n')
		f.write('\t\t\t\tfinal_step = int(line_step[1])\n')
		f.write('\t\t\telif "NSTEP" == line_step[0] and "=" in line_step[1] and len(line_step[1]) > 1:\n')
		f.write('\t\t\t\t#NSTEP =500\n')
		f.write('\t\t\t\tfinal_step = int(line_step[1][1:])\n')
		f.write('\t\t\telif "NSTEP" == line_step[0] and "=" == line_step[1]:\n')
		f.write('\t\t\t\t#NSTEP = 500\n')
		f.write('\t\t\t\tfinal_step = int(line_step[2])\n')
		f.write('\t\t\tcalc = final_step/mdStageStep\n')
		# Verification: if calc<1, restart the whole stage
		f.write('\t\t\tif calc < 1:\n')
		f.write('\t\t\t\tnotes = \'%s stage stopped due to unknown reasons.\\n\'%file[:-6]\n')
		# Creating input for restart file: inputFile[:-3]+"_rst.in"
		f.write('\t\t\t\tif \'%s\' not in inputFile:\n'%min_name) #minimization is not that long so just redo the whole stage
		f.write('\t\t\t\t\tnotes += \'Generating restart files!\\n\'\n')
		# Reading old input
		f.write('\t\t\t\t\tg_old = open(inputFile,"r")\n')
		f.write('\t\t\t\t\tg_read = g_old.readlines()\n') 
		f.write('\t\t\t\t\tg_old.close()\n')
		# Creating restart input
		f.write('\t\t\t\t\tg = open(inputFile[:-3]+"_rst.in","w")\n')
		f.write('\t\t\t\t\tfor j in range(len(g_read)):\n')
		# If MD step flag was found, modify it to the new number 
		f.write('\t\t\t\t\t\tif "nstlim=" in g_read[j]:\n')
		# Between bit1 and bit2 was the old number for MD steps
		f.write('\t\t\t\t\t\t\tbit1 = g_read[j].find("nstlim=",0,len(g_read[j]))+7\n')
		f.write('\t\t\t\t\t\t\tbit2 = g_read[j].find(",",bit1,len(g_read[j]))\n')
		# Adjusting MD steps
		f.write('\t\t\t\t\t\t\tj_temp = g_read[j][:bit1]+str(mdStageStep - final_step)+g_read[j][bit2:]\n')
		f.write('\t\t\t\t\t\t\tg_read[j] = j_temp\n')
		# Copy the current line from old input
		f.write('\t\t\t\t\t\tg.write(g_read[j])\n')
		f.write('\t\t\t\t\tg.close()\n')
		# Renaming necessary files - It doesn't matter that this py program may not be on the same directory as the files
		f.write('\t\t\t\t\tif rst_rename:\n')
		f.write('\t\t\t\t\t\tcmd("mv %s %s"%(file, file[:-6]+"_rst.mdout"))\n')
		f.write('\t\t\t\t\t\tcmd("mv %s %s"%(file[:-6]+".rst7", file[:-6]+"_rst.rst7"))\n')
		f.write('\t\t\t\t\t\tif \'%s\' not in inputFile:\n'%min_name)
		f.write('\t\t\t\t\t\t\tcmd("mv %s %s"%(file[:-6]+".nc", file[:-6]+"_rst.nc"))\n')
		f.write('\t\t\t\t\t\t\tcmd("mv %s %s"%(file[:-6]+".cpout", file[:-6]+"_rst.cpout"))\n')
		# No error -> reporting that everything went well		
		f.write('\telse:\n')
		f.write('\t\tnotes = \'%s stage finished without problems!\\n\'%file[:-6]\n')
		# For the report file we must know if the file alredy exists (from a previous stage)
		f.write('\ttry:\n')
		f.write('\t\tinfo_sh = stat(reportFile)\n')
		f.write('\t\tf = open(reportFile,"r+")\n')
		# If the file exists, seek the last bit to make sure we are not losing information!
		f.write('\t\tf.seek(info_sh.st_size)\n')
		f.write('\texcept FileNotFoundError:\n')
		f.write('\t\tf = open(reportFile,"w")\n')
		f.write('\tf.write(notes)\n')
		f.write('\tf.close()\n')
		# Creating a flag file to know if it's ok to continue to the next MD stage
		f.write('\tf = open(flag_file,"w")\n')
		f.write('\tif continue_sim:\n')
		f.write('\t\tif not all_done:\n')
		f.write('\t\t\tf.write("restart")\n')
		f.write('\t\telse:\n')
		f.write('\t\t\tf.write("yes")\n')
		f.write('\telse:\n')
		f.write('\t\tf.write("no")\n')
		f.write('\tf.close()\n')
		## Main
		f.write('if __name__ == "__main__":\n')
		f.write('\timport sys\n')
		f.write('\targ = sys.argv\n')
		# 'bit' looks for the path of this report code 
		f.write('\tbit = arg[0].find("%s")\n'%self.report_file)
		f.write('\treport_name = arg[0][:bit]+\'%s\'\n'%reportFile)
		f.write('\tflag_name   = "ContinueMD.in"\n')
		f.write('\trname_flag  = True\n') 
		f.write('\tif len(arg) == 4:\n\t\tflag_name = arg[3]\n\t\trname_flag = False\n')
		##                  arg[0]                arg[1]     arg[2]
		## Use: $ python3 ASMhome/report_temp.py Stage.mdout ASMhome/Stage.in
		'''mdStageStep => Min: fixed maxcyc=10000# Anneal-Equil-Prod:mdsteps_factor*self.steps;
		self.input_heat  -> mdsteps_factor=self.AnnealEqFactor*self.mode[self.simulation_mode][0];
		self.input_equil -> mdsteps_factor=self.mode[self.simulation_mode][0];
		self.prod_cph    -> mdsteps_factor=self.mode[self.simulation_mode][1].'''
		if not self.mode_custom:
			f.write('\tif \'%s\' in  arg[1]:\n\t\tmdstep = 10**4\n'%min_name)
			f.write('\telif \'%s\' in  arg[1]:\n\t\tmdstep = %d\n'%(heat_name, self.AnnealEqFactor*self.mode[self.simulation_mode][0]*self.steps))
			f.write('\telif \'%s\' in  arg[1]:\n\t\tmdstep = %d\n'%(eq_name, self.mode[self.simulation_mode][0]*self.steps))
			f.write('\telif \'Production\' in  arg[1] or \'CpHMD_prod\' in  arg[1]:\n\t\tmdstep = %d\n'%(self.mode[self.simulation_mode][1]*self.steps))	
		else:
			f.write('\tif \'%s\' in  arg[1]:\n\t\tmdstep = %d\n'%(min_name, self.mode_custom_min))
			f.write('\telif \'%s\' in  arg[1]:\n\t\tmdstep = %d\n'%(heat_name, self.mode_custom_a))
			f.write('\telif \'%s\' in  arg[1]:\n\t\tmdstep = %d\n'%(eq_name, self.mode_custom_e))
			f.write('\telif \'Production\' in  arg[1] or \'CpHMD_prod\' in  arg[1]:\n\t\tmdstep = %d\n'%(self.mode_custom_p))	
		f.write('\treport(file=arg[1], reportFile=report_name, inputFile = arg[2], rst_rename=rname_flag, mdStageStep = mdstep, flag_file= flag_name)\n')
		f.close()
		
	def analysis(self, inp = 'rmsdf.cpptraj', mini_name='minimization', heat_name='annealing',
	eq_name='equilibration', phs=[0.0,7.0]):
		'''Trajectory analysis.

Parameters
----------

inp: Cpptraj input file's name.

mini_name: Minimization's files name.

heat_name: Annealing's files name.

eq_name: Equilibration's files name.

phs: Ph range used to titrate.'''

		ph_i = phs[0]
		pH   = []
		while ph_i <= phs[1]:
			pH.append(ph_i)
			ph_i += self.pH_step

		f = open('simAnalysis.sh','w')
		description = 'Script - Minimization stage'
		f.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		f.write('cd Analysis\n')
		
		if not self.cph:
			for pp in pH:
				self.cpptraj_in(input_name=inp, rms_file=self.rms, rmsf_file=self.rmsf, min_name=mini_name, annealing_name=heat_name, equil_name=eq_name, phr=phs)
				# We need to look for the rmsd input file in the father directory because 'f' is writting in the simulation.sh file which will run later
				f.write('cpptraj -i ../MD_%.2f_%s &> cpptraj.log\n'%(pp,inp))
				f.write('process_mdout.perl %s.mdout %s.mdout Production.mdout\n'%(heat_name, eq_name))
		else:
			f.write('echo \"Offset: is the difference between the predicted pka and the system pH\" > ph_calcpka_populations.info\n')
			f.write('echo \"Pred: is the predicted pka\" >> ph_calcpka_populations.info\n')
			f.write('echo \"Frac Prot: is the fraction of time the residue spends protonated\" >> ph_calcpka_populations.info\n')
			f.write('echo \"Transitions:gives the number of accpeted protonations state transitions\" >> ph_calcpka_populations.info\n')
			f.write('echo \"Average total molecular protonation: is the sum of the fractional protonations\\n\" >> ph_calcpka_populations.info\n')
			f.write('echo \"populations.dat: prints the fraction of snapshots that the system spent in each state for each residue, where 1.0 means 100% of the time.\\n  OBS:The parentheses (*), indicates the number protons present in that state.\" >> ph_calcpka_populations.info\n')
			for pp in pH:
				if not self.exp_solv:
					# There is a problem in AMBERtools20 with the following command, but it works just fine on AMBERtools18
					f.write('cphstats -i ../system_%s.cpin CpHMD_prod_%.2f.cpout -o pH%.2f_calcpka.dat --population pH%.2f_populations.dat\n'%(eq_name, pp, pp, pp)) # Look in https://github.com/Amber-MD/cpptraj/pull/816/files
				else:
					f.write('cphstats -i ../system.cpin CpHMD_prod_%.2f.cpout -o pH%.2f_calcpka.dat --population pH%.2f_populations.dat\n'%(pp, pp, pp))
				self.cpptraj_in(input_name='rmsdf_CpHMD_%.2f.cpptraj'%pp, rms_file=self.rms, rmsf_file=self.rmsf, min_name=mini_name, annealing_name=heat_name, equil_name=eq_name, phr=phs)
				f.write('cpptraj -i ../rmsdf_CpHMD_%.2f.cpptraj &> cpptraj_%.2f.log\n'%(pp,pp))
		f.close()

	def hmr_run(self, inp1='parmedhmr.in', prmtop_new='systemHMR.prmtop'):
		'''Shell file which is responsible for the HMR transform in the topology file.

Parameters
----------

inp1: Name for the Parmed input.

prmtop_new: Name chosen for the new topoly file.'''

		prmtop = cp(self.prmtop)
		self.hmr_transform_file(prmtop_new=prmtop_new, input_name=inp1)
		cmd('rm hmr-in.sh')
		f = open('hmr-in.sh','w')
		description = 'Temporary component of the Amber_run class'
		f.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		# Old topology on 'prmtop'. This cleans the directory of a possible "systemHMR.prmtop" file
		f.write('rm %s\n'%self.prmtop)
		f.write('parmed %s < %s'%(prmtop,inp1))
		f.close()

	def minimization_to_analysis(self, arq_type='pmemd', min_name='minimization', heat_name='annealing',
	eq_name='equilibration', phs=[7.0,7.0]):
		'''Amber run without GPU.

Parameters
----------

arq_type: Amber binary chosen to run, between 'pmemd', 'sander' or 'pmemd.cuda'

min_name: Minimization files name (without .extension).

heat_name: Annealing files name (without .extension).

eq_name: Equilibration files name (without .extension).

phs: Initial and final ph-value for CpHMD, with self.pH_step pH-unit intervals.'''

		# MPI choice
		if self.mpi_use and arq_type.lower() != 'pmemd.cuda':
			process = '%s -n %d %s.MPI'%(self.mpi, self.cores, arq_type)
		else:
			process = arq_type.lower()

		# Verification file template
		sh_verif = '''#!/bin/bash
read simflag
if test $simflag = "yes"
then
	%s
	cd ..
	%s
elif test $simflag = "restart"
then
	%s
	%s
	while read RLINE
	do
		flag=$RLINE
	done < Continue_rst2.in
	if test $flag = "yes"
	then
		%s
		cd ..
		%s
	else
		echo "Stopping run..."
	fi
else
	echo "Stopping run..."
fi'''
		keys = []
		ph_i = phs[0]
		pH   = []
		while ph_i <= phs[1]:
			pH.append(ph_i)
			ph_i += self.pH_step
		cphmd_names          = ['CpHMD_prod_%.2f'%p for p in pH]
		productionMD_names   = ['Production_%.2f'%p for p in pH]
		cphmd_SHnames        = []
		productionMD_SHnames = []
		for pp in range(len(cphmd_names)):
			if pp == 0:
				cphmd_SHnames.append( 'sh sim%s.sh'%cphmd_names[pp] )
			else:
				cphmd_SHnames.append( 'cd .. && sh sim%s.sh'%cphmd_names[pp] )
		for pp in range(len(productionMD_names)):
			if pp == 0:
				productionMD_SHnames.append( 'sh sim%s.sh'%productionMD_names[pp] )
			else:
				productionMD_SHnames.append( 'cd .. && sh sim%s.sh'%productionMD_names[pp] )

		comp_report = 'python3 ../%s %s_rst.mdout %s_rst.in Continue_rst2.in'
		if self.cph:
			comp_report = 'python3 ../../%s %s_rst.mdout %s_rst.in Continue_rst2.in'
			stages    = [min_name, heat_name, eq_name]
			stages.extend( cphmd_names )
			stages.append('')
			stages_sh = ['sh simMin.sh', 'sh simAnneal.sh', 'sh simEquil.sh']
			stages_sh.extend( cphmd_SHnames )
			stages_sh.append('sh simAnalysis.sh')
			#stages_sh.append('') # '' := After production do nothing
		else:
			stages    = [min_name, heat_name, eq_name]
			stages.extend( productionMD_names )
			stages.append('')
			stages_sh = ['sh simMin.sh', 'sh simAnneal.sh', 'sh simEquil.sh']
			stages_sh.extend( productionMD_SHnames )
			stages_sh.append('sh simAnalysis.sh')
			#stages_sh.append('') # '' := After production do nothing

		for txt_temp in range(1,len(stages)):
			txt        = stages[txt_temp-1]
			next_stage = stages[txt_temp]
			next_sh    = stages_sh[txt_temp]
			if self.cph:
				if txt == min_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -ref ../systemLeap.rst7 -cpin ../system.cpin'%(process, txt, txt, self.prmtop, txt, txt)
					copyingTo     = 'cp %s.rst7 ../%s && cp %s.rst7 ../Analysis'%(txt, next_stage, txt)
					copyingRstTo  = copyingTo
				elif txt == heat_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc -ref %s_rst.rst7 -cpin ../system.cpin'%(process, txt, txt, self.prmtop, txt, txt, txt, txt)
					copyingTo     = 'cp %s.rst7 ../%s && cp %s.nc ../Analysis  && cp %s.mdout ../Analysis'%(txt, next_stage, txt, txt)
					copyingRstTo  = 'cp %s.rst7 ../%s && cp %s.nc ../Analysis  && cp %s.mdout ../Analysis && cp %s_rst.mdout ../Analysis && cp %s_rst.nc ../Analysis'%(txt, next_stage, txt, txt, txt, txt)
				elif txt == eq_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc -cpin system_%s.cpin -cpout system_%s.cpout -cprestrt system_%s_rstP2.cpin'%(process, txt, txt, self.prmtop, txt, txt, txt, txt, txt, txt)
					copyingTo     = 'cp %s.rst7 ../Cph && cp %s.nc ../Analysis && cp %s.mdout ../Analysis && cp %s.rst7 ../Analysis && cp system_%s.cpin ../'%(txt, txt, txt, txt, txt)
					copyingRstTo  = 'cp %s.rst7 ../Cph && cp %s.nc ../Analysis && cp %s.mdout ../Analysis && cp %s.rst7 ../Analysis && cp %s_rst.mdout ../Analysis && cp %s_rst.nc ../Analysis && cp system_%s.cpin ../'%(txt, txt, txt, txt, txt, txt, txt)
				elif 'CpHMD' in txt:
					redoing_stage = '%s -O -i ../../%s_rst.in -o %s.mdout -p ../../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc -cpin %s.cpin -cpout %s.cpout -cprestrt %s_rstP2.cpin'%(process, txt, txt, self.prmtop, txt, txt, txt, txt, txt, txt)
					copyingTo     = 'cp %s.cpout ../../Analysis && cp %s.nc ../../Analysis && cp %s.mdout ../../Analysis'%(txt,txt,txt)
					copyingRstTo  = 'cp %s.cpout ../../Analysis && cp %s.nc ../../Analysis && cp %s.mdout ../../Analysis && cp %s_rst.mdout ../../Analysis && cp %s_rst.nc ../../Analysis'%(txt,txt,txt,txt,txt)
			else:
				if txt == min_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7'%(process, txt, txt, self.prmtop, txt, txt)
					copyingTo     = 'cp %s.rst7 ../%s && cp %s.rst7 ../Analysis'%(txt, next_stage, txt)
					copyingRstTo  = copyingTo
				elif txt == heat_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc'%(process, txt, txt, self.prmtop, txt, txt, txt)
					copyingTo     = 'cp %s.rst7 ../%s && cp %s.nc ../Analysis && cp %s.mdout ../Analysis'%(txt, next_stage, txt, txt)
					copyingRstTo  = 'cp %s.rst7 ../%s && cp %s.nc ../Analysis && cp %s.mdout ../Analysis && cp %s_rst.mdout ../Analysis && cp %s_rst.nc ../Analysis'%(txt, next_stage, txt, txt, txt, txt)
				elif txt == eq_name:
					redoing_stage = '%s -O -i ../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc'%(process, txt, txt, self.prmtop, txt, txt, txt)
					copyingTo     = 'cp %s.nc ../Analysis && cp %s.mdout ../Analysis && cp %s.rst7 ../Analysis'%(txt, txt, txt)
					copyingRstTo  = 'cp %s.nc ../Analysis && cp %s.mdout ../Analysis && cp %s.rst7 ../Analysis && cp %s_rst.mdout ../Analysis && cp %s_rst.nc ../Analysis'%(txt, txt, txt, txt, txt)
				elif 'Production' in txt:
					redoing_stage = '%s -O -i ../../%s_rst.in -o %s.mdout -p ../%s -c %s_rst.rst7 -r %s.rst7 -x %s.nc'%(process, txt, txt, self.prmtop, txt, txt, txt)
					copyingTo     = 'cp %s.nc ../../Analysis && cp %s.mdout ../../Analysis'%(txt, txt)
					copyingRstTo  = 'cp %s.nc ../../Analysis && cp %s.mdout ../../Analysis && cp %s_rst.mdout ../../Analysis && cp %s_rst.nc ../../Analysis'%(txt, txt, txt, txt)
			
			keys.append( ( copyingTo, next_sh, redoing_stage, comp_report%(self.report_file, txt, txt), copyingRstTo, next_sh) )

		# Stage Minimization
		m = open('simMin.sh','w')
		description = 'Script - Minimization stage'
		m.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		m.write('export CUDA_VISIBLE_DEVICES=%s\n'%self.gpunumber)
		m.write('cd Minimization\n')
		if self.cph:
			m.write('echo \'Running %s...\\n\'\n'%min_name)
			m.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c ../systemLeap.rst7 -r %s.rst7 -ref ../systemLeap.rst7 -cpin ../system.cpin\n'%(process, min_name, min_name, self.prmtop, min_name))
		else:
			m.write('echo \'Running %s...\\n\'\n'%min_name)
			m.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c ../systemLeap.rst7 -r %s.rst7\n'%(process, min_name, min_name, self.prmtop, min_name))
		m.write('python3 ../%s %s.mdout ../%s.in\n'%(self.report_file, min_name, min_name))
		# ContinueMD.in := "restart"; "yes"; "no".
		# Verification_sh files are on the same directory of simStage_sh, 
		# but the verification_sh are being executed inside each stage's directory.
		m.write('sh ../verify_min.sh < ContinueMD.in\n')
		m.close()
		m = open('verify_min.sh','w')
		m.write(sh_verif%keys[0])
		m.close()
		
		# Stage Annealing
		a = open('simAnneal.sh','w')
		description = 'Script - Annealing stage'
		a.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		#a.write('export CUDA_VISIBLE_DEVICES=%s\n'%self.gpunumber)
		a.write('cd Annealing\n')
		if self.cph:
			a.write('echo \'Running %s...\\n\'\n'%heat_name)
			a.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c %s.rst7 -r %s.rst7 -x %s.nc -ref %s.rst7 -cpin ../system.cpin\n'%(process, heat_name, heat_name, self.prmtop, min_name, heat_name, heat_name, min_name))
		else:
			a.write('echo \'Running %s...\\n\'\n'%heat_name)
			a.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c %s.rst7 -r %s.rst7 -x %s.nc\n'%(process, heat_name,heat_name,self.prmtop,min_name,heat_name,heat_name))
		a.write('python3 ../%s %s.mdout ../%s.in\n'%(self.report_file, heat_name, heat_name))
		a.write('sh ../verify_anneal.sh < ContinueMD.in\n')
		a.close()
		a = open('verify_anneal.sh','w')
		a.write(sh_verif%keys[1])
		a.close()

		# Stage Equilibration
		e = open('simEquil.sh','w')
		description = 'Script - Equilibration stage'
		e.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
		#e.write('export CUDA_VISIBLE_DEVICES=%s\n'%self.gpunumber)
		e.write('cd %s\n'%eq_name)
		if self.cph:
			e.write('echo \'Running %s...\\n\'\n'%eq_name)
			e.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c %s.rst7 -r %s.rst7 -x %s.nc -cpin ../system.cpin -cpout system_%s.cpout -cprestrt system_%s.cpin\n'%(process, eq_name, eq_name, self.prmtop, heat_name, eq_name, eq_name, eq_name, eq_name))
			e.write('python3 ../%s %s.mdout ../%s.in\n'%(self.report_file, eq_name, eq_name))
			e.write('sh ../verify_equil.sh < ContinueMD.in\n')
			e.close()
			e = open('verify_equil.sh','w')
			e.write(sh_verif%keys[2])
			e.close()

			# Stage Production for constant pH
			cphmd_name_counter = 0
			for Name in cphmd_names:
				cph = open('sim%s.sh'%Name,'w')			
				description = 'Script - Production stage for CpHMD'
				cph.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
				#cph.write('export CUDA_VISIBLE_DEVICES=%s\n'%self.gpunumber)
				cph.write('cd Cph\n')
				cph.write('rm -r %s\nmkdir %s\ncd %s\n'%(Name, Name, Name))
				cph.write('echo \'Running %s...\\n\'\n'%Name)
				if not self.exp_solv:
					cph.write('%s -O -i ../../%s.in -o %s.mdout -p ../../%s -c ../%s.rst7 -r %s.rst7 -x %s.nc -cpin ../../system_%s.cpin -cpout %s.cpout -cprestrt %s.cpin\n'%(process, Name, Name, self.prmtop, eq_name, Name, Name, eq_name, Name, Name))
				else:
					cph.write('%s -O -i ../../%s.in -o %s.mdout -p ../../%s -c ../%s.rst7 -r %s.rst7 -x %s.nc -cpin ../../system.cpin -cpout %s.cpout -cprestrt %s.cpin\n'%(process, Name, Name, self.prmtop, eq_name, Name, Name, Name, Name))
				cph.write('python3 ../../%s %s.mdout ../../%s.in\n'%(self.report_file, Name, Name))
				cph.write('sh ../../verify_%s.sh < ContinueMD.in\n'%Name)
				cphVerif = open('verify_%s.sh'%Name,'w')
				cphVerif.write(sh_verif%(keys[3+cphmd_name_counter]))
				cphmd_name_counter += 1
				cphVerif.close()
				if cphmd_name_counter == len(cphmd_names)+1:
					cph.write('cd ..\n')
				cph.close()
		else:
			e.write('echo \'Running %s...\\n\'\n'%eq_name)
			e.write('%s -O -i ../%s.in -o %s.mdout -p ../%s -c %s.rst7 -r %s.rst7 -x %s.nc\n'%(process, eq_name, eq_name, self.prmtop, heat_name, eq_name, eq_name))
			e.write('python3 ../%s %s.mdout ../%s.in\n'%(self.report_file, eq_name, eq_name))
			e.write('sh ../verify_equil.sh < ContinueMD.in\n')
			e.close()
			e = open('verify_equil.sh','w')
			e.write(sh_verif%keys[2])
			e.close()

			# Stage Production - we are being redundant here for an easier management of the dynamics stages and the error report.
			name_counter = 0
			for Name in productionMD_names:
				prod = open('sim%s.sh'%Name,'w')
				description = 'Script - Production stage'
				prod.write('#!/bin/bash\n#\n#%s\n#\n#Autor:%s\n#Email:%s\n%s\n'%(description,self.autor,self.email,self.linebreak))
				#prod.write('export CUDA_VISIBLE_DEVICES=%s\n'%self.gpunumber)
				prod.write('cd %s\n'%eq_name)
				prod.write('rm -r %s\nmkdir %s\ncd %s\n'%(Name, Name, Name))
				prod.write('echo \'Running %s...\\n\'\n'%Name)
				prod.write('%s -O -i ../../%s.in -o %s.mdout -p ../../%s -c ../%s.rst7 -r %s.rst7 -x %s.nc\n'%(process,Name,Name, self.prmtop, eq_name,Name,Name))
				prod.write('ambpdb -p ../../%s -c %s.nc > %s.pdb\n'%(self.prmtop,Name,Name))
				prod.write('python3 ../../%s %s.mdout ../../%s.in\n'%(self.report_file,Name,Name))
				prod.write('sh ../../verify_%s.sh < ContinueMD.in\n'%Name)
				prodVerif = open('verify_%s.sh'%Name,'w')
				prodVerif.write(sh_verif%keys[3+name_counter])
				name_counter += 1
				prodVerif.close()
				if name_counter == len(productionMD_names)+1:
					prod.write('cd ..\n')
				prod.close()

class Amber_mutation(Amber_run):
	'''Shell script's manager extension for residues mutations'''
	aminoacid_atmcount = {'ARG':12, 'HIS':11, 'LYS':10, 'ASP':9, 'GLU':10, 'SER':7, 'THR': 8, 'ASN':9, 'GLN':10, 'CYS':7, 'GLY':5, 'PRO':8, 'ALA':6, 'VAL':8, 'ILE':9, 'LEU':9, 'MET':9, 'PHE':12, 'TYR':13, 'TRP':15}
	mutation_links={'ARG':['ALA','GLY'], 'HIS':['ALA','GLY'], 'LYS':['ALA','GLY'], 'ASP':['ALA','GLU','GLY'], 'GLU':['ALA','ASP','GLY'], 'SER':['ALA','THR','GLY'], 'THR': ['ALA','SER','GLY'], 'ASN':['ALA','GLN','GLY'], 'GLN':['ALA','ASN','GLY'], 'CYS':['ALA','MET','GLY'], 'GLY':['ALA'], 'PRO':['ALA','GLY'], 'ALA':['VAL','GLY'], 'VAL':['ALA','ILE','GLY'], 'ILE':['ALA','VAL','LEU','GLY'], 'LEU':['ALA','ILE','GLY'], 'MET':['ALA','CYS','LEU','GLY'], 'PHE':['ALA','TYR','GLY'], 'TYR':['ALA','PHE','GLY'], 'TRP':['ALA']}
	improbable_mutation={
'ARG': ['HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'], 
'HIS': ['ARG', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'], 
'LYS': ['ARG', 'HIS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'], 
'ASP': ['ARG', 'HIS', 'LYS', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'GLU': ['ARG', 'HIS', 'LYS', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'SER': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'THR': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'ASN': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'GLN': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'CYS': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'PRO', 'VAL', 'ILE', 'LEU', 'PHE', 'TYR', 'TRP'],
'GLY': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'PRO': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'ALA': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'VAL': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'],
'ILE': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'MET', 'PHE', 'TYR', 'TRP'],
'LEU': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'MET', 'PHE', 'TYR', 'TRP'],
'MET': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'PRO', 'VAL', 'ILE', 'PHE', 'TYR', 'TRP'],
'PHE': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'TRP'],
'TYR': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'TRP'],
'TRP': ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR']}
	restypeatoms = {
'ALA': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB1', 'HC', 'HB2', 'HB3', 'C', 'O'],
'ARG': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'N2', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O'],
'ASH': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'C', 'OD1', 'O', 'OD2', 'OH', 'HD2', 'HO'],
'ASN': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'C', 'OD1', 'O', 'ND2', 'HD21', 'HD22'],
'ASP': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'C', 'OD1', 'O2', 'OD2', 'O'],
'CYM': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB3', 'HB2', 'SG', 'SH', 'C', 'O'],
'CYS': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HB3', 'SG', 'SH', 'HG', 'HS', 'C', 'O'],
'CYX': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HB3', 'SG', 'S', 'C', 'O'],
'GLH': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'C', 'OE1', 'O', 'OE2', 'OH', 'HE2', 'HO'],
'GLN': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'C', 'OE1', 'O', 'NE2', 'HE21', 'HE22'],
'GLU': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'C', 'OE1', 'O2', 'OE2', 'O'],
'GLY': ['N', 'H', 'CA', 'CX', 'HA2', 'H1', 'HA3', 'C', 'O'],
'HID': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'CC', 'ND1', 'NA', 'HD1', 'CE1', 'CR', 'HE1', 'H5', 'NE2', 'NB', 'CD2', 'CV', 'HD2', 'H4', 'C', 'O'],
'HIE': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'CC', 'ND1', 'NB', 'CE1', 'CR', 'HE1', 'H5', 'NE2', 'NA', 'HE2', 'CD2', 'CW', 'HD2', 'H4', 'C', 'O'],
'HIP': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'CC', 'ND1', 'NA', 'HD1', 'CE1', 'CR', 'HE1', 'H5', 'NE2', 'HE2', 'CD2', 'CW', 'HD2', 'H4', 'C', 'O'],
'HYP': ['N', 'CD', 'CT', 'HD2', 'H1', 'HD3', 'CG', 'HG2', 'OD1', 'OH', 'HD1', 'HO', 'CB', 'HB2', 'HC', 'HB3', 'CA', 'CX', 'HA', 'C', 'O'],
'ILE': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB', 'HC', 'CG2', 'HG21', 'HG22', 'HG23', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'C', 'O'],
'LEU': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG', 'CD1', 'HD11', 'HD12', 'HD13', 'CD2', 'HD21', 'HD22', 'HD23', 'C', 'O'],
'LYN': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HP', 'HE3', 'NZ', 'N3', 'HZ2', 'HZ3', 'C', 'O'],
'LYS': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HP', 'HE3', 'NZ', 'N3', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'],
'MET': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'HG2', 'HG3', 'SD', 'S', 'CE', 'HE1', 'HE2', 'HE3', 'C', 'O'],
'PHE': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2', 'C', 'O'],
'PRO': ['N', 'CD', 'CT', 'HD2', 'H1', 'HD3', 'CG', 'HG2', 'HC', 'HG3', 'CB', 'HB2', 'HB3', 'CA', 'CX', 'HA', 'C', 'O'],
'SER': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HB3', 'OG', 'OH', 'HG', 'HO', 'C', 'O'],
'THR': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB', 'CG2', 'HG21', 'HC', 'HG22', 'HG23', 'OG1', 'OH', 'HG1', 'HO', 'C', 'O'],
'TRP': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'C*', 'CD1', 'CW', 'HD1', 'H4', 'NE1', 'NA', 'HE1', 'CE2', 'CN', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'HZ3', 'CE3', 'HE3', 'CD2', 'C', 'O'],
'TYR': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'C', 'OH', 'HH', 'HO', 'CE2', 'HE2', 'CD2', 'HD2', 'O'],
'VAL': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB', 'HC', 'CG1', 'HG11', 'HG12', 'HG13', 'CG2', 'HG21', 'HG22', 'HG23', 'C', 'O']}
	trash = {}
	res_atm = {}
	mutations_restricted = {}
	mutation_temp_prob = {}
	mutation_temp_improb = {}

	def adjusting_firstep(self, chosen_res = {'SER':[160],'HIS':[237],'ASP':[206]}):
		'''Returns a list of coords. of the atoms from the chosen_res variable. 
		Ex:{'ARG': [ [34, [coords1,coords2,...]], [53, [coords1,coords2,...]], ... ], ... } '''

		f = open(self.pdb,'r')
		temp = f.readlines()
		f.close()
		found_res=False
		res_newid={}
		for i in temp:
			data = i.split()
			if 'ATOM' in data[0] or 'ANISOU' in data[0] or 'HETATM' in data[0]:
				if not found_res:
					if len(data[3]) == 1 and data[2][-3:] in chosen_res and int(data[4]) in chosen_res[ data[2][-3:] ]:
						# ATOM    356  NH1AARG A  53 ..coords..
						found_res = True
						coord_pos = i.find(data[5],0,len(i))
						res_newid[data[2][-3:]]=[[ int(data[4]), [i[coord_pos:]] ]]
					elif len(data[3]) == 3 and data[3] in chosen_res and int(data[5]) in chosen_res[ data[3] ]:
						# ANISOU   83  NH2 ARG A  34 ..coords..
						found_res = True
						coord_pos = i.find(data[6],0,len(i))
						res_newid[data[3]]=[[int(data[5]), [ i[coord_pos:] ]]] # res_newid := [ [34, [coords1]],...]
				else:
					if len(data[3]) == 1:
						# ATOM    356  NH1AARG A  53 ..coords....
						if data[2][-3:] in res_newid:
							# If we are complementing a residue already initiated above
							for res_p in res_newid[data[2][-3:]]:
								if res_p[0] == int(data[4]):
									# res_p := [53, [coords1]]
									res_p[1].append( i[i.find(data[5],0,len(i)):] )
									# res_p := [53, [coords1,coords2]] 
									break
						elif data[2][-3:] in chosen_res and int(data[4]) in chosen_res[ data[2][-3:] ]:
							# The current line is a different aminoacid, but still one from the chosen_res list
							coord_pos = i.find(data[5],0,len(i))
							res_newid[data[2][-3:]]=[[ int(data[4]), [i[coord_pos:]] ]]
						else:
							# The current line is a different aminoacid and not one from the chosen_res list
							found_res = False
					elif len(data[3]) == 3:
						# ANISOU   83  NH2 ARG A  34 #coords....
						if data[3] in res_newid:
							# If we are complementing a residue already initiated above
							for res_p in res_newid[data[3]]:
								if res_p[0] == int(data[5]):
									# res_p := [53, [coords1]]
									res_p[1].append( i[i.find(data[6],0,len(i)):] )
									# res_p := [53, [coords1,coords2]] 
									break
						elif data[3] in chosen_res and int(data[5]) in chosen_res[ data[3] ]:
							# The current line is a different aminoacid, but still one from the active site
							coord_pos = i.find(data[6],0,len(i))
							res_newid[data[3]]=[[int(data[5]), [ i[coord_pos:] ]]]
						else:
							# The current line is a different aminoacid and not one from the active site
							found_res = False

		# res_newid := {'ARG': [ [34, ['coords1','coords2',...]], [53, ['coords1','coords2',...]], ... ], ... }
		return res_newid

	def adjusting_finalstep(self, coords={'SER': [ [160, ['-17.613  -7.764   7.103  1.00  4.45           N  \n']] ]}, id_only = False):
		'''Returns the new id of the residues in the 'coords' list.
		Returns -1 if the coords for the aminoacid are broken for some reason.

Parameters
----------

coords: List of coords. from self.adjusting_firstep.

id_only: If True, self.chosen_mutation won't change.
'''
		
		f = open(self.pdb,'r')
		temp = f.readlines()
		f.close()
		new_res_id = {}
		new_id = ''
		mut_temp = {}
		for i in coords:
			# i := residue type, eg. SER
			found = False
			for j in coords[i]:
				# j := residue id and coords, eg. [160, ['-17.613  -7.764   7.103  1.00  4.45           N  \n',...]]
				sure_flag = 0
				next_flag = False
				kpointer = 0
				for x in j[1]:
					# x:= '-17.613  -7.764   7.103  1.00  4.45           N  \n'
					for k in range(len(temp)):
						if next_flag:
							k = kpointer

						if x in temp[k] and i in temp[k]:
							coord_pos = temp[k].find(x,0,len(temp[k]))
							data      = temp[k][:coord_pos].split()
							sure_flag += 1
							if sure_flag == 2 and new_id == data[len(data)-1]:
								found = True
								break
							elif sure_flag >=2:
								print('\nA problem happened, please look for res id %s on the system-amber.pdb!\n'%data[len(data)-1])
								return -1
							new_id    = data[len(data)-1]
							next_flag = True
							kpointer  = k + 1
							break
					if found:
						if self.not_rand_mut:
							# self.chosen_mutation := {('HIS', '237'): 'ALA'} # i := 'HIP'
							i_temp = i
							if i_temp == 'AS4':
								i_temp = 'ASP'
							elif i_temp == 'HIP':
								i_temp = 'HIS'
							elif i_temp == 'GL4':
								i_temp = 'GLU'
							mut_temp[(i,str(new_id))] = self.chosen_mutation[(i_temp, str(j[0]))]

						if i not in new_res_id:
							new_res_id[i] = [int(new_id)]
						else:
							new_res_id[i].append( int(new_id) )
						break

		if self.not_rand_mut and not id_only:
			# mut_temp := {('SER', 'newid'): 'MET', ('SER', 'O1newid'): 'MET'} 
			self.chosen_mutation = mut_temp

		# new_res_id := {'SER':[newid, ...], ...}
		for i in new_res_id:
			new_res_id[i].sort()

		return new_res_id

	def res_detail(self, res_ids = {'SER':[160]}):
		'''Searches for all atoms in the res_ids and creates an attribute for "mutation possibility". Returns the tuple: (how many aminoacids are in the res_ids, how many atoms are in the residues from the list res_ids). '''

		f = open(self.pdb,'r')
		self.pdb_string = f.readlines()
		f.close()

		# mut_temp holds the lines of self.pdb_string related to the list res_ids
		mut_temp = []
		for i in self.pdb_string:
			line = i.split()
			bool1 = 'ATOM' in line[0] or 'ANISOU' in line[0] or 'HETATM' in line[0]
			bool2 = len(line) > 3
			if bool1 and bool2:
				# Now we are looking only at the atoms coords info
				numbs = re.findall(r'[0-9]+',i) # list of numbers in string format
				debug = '' # res id

				# Correcting the numbs problem that occurs when we have:
				if 'ATOM' in line[0] or 'HETATM' in line[0]:
					if len(numbs) == 14:
						# ATOM   2530  OD1 AS4 A 178  -24 969  -8 049  11 709  1 00  6 51  O
						debug = numbs[3]
					elif len(numbs) == 13:
						# ATOM   2520  HB2 ASN A 177  -29 165  -2 583  15 098  1 00 11 26    H
						debug = numbs[2]
					elif len(numbs) == 12:
						# ATOM   2518  H   ASN A 177   -26 911  -2 596  15 780  1 00  8 29   H 
						debug = numbs[1]
				else:
					# ANISOU    1  N   THR A  29     2101   4075   3166  -1332     57  -1083 N
					if len(numbs) == 8:
						debug = numbs[1]
					else:
						debug = numbs[2]

				if debug == line[5]:
					# ATOM   3389  H  BSER A 245
					bool4 = len(line[3]) > 3 and line[3][-3:] in res_ids and int(debug) in res_ids[ line[3][-3:] ]
					bool5 = len(line[3]) == 3 and line[3] in res_ids and int(debug) in res_ids[ line[3] ]  
					if bool4 or bool5:
						mut_temp.append( i )
				elif debug == line[4]:
					# ATOM   3392  HB2ASER A 245
					bool4 = line[2][-3:] in res_ids and int(debug) in res_ids[ line[2][-3:] ]
					if bool4:
						mut_temp.append( i )

		# self.pdb_str_not2Mut is every line from self.pdb_string wich is not related to mutation possibilities 
		self.pdb_str_not2Mut = []
		for i in self.pdb_string:
			if i not in mut_temp:
				self.pdb_str_not2Mut.append(i)
		self.pdb_str_2Mut = mut_temp
		res_Flag=''
		prev_res = res_Flag
		res_number='0'
		for i in range(len(self.pdb_str_2Mut)):
			line=self.pdb_str_2Mut[i].split()
			bool1 = 'ATOM' in line[0] or 'ANISOU' in line[0] or 'HETATM' in line[0]
			if len(line) > 3:
				bool2 = 'HOH' not in line[3] and 'na' not in line[3].lower() and 'cl' not in line[3].lower()
			else:
				bool2 = True
			if bool1 and bool2:
				numbs = re.findall(r'[0-9]+',self.pdb_str_2Mut[i])
				debug = '' # debug is the res id counter
				# Correcting the numbs problem that occurs when we have:
				if 'ATOM' in line[0] or 'HETATM' in line[0]:
					if len(numbs) == 14:
						# ATOM   2530  OD1 AS4 A 178  -24 969  -8 049  11 709  1 00  6 51  O
						debug = numbs[3]
					elif len(numbs) == 13:
						# ATOM   2520  HB2 ASN A 177  -29 165  -2 583  15 098  1 00 11 26    H
						debug = numbs[2]
					elif len(numbs) == 12:
						# ATOM   2518  H   ASN A 177   -26 911  -2 596  15 780  1 00  8 29   H 
						debug = numbs[1]
				else:
					# ANISOU    1  N   THR A  29     2101   4075   3166  -1332     57  -1083 N
					if len(numbs) == 8:
						debug = numbs[1]
					else:
						debug = numbs[2]
				bool3 = debug == line[5]
				#'debug' is the res counter
				if bool3:
					# ATOM   3389  H  BSER A 245
					if debug != res_number:
						# First residue of pdb or a different residue from the line above
						res_number = debug
						prev_res = res_Flag
						res_Flag = line[3]
						if res_Flag not in self.mut_pos:
							# Residuo first appearance
							# Phase A (new residue)
							# self.mut_pos := {'residueA': [(pos_i,pos_i)]} # pos_i is the first position, in self.pdb_str_2Mut, of residueA 
							self.mut_pos[res_Flag] = [(i,i)]
						else:
							# Phase C (not the same residue as previous iteration)
							# prev_res in self.aminoacid_atmcount is False in the first iteration and True elsewhere
							if prev_res in self.aminoacid_atmcount:
								at1, at2 = self.mut_pos[prev_res][len(self.mut_pos[prev_res])-1]
								# 'at1' is the first position of 'prev_res' on self.pdb_str_2Mut and 'at2' is the LAST KNOWN position 
								if at2 - at1 +1 < self.aminoacid_atmcount[prev_res]:
									# Doesn't have the minimal number of atoms
									if prev_res not in self.trash:
										# self.trash := {'prev_res': [(pos_iniA, pos_finA)]}
										self.trash[prev_res] = [self.mut_pos[prev_res].pop()]
									else:
										# self.trash := {'prev_res': [(pos_iniA, pos_finA),...,(pos_iniX, pos_finX)]}
										self.trash[prev_res].append(self.mut_pos[prev_res].pop())
							# Another (or "new") residue to consider from self.pdb_str_2Mut
							self.mut_pos[res_Flag].append((i,i))
					else:
						# Phase B (in the same residue as previous iteration)
						x,y = self.mut_pos[res_Flag][len(self.mut_pos[res_Flag])-1] #Tuples are immutable
						# y must be the LAST KNOWN position of the residue
						y=i
						# self.mut_pos := {'residueA': [(x,y)]}
						self.mut_pos[res_Flag][len(self.mut_pos[res_Flag])-1]=(x,y)
				elif debug == line[4]:
					# ATOM   3392  HB2ASER A 245
					if debug != res_number:
						res_number = debug
						res_Flag = line[2][3:]
						if res_Flag not in self.mut_pos:
							# Phase A
							self.mut_pos[res_Flag]=[(i,i)]
						else:
							# Phase C
							if prev_res in self.aminoacid_atmcount:
								at1, at2 = self.mut_pos[prev_res][len(self.mut_pos[prev_res])-1]
								if at2 - at1 +1 < self.aminoacid_atmcount[prev_res]:
									if prev_res not in self.trash:
										self.trash[prev_res] = [self.mut_pos[prev_res].pop()]
									else:
										self.trash[prev_res].append(self.mut_pos[prev_res].pop())
							self.mut_pos[res_Flag].append((i,i))
					else:
						# Phase B
						x,y = self.mut_pos[res_Flag][len(self.mut_pos[res_Flag])-1]
						y=i
						self.mut_pos[res_Flag][len(self.mut_pos[res_Flag])-1]=(x,y)

		# x, y positions from self.pdb_str_2Mut, which is given by 'res_ids'
		# self.mut_pos := {'residueA': [(xa1,ya1), ..., (xan,yan)], ..., 'residueZ': [(xz1,yz1), ..., (xzn,yzn)]}
		# Filling self.res_atm
		self.atm_detail()
		# self.res_atm := {'MET': [(0th appearance,{'C':[line 32 of self.pdb_str_2Mut, line 35 of self.pdb_str_2Mut,...]})]}
		
		return (len(self.mut_pos), len(self.pdb_str_2Mut))

	def atm_detail(self):
		'''Identify the pdb_str_2Mut data position of each atom in every residue in mut_pos (this makes easier to count atoms for each residue separately).'''

		for i in self.mut_pos:
			# i := a residue type (GLU,ALA,TRP,...)
			for j in range(len(self.mut_pos[i])):
				# j := jth appearence of residue 'i' in the system
				res_atm1, res_atm2 = self.mut_pos[i][j]
				for z in range(res_atm1,res_atm2 + 1):
					line = self.pdb_str_2Mut[z].split()
					# line[len(line)-1] is the atom and line[1] the atom number in the file. 'z' represents the line which you'd find this if the PDB had no comments.
					if i not in self.res_atm:
						# res_atm format := {'MET': [(0th appearence,{'C':[line 32 of self.pdb_str_2Mut, line 35 of self.pdb_str_2Mut,...]})]}
						self.res_atm[i]=[(j,{line[len(line)-1]:[z]})]
					elif j == self.res_atm[i][len(self.res_atm[i])-1][0]:
						if line[len(line)-1] not in self.res_atm[i][len(self.res_atm[i])-1][1]:
							# "new" atom of res i
							self.res_atm[i][len(self.res_atm[i])-1][1][line[len(line)-1]]=[z]
						else:
							self.res_atm[i][len(self.res_atm[i])-1][1][line[len(line)-1]].append(z)
					else:
						# Residuo 'i' is in res_atm and last position of self.res_atm[i] is not about the residue you are treating, i.e., the res type "i" already appeared somewhere before.
						self.res_atm[i].append( (j,{line[len(line)-1]:[z]}) )

	def choose_mutation(self, mutation_type = 1):
		'''Choose randonly an aminoacid from the active site and make a mutation based on the flag 'mutation_type' (that was not already made before).

Returns -1 if all probable mutations were already made for the chosen mutation_type.

Returns -2 if you choose an unknwon mutation_type.

Returns -3 if something went wrong with the mutation (You'll need to evaluate if its ok to just continue without this mutation or not - the problem may or not happen with other combinations so, it's always better to look what happened).

Returns 0 if everything is OK.

Parameters
----------

mutation_type:
	1 - One mutation from the attribute 'mutation_links'.
	2 - One mutation from the attribute 'improbable_mutation'. 
	3 - Mutations of choice (not restricted to active site).'''
		
		# Filling self.mut_pos and self.res_atm
		if mutation_type == 1 or mutation_type == 2: #!= 3
			# self.active_site := {'SER':[newid, ...], ...}
			self.res_detail(res_ids=self.active_site)
			self.mutation_temp_prob   = cp(self.mutation_links)
			self.mutation_temp_improb = cp(self.improbable_mutation)
		elif mutation_type == 3:
			# self.chosen_mutation := {('SER', 'newid'): 'MET', ('SER', 'O1newid'): 'MET'}
			# self.chosen_residues := {'SER':[newid, ...], ...} # directly related to self.chosen_mutation 
			self.res_detail(res_ids=self.chosen_residues)
		else:
			print("Please use 1, 2 or 3 for the argument mutation_type!")
			return -2

		self.mut_count = 0
		if mutation_type == 3:
			multiple = False
			if len(self.chosen_mutation) > 1:
				multiple = True
			for t in self.chosen_mutation:
				# self.chosen_mutation := {('SER', 'newid'): 'MET', ('SER', 'O1newid'): 'MET', ('AS4', 'newid'): 'ALA'}
				# self.chosen_residues := {'SER':[newid, ...], ...} # directly related to self.chosen_mutation 
				# self.beforeAmber_res2mut := {'SER':[old_id, ...], 'ASP':[old_id, ...], ...}
				# both self.chosen_residues and self.beforeAmber_res2mut are ordered
				# self.mut_pos := {'AS4': [(xa1,ya1), ..., (xan,yan)], ..., 'residueZ': [(xz1,yz1), ..., (xzn,yzn)]}
				self.mut_count += 1
				if self.mut_count == len(self.chosen_mutation):
					# This triggers if only one mutation will be made
					#  or if its the last of multiple multations
					multiple = False
				for t_r in range( len(self.chosen_residues[t[0]]) ):
					# t_r := position in the list of aminoacids t[0]
					if self.chosen_residues[t[0]][t_r] == int(t[1]):
						break
				olddies = cp(t[0])
				if olddies == 'AS4':
					olddies = 'ASP'
				elif olddies == 'HIP':
					olddies = 'HIS'
				elif olddies == 'GL4':
					olddies = 'GLU'

				verif = self.mutation(res=t[0], old_res_id="%s_%s"%(olddies, str(self.beforeAmber_res2mut[olddies][t_r])),
				pos_id=t_r, pos=self.mut_pos[t[0]][t_r], mut_res=self.chosen_mutation[t], mult_flag=multiple)

				if verif == -1:
					print("Go check what went wrong with the mutation %s to %s"%(t[0], self.chosen_mutation[t]))
					return -3
		else:
			# mutation_type != 3
			# Choosing a random aminoacid from the active site to mutate
			res_chosen = random.sample(self.mut_pos.keys(),1)[0]
			res_not_amber = res_chosen
			if res_chosen == 'AS4':
				res_not_amber = 'ASP'
			elif res_chosen == 'HIP':
				res_not_amber = 'HIS'
			elif res_chosen == 'GL4':
				res_not_amber = 'GLU'

			if mutation_type == 1:
				# Making a "probable" mutation
				temp_mutations = self.mutation_temp_prob
			elif mutation_type == 2:
				# Making an "improbable" mutation
				temp_mutations = self.mutation_temp_improb
		
			if res_not_amber in temp_mutations and len(temp_mutations[res_not_amber]) > 0:
				# There is a mutation for the res_chosen. So, for simplicity, take the last mutation on the list (If you are going to test multiple mutations, it's a good thing not to try the same thing more than once).
				# x, y positions from self.pdb_str_2Mut, which is given by 'res_ids'
				# self.mut_pos := {'residueA': [(xa1,ya1), ..., (xan,yan)], ..., 'residueZ': [(xz1,yz1), ..., (xzn,yzn)]}
				if len(self.mut_pos[res_chosen]) == 1:
					# If there is only one res_chosen in the self.pdb_str_2Mut
					# Popping the mutation from the list and applying it (next time the process won't see this mutation as possible).
					for j in range(len(temp_mutations[res_not_amber])):
						mut_chosen = temp_mutations[res_not_amber].pop()
						# REMINDER - AFTER AMBER:
						# self.chosen_residues #self.active_site := {'SER':[newid, ...], ...}
						# mut_chosen might have been restricted, so this loop will go on until it finds a mutation not prohibited (by self.mutations_restricted) or stop after finding none
						res_chosen_flag = '%s_%d'%(res_not_amber,self.beforeAmber_res2mut[res_not_amber][0])
						if res_chosen_flag not in self.mutations_restricted:
							# self.mutations_restricted eg:= {'SER_160': ['GLY']}
							# This residue was never chosen for a mutation
							verif = self.mutation(res=res_chosen, old_res_id=res_chosen_flag, pos_id=0, pos=self.mut_pos[res_chosen][0], mut_res=mut_chosen)
							if verif == -1:
								print("Go check what went wrong with the mutation %s to %s"%(res_chosen,mut_chosen))
								return -3
							else:
								# Taking notes on which mutations were made
								# self.mutations_restricted[res_not_amber] = [(0,mut_chosen)]
								self.mutations_restricted[res_chosen_flag] = [mut_chosen]
								return 0
						elif mut_chosen not in self.mutations_restricted[res_chosen_flag]:
							# This residue was chosen before, but not with the mutation mut_chosen 
							verif = self.mutation(res=res_chosen, old_res_id=res_chosen_flag, pos_id=0, pos=self.mut_pos[res_chosen][0], mut_res=mut_chosen)
							if verif == -1:
								print("Go check what went wrong with the mutation %s to %s"%(res_chosen,mut_chosen))
								return -3
							else:
								# Taking notes on which mutations were made
								# self.mutations_restricted[res_not_amber].append( (0,mut_chosen) )
								self.mutations_restricted[res_chosen_flag].append( mut_chosen )
								return 0

				elif len(self.mut_pos[res_chosen]) > 1:
					'''This block seems fine and we could use it for the case len(self.mut_pos[res_chosen]) == 1, but since it wasn't tested yet it is too risky to assemble this if-elif block.'''
					# There may be more than one appearance of res_chosen in the active site in some families of enzymes
					mut_chosen = temp_mutations[res_not_amber].pop()
					for i in range(len(self.mut_pos[res_chosen])):
						appearance_chosen = self.mut_pos[res_chosen][i] #(xa,ya)
						# We may want to try the "same" mut SER_ALA (this "new" SER has a different res id ), so we need to deepcopy a temporary variable to make sure all combinations were done 
						temp_jmut = cp(temp_mutations[res_not_amber])
						for j in range(len(temp_mutations[res_not_amber])):
							mut_chosen = temp_jmut.pop()
							# This mut_chosen might have been restricted, so this loop will go on until it finds a mutation not prohibited (by self.mutations_restricted) or stop after finding none
							# self.beforeAmber_res2mut := {'SER':[old_id, ...], ...}
							res_chosen_flag = '%s_%d'%(res_not_amber,self.beforeAmber_res2mut[res_not_amber][i])
							# self.mutations_restricted eg:= {'SER_160': ['GLY']}
							if res_chosen_flag not in self.mutations_restricted:
								# This residue was never chosen for a mutation
								verif = self.mutation(res=res_chosen, old_res_id=res_chosen_flag, pos_id=i, pos=appearance_chosen, mut_res=mut_chosen)
								if verif == -1:
									print("Go check what went wrong with the mutation %s%d to %s"%(res_chosen,i,mut_chosen))
									return -3
								else:
									# Taking notes on which mutations were made
									self.mutations_restricted[res_chosen_flag] = [mut_chosen]
									# Remember, the original system will suffer only one mutation # if it's need to do multiple mutations use mutation_type == 3!
									return 0
							elif mut_chosen not in self.mutations_restricted[res_chosen_flag]:
								# {res_chosen_flag: []}
								verif = self.mutation(res=res_chosen, old_res_id=res_chosen_flag, pos_id=i, pos=appearance_chosen, mut_res=mut_chosen)
								if verif == -1:
									print("Go check what went wrong with the mutation %s%d to %s"%(res_chosen,i,mut_chosen))
									return -3
								else:
									# Taking notes on which mutations were made
									self.mutations_restricted[res_chosen_flag].append( mut_chosen )
									return 0

						# If the process ends up here at the final line of the "i-loop" block then it will look for another residuo to mutate
			else:
				# If the process ends here it's because all probable mutations were already made
				return -1

		return 0

	def mutation(self, res = 'SER', old_res_id = 'SER_160', pos_id = 0, pos = (0,10), mut_res = 'ALA', mult_flag = False):
		'''Changes the residues in 'res' to the one in 'mut_res' in the system (self.pdb).
		Returns 0 if everything seems ok (after using Leap on mult_flag == False).
		Returns 1 for multiple mutation working correctly.
		Returns -1 if there is any fatal error.

Parameters
----------

res: Residue type to mutate.

old_res_id: Residue id from the original PDB file (for naming purposes only)

pos_id: Order of appearence in the self.pdb_str_2Mut for the exact residue you want to mutate (in relation to its residue type only, ex: 0 -> first one; 1 -> second ...).

pos: Sequence of positions in self.pdb_str_2Mut which corresponds to the chosen residue atoms.

	OBS: Both 'pos' and 'pos_id' can be found in self.res_atm

mut_res: Residue you wish for the mutation.

mult_flag: (Boolean) If True it won't pass through Leap and just modify self.pdb_str_2Mut. It must use LEap (False option) only after all wanted modifications.
'''

		# Lines to be modified
		# self.pdb_str_2Mut is about the system-amber.pdb
		atom_list = cp(self.pdb_str_2Mut[pos[0]:pos[1]+1])
		temporary_files = []
		# self.res_atm := {'MET': [(0th appearance,{'C':[line 32 of self.pdb_str_2Mut, line 35 of self.pdb_str_2Mut,...]})]}
		# Can't modify directly the self.res_atm because it uses tuples..so you need to change a tuple for another
		temp_flag = False
		for i in self.res_atm[res]:
			if i[0] == pos_id:
				# To be sure it is the correct appearance
				temp_flag = True
				atom_pos = cp(self.res_atm[res][pos_id][1])
				break
		if not temp_flag:
			return -1
		for i in atom_pos:
			for j in range(len(atom_pos[i])):
				atom_pos[i][j] -= pos[0]
		# Now atom_pos holds the correct positioning of atoms in the atom_list
		
		modified_text = self.aminoacid_mutation(transform=mut_res, atm_list = atom_list, atm_pos= atom_pos)
		
		# modified_text probably hasn't the same size as the text before the mutation change!
		# It's necessary to update self.mut_pos and self.res_atm for multiple mutations!
		orig_size = pos[1] - pos[0] +1
		new_size = len(modified_text)
		diff_size = orig_size - new_size
		for i in self.mut_pos:
			if i != res:
				temp_i, temp_f = self.mut_pos[i][pos_id]
				temp_i -= diff_size
				temp_f -= diff_size
				self.mut_pos[i][pos_id] = (temp_i, temp_f)
		self.res_atm = {}
		self.atm_detail()

		# temp holds all the info from pdb_str_2Mut, mutated and non-mutated parts (in the random active site mutation there are residues which will not be modified) # see self.res_detail
		temp = cp(self.pdb_str_2Mut[:pos[0]])
		# In case pos[0] ==0 -> self.pdb_str_2Mut[:0] := []
		for i in modified_text:
			temp.append( i )
		# In case pos[1]+1 is out of range -> self.pdb_str_2Mut[pos[1]+1:] := []
		temp.extend(self.pdb_str_2Mut[pos[1]+1:])

		# If there are more mutations to make, it is not recommended to pass through Leap before making all modifications
		#  since if we use Leap in every iteration we'll acumulate erros from numerical aproximations.
		if mult_flag:
			self.pdb_str_2Mut = temp
			return 1
	
		# self.pdb_str_not2Mut is every line from self.pdb_string wich is not related to mutation possibilities 
		count_mut = 0
		# This is done to maintain the residue order
		for i in range(len(self.pdb_string)):
			if  self.pdb_string[i] != self.pdb_str_not2Mut[i]:
				#  self.pdb_string[i] is looking at the old non-mutated residuo
				count_mut = i
				break
		
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################

		# FAZER DEPOIS DA QUALIFICACAO!
		# numbs = re.findall(r'[0-9]+',self.pdb_str_2Mut[i])
		# len(numbs) == 14: 
		#  ATOM   2530  OD1 AS4 A 178  -24 969  -8 049  11 709  1 00  6 51  O
		# self.chosen_residues := {'SER':[newid, ...], ...} eg: {'ASP':[178], 'HIP': [209]}

		new_pdb_string = self.pdb_str_not2Mut[:count_mut] #self.pdb_string[:count_mut]
		new_pdb_string.extend(temp)
		new_pdb_string.extend(self.pdb_str_not2Mut[count_mut:])

		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################
		####################################################################################

		self.pdb_string = new_pdb_string 

		# Generating a new pdb with the mutation
		if not self.not_rand_mut:
			new_filename = '%s_mut_%s-%s.pdb'%(self.pdb[:-4],old_res_id,mut_res)	
		else:
			naming = ''
			name_count = 0
			for n in self.oldid_mut:
				name_count += 1
				if name_count > 1:
					naming += '__'
				self.newid_info.append( '\nOld:\n%s:\t%s\n'%(n[0],n[1]) ) 
				naming += '%s_%s-%s'%(n[0],n[1],self.oldid_mut[n])
			new_filename = '%s_mut_%s.pdb'%(self.pdb[:-4],naming)

		cmd('rm %s'%new_filename)
		f = open(new_filename,'w')
		for i in self.pdb_string:
			f.write(i)
		f.close()

		# Renumbering the atoms and preparing for constant pH MD #--reduce
		# Even though it looks unnecessary because we are not modifying the residue order, this prevents a huge error!
		cmd('pdb4amber -i %s -o %s_renum.pdb --constantph'%(new_filename,new_filename[:-4]))
		temporary_files.append( self.pdb )
		self.pdb = '%s_renum.pdb'%new_filename[:-4]

		''' In the following lines of this code we'll treat the 'charge' and 'check' commands of Leap, since with 'charge' we can count the charges and neutralize the system correctly after and 'check' verify if leap can build the topology (unknwon things like atoms, residues, or even numerical bugs can give FATAL errors which stops Leap from building the topology.)
		If in the "charge step" we save the system; and load and save again in the "check step" to finally load and "build" the topology, Leap won't work properly. Leap will give a FATAL error, even if there is none, it's interpretting whatever it reads like an unknown residue somewhere.

		This is propably happening because it uses a numerical approximation everytime a system is loaded and needs a check or whatever, so when you save the system, you save this approximation. If you reopen to check other things you'll approximate an approximation...numerical problems like this gives birth to random abnormalities.

		To overcome this we are not saving any leap-process untill everything is checked and saving only in the topology building step.''' 
		self.leap_in(checking=True, checking_cmd='charge', charge_fix=False, add_Na=0, add_Cl=0)
		cmd('rm leap.log')
		cmd('tleap -f tleap.in')
		
		# Checking the total charge, if it's not zero, the code below will neutralize the system adding ions
		tt = open('leap.log','r')
		temp = tt.readlines()
		tt.close()
		charge = 0.0
		for i in range(len(temp)):
			if 'Total unperturbed' in temp[-i]:
				data = temp[-i].split()
				charge = float(data[len(data)-1])
				break
		cha = int(charge)
		if charge - cha > 0.05:
			print("\nError with the charge\n")
			return -1
		elif cha != 0:
			if cha > 0:
				self.leap_in(checking=True, checking_cmd='check', charge_fix=True, add_Na=3, add_Cl=cha+3)
			else:
				self.leap_in(checking=True, checking_cmd='check', charge_fix=True, add_Na=cha+3, add_Cl=3)

			# This will check for any fatal errors, proximity warnings don't matter (minimization can fix it). And in order to check for the true error we need to delete all leap log before attempting the leapfix
			cmd('rm leap.log')
			cmd('tleap -f tleap.in')
			# In the check system, we are looking for something troubling like this:
			# FATAL:  Atom .R<ALA 100>.A<HG2 11> does not have a type.
			# With this kind of error it's impossible to build the topology 		
			f = open('leap.log','r')
			temp = f.readlines()
			error = []
			for i in temp:
				if 'FATAL:' in i:
					error.append(i)
			if len(error) >0:
				for i in error:
					print(i[:-1])
				return -1

		# If the process survived until here, then everything seems ok and we now need to create the topology
		if cha > 0:
			self.leap_in(checking=False, charge_fix=True, add_Na=3, add_Cl=cha+3)
		else:
			self.leap_in(checking=False, charge_fix=True, add_Na=cha+3, add_Cl=3)

		cmd('tleap -f tleap.in')
		temporary_files.append( 'tleap.in' )
		self.pdb = 'systemLeap.pdb'

		# Now we can delete all transitory files
		cmd('ls > temporary_stage_ls.out')
		for i in temporary_files:
			if '.py' in i:
				continue
			else:
				cmd('rm %s'%i)

		return 0

	def aminoacid_mutation(self, transform='MET', atm_list = ['ATOM 1432 N MET A ...','ATOM 1433 CA MET A ...'],
	atm_pos={'N':[0], 'C': [1]}):
		'''Edits the atm_list.

Parameters
----------
transform: Any aminoacid.

atm list: Sequence of lines, from the pdb file, to be modified for the mutation.

atm pos: Atoms' positions respectively to the sequence in atm_list.
'''
		# self.restypeatoms format:
		#'ALA': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB1', 'HC', 'HB2', 'HB3', 'C', 'O']
		#'GLY': ['N', 'H', 'CA', 'CX', 'HA2', 'H1', 'HA3', 'C', 'O']
		#'ASP': ['N', 'H', 'CA', 'CX', 'HA', 'H1', 'CB', 'CT', 'HB2', 'HC', 'HB3', 'CG', 'C', 'OD1', 'O2', 'OD2', 'O']

		new_list = []
		for i in atm_pos:
			# Loop chose an element
			elem_pos = []
			for j in atm_pos[i]:
				# Loop chose a position
				temp = atm_list[j].split()
				if temp[2] in self.restypeatoms[transform]:
					# Ex: if 'HB2' in restypeatoms['MET'], then add to the list
					# Before confirming this atom to new_list, we'll change the residue, temp[3], to the 'transform' variable residue
					data_pos = atm_list[j].find(temp[3],0,len(atm_list[j]))
					temp2 = atm_list[j][:data_pos]
					temp2 += transform
					temp2 += atm_list[j][data_pos+3:]
					elem_pos.append( temp2 )

			# Maybe now there are some missing atoms, but there is no strange atom names
			new_list.extend( elem_pos )
 
		if self.mut_count > 1:
			info = stat('mutation_made.out')
			f    = open('mutation_made.out', 'r+')
			f.seek(info.st_size)
			f.write('\n')
			f.write('-'*50)
			f.write('\nOld residue:\n\n')
		else:
			cmd('rm mutation_made.out')
			f = open('mutation_made.out', 'w')
			f.write('Old residue:\n\n')
		for i in atm_list:
			f.write(i)
		f.write('\nModification before leap:\n\n')
		for i in new_list:
			f.write(i)
		f.close()

		return new_list


if __name__ == "__main__":
	import sys
	arg = sys.argv
	version = '''Manager Alpha - Molecular Dynamics simulation manager (using Amber and Ambertools).

Manager Beta:
	*Information_cycle bug corrected.
	*Tutorial guide made which show how the code works and asks for each variable specifying the formats.
	*Simulation can now be made directly on Terminal console.

Manager Version 1.0:
	*Added option for explicit solvent.
	*Method 'Amber_run.simulation' now verifies titratable residues and applies cpinutil correctly.
	*'Amber_run' class rewritten for CpHMD.
	*Added mpi options.
	*Code better commented.
	*__main__ compilation fixed.
	*Tutorial guide removed and created a better help message.
	*Added method 'Amber_run.leap_exec' so the code now tries to create mol2 files if the Leap libraries can't deal with your system.
	*Fixed bug in the explicit solvent MD run.
	*New mutation option:
		**Added attributes to constructor to be able to mutate.
		**Amber_mutation class added
		**Added random mutation method at the active site (for enzymes).
		**Added a mutation of choice option.
		OBS: Simulation with mutation needs to run with explicit solvent, if none given the code will use water by default.
	*Corrected methods:
		**'Amber_par.cpptraj_in';
		**'Amber_run.analysis'.
		**Methods 'Amber_run.sander_run', 'Amber_run.gpu_run' and 'Amber_run.pmemd_only' merged together in the 'Amber_run.minimization_to_analysis'.
	*Fixed explicit solvent CpHMD BUG, by modifying the following methods from 'Amber_par':
		**'input_min';
		**'input_heat';
		**'input_equil';
		**'cpptraj_in' - RMSD, RMSF and RADGYR are correctly done now for CpHMD.

Manager Version 1.1:
	*Separation of the simulation.sh in stages to minimize errors from external sources. 
	*Annealing duration shortened because it was unnecessarily long.
	*Mutation method fixed.
	*Error report file created.
	*An easy restart .sh of current simulation stage created.

Manager Version 1.2:
	*Added a custom option for the duration of each stage of simulation 


Copyright (C) 2021  Braga, B. C. 
e-mail: bruno.braga@ufms.br 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

	# Default keys
	version_n     = '1.2'
	inst_only     = False
	version_only  = False
	arqui         = 'gpu'
	mode          = 'L'
	custom_mode   = False
	custom_min    = 10**4
	custom_a      = 10**6
	custom_e      = 10**6
	custom_p      = 10**6
	gpu_unit      = 0
	icyc          = 200
	sim_goal      = 'MD'
	sim_pdb       = 'WillThisWork.pdb'
	mpi_flag      = False
	cores_        = 2
	explicit_flag = False
	solvent       = 'water'
	cuttoff_      = 12.0
	prt_res       = 'AS4 SER HIP'
	atv_site      = {'SER':[160],'HIS':[237],'ASP':[206]}
	mut_rand      = False
	mut_choice    = False
	mut           = []
	mut_type      = 1
	mut_done      = {}
	prep_only     = False
	multiple_pH   = False
	pH_ini        = 7.0
	initial_ph    = pH_ini
	ini_pH        = pH_ini
	fin_pH        = ini_pH
	pHInterval    = 1.0

	# If for some reason you don't want to use the HMR method, set the attribute "self.hmr" to False

	flags = ["&","-v","--version","-h","--help",'-arq','gpu','-mode','-res','-mut','-rdmut','-atsite','-rdydone','-icyc','-solv','-g','ph','-i','-mpi','-explicit','-prepstp','-cut','-phdset','MIN','A','E','P']
	
	# Flag verification
	for i in arg:
		if i[0] == '-' and i.lower() not in flags:
			print("Unkown Flag used: ", i)
			inst_only = True

	cut   = 0 # Counter for input flags
	for i in range(len(arg)):
		if inst_only:
			break
		if cut == i:
			if arg[i].lower() == "&":
				break
			elif arg[i].lower() == "-v" or arg[i].lower() == "--version":
				version_only = True
				print("Current version: %s\n"%version_n)
				print(version)
				break
			elif arg[i].lower() == "-h" or arg[i].lower() == "--help":
				inst_only = True
				print("Welcome to Amber Simulation Manager %s:\n"%version_n)
				print("Copyright (C) 2021  Braga, B. C.\nThis program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; use option '-v' for details.\n")
				print("\nUsage:\n\t-v or --version\t\tPrints current version and its corrections relating previous versions.\n")
				print("\t-h or --help\t\tPrints this message.\n")
				print("\t-i\t\tInput PDB file.\n")
				print("\t-g\t\tGoal as MD or CpHMD. Default: MD.\n\t\t\tOBS: If you choose CpHMD and don't use flag 'res' the code will choose by default to tritate %s.\n"%prt_res)
				print("\t-ph\t\tProduction pH. Default: 7.0.\n")
				print("\t-phdset\t\t(This option has priority over the -ph option) Defines the pH list (only if you'll do multiple pH CpHMD production). After this flag an initial pH, a final pH and an interval unit must be given in this order (Ex: -phdset 4.0 7.0 0.5). Which represents pH:[4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0].\n")
				print("\t-res\t\tResidues to tritate, for CpHMD (which must be informed right after this flag).\n")
				print("\t-rdmut\t\tRandom active site mutation.\n\t\t\t1:random aminoacid mutates to a similar one (eg. SER_160-THR). \n\t\t\t2:random aminoacid mutates to one not so similar (eg. SER_160-MET).\n")
				print("\t-atsite\t\t(Auxiliary option for -rdmut) Active site - for enzymes only. Must be informed right after this flag as eg: -atsite SER_160 HIS_237 ASP_206. Necessary only if a mutation will be done. Defaut: set as petase's active site.\n")
				print("\t-rdydone\t(Auxiliary option for -rdmut) Informs which mutation you're restricting from the \"random\" choice. Eg. -rdydone SER_160-MET ASP_206-GLY ASP_206-ILE.\n")
				print("\t-mut\t\tMutation of choice. Eg: -mut SER_160-MET.\n\t\t\t\tObs: For the mutation options (-mut and -rdmut) a solvent must be used. If none given, water will be used by default.\n")
				print("\t-mode\t\tDuration of the simulation's stages (in time steps).\n\t\t\tL: Annealing: 10**4 | Equilibration: 10**4 |  Production: 10**5\n\t\t\tM: Annealing: 10**4 | Equilibration: 10**5 |  Production: 10**6\n\t\t\tGM: Annealing: 10**5 | Equilibration: 10**6 |  Production: 5*10**6\n\t\t\tGH: Annealing: 5*10**5 | Equilibration: 5*10**6 |  Production: 10**7\n\t\t\tGU: Annealing: 10**6 | Equilibration: 10**7 |  Production: 5*10**7\n\t\t\tCT: Minimization: MIN | Annealing: A | Equilibration: E |  Production: P | Ex: -mode CT MIN 10**4 A 10e6 E 10e7 P 10e9\n")
				print("\t-explicit\tExplicit solvent will be used. One of the following options must be given after this flag:\n\t\t\twater\n\t\t\tmethanol\n\t\t\tchloroform\n\t\t\tN-methyacetamide\n\t\t\turea\n")
				print("\t-icyc\t\t(Integer between 50 and 1000) Number of times information regarding energy, rmsd, rmsf, coordinates, etc, are saved for each simulation stage (Default: 200).\n\t\t\t\tObs: For -mode GU, this is set between 800 and 1000.\n")
				print("\t-arq\t\tChoice in architecture\n\t\t\tgpu (if chosen, the GPU id must be informed right after this flag)\n\t\t\tsander\n\t\t\tpmemd\n")
				print("\t-mpi\t\tMulticore run (the number of cores must be informed right after this flag). Default compiler: mpiexec.\n\t\t\t\tObs: Not implemented for -arq gpu.\n")
				print("\t-cut\t\tNonbonded cutoff in angstrom (Default: 12).\n")
				print("\t-prepstp\t(Option for Devs.) Stops the run after all preparations are complete (right before minimization).\n\t\t\t\tObs: The shell script, to run the simulation, will still be created for you.\n")
				print("\nExamples:\n\t$ python3 ASM.py -i 1UBQ.pdb -arq sander -mpi 4 -mode L -g CpHMD -phdset 6.0 7.0 0.5\n")
				print("\n\t$ python3 ASM.py -i 6eqe.pdb -arq pmemd -mpi 4 -mode L -g CpHMD -ph 6.5 -explicit water\n")
				print("\n\t$ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode GM -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -rdmut 1\n")
				print("\n\t$ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode GM -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -mut SER_160-MET\n")
				print("\n\t$ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode GM -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -mut SER_160-MET HIS_237-ALA\n")
				print("\n\t$ python3 ASM.py -i paracetamol.pdb -arq sander -mode M\n")
				print("\n\t$ python3 ASM.py -i 6eqe.pdb -arq sander -mpi 4 -mode L -g MD -explicit water\n")
				break

			# Key verifications
			elif arg[i].lower() == '-arq':
				arqui = arg[i+1]
				if arqui.lower() == 'gpu':
					gpu_unit = arg[i+2]
				continue
			elif arg[i].lower() == '-mode':
				mode = arg[i+1]
				if mode.upper() == 'CT':
					custom_mode = True
					cc = i+2
					while cc < len(arg):
						if arg[cc].upper() == 'MIN':
							cc += 1
							custom_min = IntTransNumb(arg[cc])
						elif arg[cc].upper() == 'A':
							cc += 1
							custom_a = IntTransNumb(arg[cc])
						elif arg[cc].upper() == 'E':
							cc += 1
							custom_e = IntTransNumb(arg[cc])
						elif arg[cc].upper() == 'P':
							cc += 1
							custom_p = IntTransNumb(arg[cc])
						else:
							i = cc-1 
							break
						cc += 1
				continue
			elif arg[i].lower() == '-res':
				res_i = arg[i+1]
				cc = i+2
				while cc < len(arg):
					if arg[cc] not in flags:
						res_i += arg[cc]
					else:
						# This is a 'for-loop' so 'i' will increase at the end of iteration
						i = cc-1 
						break
					cc += 1
				prt_res = res_i
				continue
			elif arg[i].lower() == '-atsite':
				site_i = [ arg[i+1] ]
				cc = i+2
				while cc < len(arg):
					if arg[cc] not in flags:
						site_i.append( arg[cc] )
					else:
						# This is a 'for-loop' so 'i' will increase at the end of iteration
						i = cc-1 
						break
					cc += 1
				# site_i := ['SER_160', 'HIS_237', 'ASP_206']
				atv_site= {}
				for ji in site_i:
					pos_ji = ji.find('_',0,len(ji))
					# {'SER':[160],'HIS':[237],'ASP':[206]}
					if ji[:pos_ji].upper() not in atv_site:
						atv_site[ ji[:pos_ji].upper() ] = [int(ji[pos_ji + 1:])]
					else:
						atv_site[ ji[:pos_ji].upper() ]. append( int(ji[pos_ji + 1:]) )	
				continue
			elif arg[i].lower() == '-rdmut':
				mut_rand = True
				mut_type = int(arg[i+1])
				if mut_type not in [1,2]:
					print("\nINVALID ARGUMENT FOR MUTATION TYPE!\n")
					inst_only = True
					break
				continue
			elif arg[i].lower() == '-mut':
				mut_choice = True
				mut = [ arg[i+1] ]
				# mut := ['SER_160-MET']
				cc = i+2
				while cc < len(arg):
					if arg[cc] not in flags:
						mut.append( arg[cc] )
						cc += 1
					else:
						i = cc-1 
						break
				continue
			elif arg[i].lower() == '-rdydone':
				r_done = [ arg[i+1] ]
				cc = i+2
				while cc < len(arg):
					if arg[cc] not in flags:
						r_done.append( arg[cc] )
						cc += 1
					else:
						i = cc-1 
						break
				# r_done := ['SER_160-MET', 'ASP_206-GLY', 'ASP_206-ILE']
				# mut := {'SER_160': ['GLY']}
				mut_done = {}
				for ji in r_done:
					pos_ji = ji.find('-',0,len(ji))
					if ji[:pos_ji] not in mut_done:
						mut_done[ ji[:pos_ji] ] = [ji[pos_ji + 1:]]
					elif ji[pos_ji + 1:] not in mut_done[ ji[:pos_ji] ]:
						mut_done[ ji[:pos_ji] ].append( ji[pos_ji + 1:] )
				continue
			elif arg[i].lower() == '-icyc':
				icyc = int(arg[i+1])
				continue
			elif arg[i].lower() == '-solv':
				solvent = arg[i+1]
				continue
			elif arg[i].lower() == '-g':
				sim_goal = arg[i+1]
				continue
			elif arg[i].lower() == '-ph':
				initial_ph = float(arg[i+1])
				continue
			elif arg[i].lower() == '-phdset':
				multiple_pH = True
				ini_pH = float(arg[i+1])
				fin_pH = float(arg[i+2])
				pHInterval = float(arg[i+3])
				continue
			elif arg[i].lower() == '-i':
				sim_pdb = arg[i+1]
				continue
			elif arg[i].lower() == '-mpi':
				mpi_flag = True
				cores_ = arg[i+1]
				continue
			elif arg[i].lower() == '-explicit':
				explicit_flag = True
				solvent_ = arg[i+1].lower()
				if solvent_ not in ['water','methanol','chloroform','n-methyacetamide','urea']:
					print('Error: Unknown solvent chosen. Please use one of the following:')
					for solv in ['water','methanol','chloroform','N-methyacetamide','urea']:
						print(solv)
					inst_only = True
					break
				continue
			elif arg[i].lower() == '-cut':
				cuttoff_ = arg[i+1]
				continue
			elif arg[i].lower() == '-prepstp':
				prep_only = True
				continue
			cut +=1
		else:
			# cut!= i means that the current arg[i] was used in the previous iteration
			cut = i+1

	if mode == 'GU' and icyc < 800:
		icyc = 800
				
	if not inst_only and not version_only:
		if sim_pdb == 'WillThisWork.pdb':
			print("No system choosen!\nUse: python3 ASM.py --help\n")
		else:
			if not multiple_pH:
				pH_ini = initial_ph
				pH_fin = pH_ini 
			else:
				pH_ini = ini_pH
				pH_fin = fin_pH 

			if not mut_rand and not mut_choice:
				objeto = Amber_run(system=sim_pdb, simulation=sim_goal, pHstep=pHInterval, simulation_length=mode, custom_mode = custom_mode, custom_min = custom_min, custom_a = custom_a, custom_e = custom_e, custom_p = custom_p, pH=pH_ini, exp_solv= explicit_flag, solvent_in=solvent, prot_res=prt_res, box_cuttoff=cuttoff_, mpi_use= mpi_flag, mpicores=cores_, information_cycles=icyc, gpu=gpu_unit, prep_stop= prep_only)
				objeto.simulation(arq=arqui, ph_range=[pH_ini,pH_fin])
			elif mut_rand and not mut_choice:
				objeto = Amber_mutation(system=sim_pdb, simulation=sim_goal, pHstep=pHInterval, simulation_length=mode, custom_mode = custom_mode, custom_min = custom_min, custom_a = custom_a, custom_e = custom_e, custom_p = custom_p, pH=pH_ini, exp_solv= explicit_flag, solvent_in=solvent, prot_res=prt_res, box_cuttoff=cuttoff_, mpi_use= mpi_flag, mpicores=cores_, information_cycles=icyc, gpu=gpu_unit, prep_stop= prep_only, active_site=atv_site)
				objeto.mutations_restricted = mut_done
				if objeto.chosen_residues != -1:
					# There wasn't a problem with the mutation residues
					# objeto.chosen_mutation != {} 
					okMutation_flag = objeto.choose_mutation(mutation_type = mut_type)
					if okMutation_flag == 0:
						objeto.simulation(arq=arqui, ph_range=[pH_ini,pH_fin])
					elif okMutation_flag == -3:
						print("\nSOMETHING WENT WRONG WITH THE CHOSEN MUTATION! PLEASE CHECK THE FILE LEAP.LOG!\n")
					elif okMutation_flag == -1:
						print("\nAll possible mutations were already made for the chosen mutation type.\n")
			elif not mut_rand and mut_choice:
				# mut := ['SER_160-MET']
				objeto = Amber_mutation(chosen_mut_flag = True, chosen_mut=mut, system=sim_pdb, simulation=sim_goal, pHstep=pHInterval, simulation_length=mode, custom_mode = custom_mode, custom_min = custom_min, custom_a = custom_a, custom_e = custom_e, custom_p = custom_p, pH=pH_ini, exp_solv= explicit_flag, solvent_in=solvent, prot_res=prt_res, box_cuttoff=cuttoff_, mpi_use= mpi_flag, mpicores=cores_, information_cycles=icyc, gpu=gpu_unit, prep_stop= prep_only)
				objeto.mutations_restricted = {}
				if objeto.chosen_mutation != {}:
					# There wasn't a problem with the mutation residues
					okMutation_flag = objeto.choose_mutation(mutation_type = 3)
					if okMutation_flag == 0:
						objeto.simulation(arq=arqui, ph_range=[pH_ini,pH_fin])
					elif okMutation_flag == -3:
						print("\nSOMETHING WENT WRONG WITH THE CHOSEN MUTATION! PLEASE CHECK THE FILE LEAP.LOG!\n")
					elif okMutation_flag == -1:
						print("\nAll possible mutations were already made for the chosen mutation type.\n")
			else:
				print("\nOptions -mut and -rdmut can not be used in the same instance!\n")
