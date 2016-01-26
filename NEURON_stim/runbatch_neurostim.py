# runbatch.py -- Python script for running a batch of simulation 
#   multi-seed runs on a single machine. 
#
# Use
#   python runbatchpbs_mist_20params.py 
# from the command line to start the batch.
#
# eg. python runbatchpbs_mist_20params.py 14sep16b data/14sep16b_evol/gen_86_cand_151 
# Last update: 13/09/01 (salvadord)

# Set whether the commands are actually done (as opposed to just displayed).
do_commands = True

# Import packages needed for commands.
import os, sys
from numpy import *
from pylab import *
import subprocess
from popen2 import popen2
import time
import pickle
import random as rnd
from shutil import copyfile

# save list of perturbed cells in file
def makePerturbCells (filestem, percent):
	cellsTot = 191
	cellsPerturbed = int(round(percent/100.0 * cellsTot))
	randCells = rnd.sample(xrange(cellsTot), cellsPerturbed)
	outfilename = "%s-cellsPerturbed.txt"%(filestem) 
	savetxt(outfilename,randCells,fmt='%11.3f')
	print "Cell perturbation file %s saved." % outfilename

def copyPerturbCellsFile(fileorigin,filestem):
	copyfile("%s-cellsPerturbed.txt"%(fileorigin), "%s-cellsPerturbed.txt"%(filestem))

# save list of perturbed synapses in file
def makePerturbSyns (filestem, percent):
	# check max number of connections
	synsTot = 21588
	synsPerturbed = int(round(percent/100.0 * synsTot))
	randCells = rnd.sample(xrange(synsTot), synsPerturbed)
	outfilename = "%s-synsPerturbed.txt"%(filestem) 
	savetxt(outfilename,randCells,fmt='%11.3f')
	print "Synapse perturbation file %s saved." % outfilename

def copyPerturbSynsFile(fileorigin,filestem):
	copyfile("%s-synsPerturbed.txt"%(fileorigin), "%s-synsPerturbed.txt"%(filestem))

# save list of stimulated cells in file
def makeProbingStims (filestem,mistCell,mistStart,mistDuration,mistRate):
	 # Create full stimulus
	event1=[mistCell,mistStart,mistRate] # create data to write
	event2=[mistCell,mistStart+mistDuration,0] # 
	events=vstack([event1,event2])
	#event1=arange(1,10)
	outfilename="%s-in.txt"%(filestem) 
	savetxt(outfilename,events,fmt='%11.3f')
	print "Stimulus file %s saved." % outfilename
	
# save list of randomly stimulated cells in file
def makeMultProbingStims (filestem, mistCell,mistStart,mistDuration,mistRate):
	maxCells = 10 # max number of cells stim simultaneously
	rateRange = [100,500]  # range of random stimulation rate (Hz)
	durRange = [200,201]  # rate of random durations (ms)
	event_list = []  
	event1=[mistCell,mistStart,mistRate]  # first stim cell fixed
	event2=[mistCell,mistStart+mistDuration,0] # 
	event_list.append(event1) 
	event_list.append(event2) 
	for i in range(rnd.randrange(maxCells)):
		randCell = rnd.randrange(0,192,1)
		event = [randCell, mistStart, mistRate] #randrange(rateRange[0], rateRange[1])]
		event_list.append(event)
		event = [randCell, mistStart+mistDuration, 0]#randrange(durRange[0], durRange[1]), 0]
		event_list.append(event)
	
	events=vstack(event_list)
	
	outfilename="%s-in.txt"%(filestem) 
	savetxt(outfilename,events,fmt='%11.3f')
	print "Stimulus file %s saved." % outfilename

# save list of stimulated cells in file
def makeRepairStims (filestem, target, param1, param2):
	fileorigin = '%s-repair_target-%d_ptype-%d_pperc-%d-in.txt' % (filestem, target, param1, param2)
	filedest = "%s-in.txt" % (filestem)
	copyfile(fileorigin, filestem)
	print "Repair stimulus file %s copied to %s." % (fileorigin, filedest)

# Load params 
paramFile = 'stimdata/gen_86_cand_151'
with open('%s_params'% (paramFile)) as f:
	param = pickle.load(f) # read initial params from file
print param

# Set sim duration (in s)
simdur = 2

# number of different iseeds required
iseedstart = 1235
iseedstep = 17

# seed values
niseeds=50
wseedvals = 120456
iseedvals = arange(iseedstart, iseedstart+(niseeds*iseedstep), iseedstep) 

target_range = [1,3] #arange(4) # only four targets
wseed_range = arange(1)
inseed_range = arange(niseeds)

probing = 0 # 0 = Repair, 1 = single cell, 2 = multiple cell, 3 = special cases (multiple test seeds)

param1_range = [0,1] #kill cells or synapses
param2_range = [5,10] # kill 5% or 10%

if probing == 3:
	param3_range = [iseedvals[i] for i in arange(0,19)]#arange(0,191,1) # mistCell
	param4_range = [999]#[400] # mistStart
	param5_range = [1]#[200] # mistDuration
	param6_range = [1]#[250] # mistRate
	repair_range = [''];
elif probing:
	param3_range = arange(0,191)#arange(0,191,1) # mistCell
	param4_range = [100]#[400] # mistStart
	param5_range = [200]#[200] # mistDuration
	param6_range = [250]#[250] # mistRate
	repair_range = [''];
else:
	param3_range = [1000] #arange(0,191)#arange(0,191,1) # mistCell
	param4_range = [100]#[400] # mistStart
	param5_range = [200]#[200] # mistDuration
	param6_range = [250]#[250] # mistRate
	repair_range = ['both']; #['single', 'multi'];#, 

selectCellPerturbs = [[2,2],[30,9]] 
selectSynPerturbs = [[43,45],[8,34]] 

# save params
with open('stimdata/params', 'w') as f:
	pickle.dump([param, wseed_range, inseed_range], f)

# to save LFPs and EMGs"
calcErr = -1

# Loop over repair types (using single, multi or both probing seqs)
for irepair in repair_range:
	targetCounter = -1
	# Loop over target vals
	for target in target_range:
		targetCounter += 1
		# Loop over wiring seed...
		if target == 1:
			iseed = iseedvals[12]
		elif target == 3:
			iseed = iseedvals[22]
		mistparam1 = 'stimdata/target-%d_i-%d_w-120456_plastnq.nqs' % (target,iseed)
		for wseed in wseed_range:
			wseedval = wseedvals
			# Loop over param1 vals
			for param1 in param1_range:
				param2Counter = -1
				# Loop over param2 vals
				for param2 in param2_range:
					param2Counter += 1
					# Loop over param3 vals
					for param3 in param3_range:
						# Loop over param4 vals
						for param4 in param4_range:
							# Loop over param5 vals
							for param5 in param5_range:
								# Loop over param6 vals
								for param6 in param6_range:
									if (param5 > 0 and param6 > 0):# generate microstim
										if probing == 1: # probing data with single cell stim (filename includes mist details)
											outfilestem = 'stimdata/target-%d_ptype-%d_pperc-%d_cell-%d_start-%d_dur-%d_rate-%d' % \
											(target, param1, param2, param3, param4, param5, param6)
											makeProbingStims(outfilestem, param3,param4,param5,param6) # generate .txt file probing stim
										elif probing == 2: # probing data with multiple cell stim (filename includes mist details and 'multi')
											outfilestem = 'stimdata/target-%d_ptype-%d_pperc-%d_cell-%d_start-%d_dur-%d_rate-%d_multi' % \
											(target, param1, param2, param3, param4, param5, param6)
											makeMultProbingStims(outfilestem, param3,param4,param5,param6) # generate .txt file probing stim
										elif probing == 3: # special case: generating multiple test iseeds (no probing or repair)
											outfilestem = 'stimdata/target-%d_ptype-%d_pperc-%d_iseed-%d_start-%d_dur-%d_rate-%d' % \
											(target, param1, param2, param3, param4, param5, param6)
										elif probing == 0: # repair data (file name says 'repair')
											outfilestem = 'stimdata/target-%d_ptype-%d_pperc-%d_%s_repair' % (target, param1, param2, irepair)				   

									# generate perturbations (saves data to file so can be read by sim)
									# or copies over specific perturbation file previously selected
									if param1 == 0:
										iPerturb = selectCellPerturbs[param2Counter][targetCounter]
										filename = 'stimdata/target-%d_ptype-%d_pperc-%d_cell-%d_start-%d_dur-%d_rate-%d' % \
									(target, param1, param2, iPerturb, 0, 0, 0)
										copyPerturbCellsFile(filename, outfilestem)
									elif param1 == 1:
										iPerturb = selectSynPerturbs[param2Counter][targetCounter]
										filename = 'stimdata/target-%d_ptype-%d_pperc-%d_cell-%d_start-%d_dur-%d_rate-%d' % \
									(target, param1, param2, iPerturb, 0, 0, 0)
										copyPerturbSynsFile(filename, outfilestem)

									print filename

									# launch sim
									if probing == 3: 
										sys_str = './runsim_neurostim %s %.2f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %s' % \
										(outfilestem, param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], \
										 param[16], param[17], param[18], param[19], param[20], target, param3, wseedval, calcErr, mistparam1)
									else:
										sys_str = './runsim_neurostim %s %.2f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %s' % \
										(outfilestem, param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], \
										 param[16], param[17], param[18], param[19], param[20], target, iseed, wseedval, calcErr, mistparam1)
										
									print sys_str

									# Open a pipe to the qsub command.
									# output, input = popen2('qsub')
									
									# # Customize your options here
									# job_name = outfilestem
									# walltime = "00:30:00"
									# processors = "nodes=1:ppn=2"
									# command = sys_str

									# job_string = """#!/bin/bash
									# #PBS -N %s
									# #PBS -l walltime=%s
									# #PBS -q longerq
									# #PBS -l %s
									# #PBS -o %s.out
									# #PBS -e %s.err
									# cd $PBS_O_WORKDIR
									# echo $PBS_O_WORKDIR
									# %s""" % (job_name, walltime, processors, job_name, job_name, command)

									# # Send job_string to qsub
									# input.write(job_string)
									# input.close()
									os.system(sys_str)

									# Print your job and the response to the screen
									# print job_string
									# print output.read()

									time.sleep(0.1)
									

