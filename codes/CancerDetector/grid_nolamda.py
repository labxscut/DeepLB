import sys
import os
import configparser
import re
import operator
import numpy as np
import math
from copy import deepcopy

def ismember(A, B): # A and B are numpy.array
	return [np.sum(B==a) for a in A]

def loadReadProbilities(fileName):
	markerIdReads = []
	p = []
	with open(fileName) as f:
		next(f) # skip the first header line
		for line in f:
			items = line.rstrip().split('\t')
			markerIdReads.append( int(items[0]) )
			p.append( np.array( list(map(float,items[1:])) ) )
	f.close()
	if not p: # p is an empty list
		return (None, None)
	else:
		return (np.array(markerIdReads), np.matrix(p))

def grid_search( p, nGrid=1000):
	(N, M) = p.shape # N=#reads, M=2 (#class)
	theta = np.linspace(0, 1, num=nGrid+1)
	theta = np.matrix( np.vstack( (theta, 1-theta) ) )
	obj = np.sum( np.log( p*theta ), axis=0 ) # column sums
	indMax = np.nanargmax( obj )
	theta_best = [float(theta[:,indMax][0]), float(theta[:,indMax][1]) ]
	return ( theta_best, float(obj[:,indMax]) )

def grid_bins( markerIdReads, p ):
	theta_bins = []
	max_obj_bins = []
	marker_ids = np.unique( markerIdReads )
	[N, M] = p.shape
	for i in range(len(marker_ids)):
		marker_id = marker_ids[i]
		p_marker = p[markerIdReads==marker_id, :]
		(theta_best, obj_max) = grid_search( p_marker )
		theta_bins.append( theta_best )
		max_obj_bins.append( obj_max )
	return (np.array(theta_bins), np.array(max_obj_bins), marker_ids)

readProbFile = sys.argv[1]
outFilePrefix = sys.argv[2]

print('Load read prob file ...')
sys.stdout.flush()
markerIdReads, p = loadReadProbilities( readProbFile )
if p is None:
	sys.stderr.write('Warning: read prob file (%s) is empty!\nExit.\n'%readProbFile)
	sys.exit()

# Remove the reads whose probabilities to each class is zero, especially when #class=2, they will lead to NAN in computation and so need to be removed before computation
goodReadsIdx = np.where(np.sum(p, axis=1) !=0)[0]
p = p[goodReadsIdx, :]
markerIdReads = markerIdReads[goodReadsIdx]
if markerIdReads.size==0:
	sys.stderr.write('Warning: read prob file (%s) is empty, after removing reads with class-specific probabilities=0!\nExit.\n'%readProbFile)
	sys.exit()

# step 1. infer tumor burden of each marker
print('step 1. infer tumor burden of each marker ...')
sys.stdout.flush()
theta_bins_both, max_obj_bins, marker_ids = grid_bins( markerIdReads, p )
theta_bins = theta_bins_both[:,1] # the 2nd column is theta of normal samples
theta_std = np.std(theta_bins, ddof=1)
thetaBinsFile = outFilePrefix + '.theta_marker'
with open(thetaBinsFile, 'w') as f:
	for i in range(len(list(marker_ids))):
		f.write('%d\t%6.10f\t%6.10f\t%e\n'%(marker_ids[i], theta_bins_both[i,0], theta_bins_both[i,1], max_obj_bins[i]))
f.close()


# step 2. infer tumor burden over all markers
# after 20 round, the tumor burden should converge
print('step 2. infer tumor burden over all markers ...')
sys.stdout.flush()
x_list = [1000] #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
nRound = 20
logFile = outFilePrefix + '.theta.log'
outputFile = outFilePrefix + '.theta'
with open(logFile, 'w') as f1, open(outputFile, 'w') as f2:
	for x in x_list:
		f1.write('++++++ x = %g ++++++\n'%x)
		f1.flush()
		print('++++++ x = %g ++++++'%x)
		sys.stdout.flush()
		p_goodmarkers = p
		good_marker_ids = marker_ids

		n_good_marker_prev = 0
		for round in range(nRound):
			print('round %d. run grid search for all markers'%round)
			sys.stdout.flush()
			thetaVectorUpdateBest, max_objective_value_best = grid_search( p_goodmarkers )
			theta_cutoff_biomarker = thetaVectorUpdateBest[1] - theta_std * x # x is a parameter to adjust the normal cfDNA fraction. We suppose the 2nd element of 'thetaVectorUpdateBest' is the fraction of normal cfDNA
			f1.write('%d\t%6.10f\t%6.10f\t%d\t%6.10f\n'%(round,thetaVectorUpdateBest[0],thetaVectorUpdateBest[1],len(good_marker_ids),theta_cutoff_biomarker))
			f1.flush()
			print('round %d. Theta_cutoff_marker %f'%(round,theta_cutoff_biomarker))
			sys.stdout.flush()
			good_marker_ids = marker_ids[theta_bins >= theta_cutoff_biomarker]
			if len(good_marker_ids) != n_good_marker_prev:
				n_good_marker_prev = len(good_marker_ids)
			else:
				break
			good_reads_indexes = np.where( ismember(markerIdReads, good_marker_ids) )[0];
			p_goodmarkers = p[good_reads_indexes, :]
			print('round %d. %%good_markers: %.4f = %d / %d'%(round,len(good_marker_ids)/float(len(marker_ids)), len(good_marker_ids), len(marker_ids)))
			print('round %d. %%good_reads: %.4f = %d / %d'%(round,sum(good_reads_indexes)/float(len(markerIdReads)),sum(good_reads_indexes),len(markerIdReads)))
			sys.stdout.flush()

		f2.write('%g\t%6.10f\t%6.10f\t%d\t%.4f\t%d\t'%(x, thetaVectorUpdateBest[0], thetaVectorUpdateBest[1], len(good_marker_ids), len(good_marker_ids)/float(len(marker_ids)), len(marker_ids)))
		f2.flush()
		for i in range(len(good_marker_ids)-1):
			f2.write('%d,'%good_marker_ids[i])
		f2.write('%d\n'%good_marker_ids[-1])
		f2.flush()
f1.close()
f2.close()

# with open('ttt1','w') as f:
	# [m,n] = p.shape
	# for i in range(m):
		# f.write('%g\t%g\n'%(p[i,0],p[i,1]))
# f.close()
