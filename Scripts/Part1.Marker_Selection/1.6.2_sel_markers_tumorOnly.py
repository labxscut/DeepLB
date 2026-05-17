"""Select markers for the tumor-only method using difference criteria."""
import sys
import os
import re
import numpy as np
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests as mt
import math

def readData(fileName,type):
    data = []
    with open(fileName) as f:
        for line in f:
            line = re.sub('NA','nan',line)
            info = line.rstrip().split('\t')
            sampleType = info.pop(0)
            if sampleType==type:
                sampleProfile = list(map(float,info))
                data.append( sampleProfile )
    return np.array(data)

# This function applies to read betas file as well (alpha and beta are shape parameters of Beta-distribution)
def readAlphas(fileName):
        data = {}
        with open (fileName) as f:
                for line in f:
                        info = line.rstrip().split('\t')
                        sampleType = info.pop(0)
                        data[sampleType] = info
        return data


method = sys.argv[1]
tumor_type = sys.argv[2] # for example, 'lihc', 'luad', 'lusc', 'coad', 'read', 'brca'
file = sys.argv[3]
file_matched = sys.argv[4]
file_marker_def = sys.argv[5]
file_alphas = sys.argv[6]
file_betas = sys.argv[7]
out_file_prefix = sys.argv[8]

normalType = 'plasma_background'
tumorType = tumor_type + '_tumor'

dataPlasma = readData(file, normalType)
dataTumor = readData(file, tumorType)
print(len(dataTumor))
print(len(dataPlasma))
(nTumor, n_marker) = dataTumor.shape
(nPlasma, _) = dataPlasma.shape
# read start/end position of all markers from marker definition file
all_markers = {}
with open(file_marker_def) as f:
	next(f) # skip the first header line
	for line in f:
		(ind, chr, start, end) = line.rstrip().split('\t')
		ind = int(ind) # 1-base
		all_markers[ind] = {'index':ind, 'chr':chr, 'start':start, 'end':end}

# read alphas and betas file
alphas = readAlphas(file_alphas)
betas = readAlphas(file_betas)

if method.startswith('median-diff'):
	# We use this method
	params = method.split(',') # for example, median-diff,0.3
	print(method)
	diff_min = float(params[1])
	# step 1: modified from the code sel_markers_paired.py for CancerDetector_NAR paper
	selected_hyper  = list(range(n_marker))
	selected_hypo  = list(range(n_marker))
	
	# step 2: markers that the median difference between tumor and normal plasma >= cutoff
	NAN = float('nan')
	d_hyper = []
	for i in selected_hyper:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hyper.append( np.nanmedian(dataTumor[:,i]) - np.nanmedian(dataPlasma[:,i]) )
		else:
			d_hyper.append( NAN )
	final_selected_hyper = np.where( np.array(d_hyper) >= diff_min )
	final_selected_hyper = list( np.array(selected_hyper)[final_selected_hyper[0]] )
	d_hypo = []
	for i in selected_hypo:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hypo.append( np.nanmedian(dataPlasma[:,i]) - np.nanmedian(dataTumor[:,i]) )
		else:
			d_hypo.append( NAN )
	final_selected_hypo = np.where( np.array(d_hypo) >= diff_min )
	final_selected_hypo = list( np.array(selected_hypo)[final_selected_hypo[0]] )
	print('#final_hyper:', len(final_selected_hyper))
	print('#final_hypo:', len(final_selected_hypo))

	# step 3: output two formats
	final_selected_hyper = list( np.array(final_selected_hyper) +1 ) # convert index from 0-base to 1-base
	final_selected_hypo = list( np.array(final_selected_hypo) +1 )
	selected = [(i, 'hyper') for i in final_selected_hyper] + [(i, 'hypo') for i in final_selected_hypo]
	selected.sort(key=lambda tup: tup[0])
	normal_type = normalType
	# format 1: tumor_type list_of_selected_marker_indexes
	out_file = out_file_prefix + '.sk'
	with open(out_file,'w') as f:
		f.write(tumor_type)
		for (ind, type) in selected:
			if alphas[tumor_type][ind-1]!='NA' and betas[tumor_type][ind-1]!='NA' \
				and alphas[normal_type][ind-1]!='NA' and betas[normal_type][ind-1]!='NA': 
				f.write('\t'+str(ind))
		f.write('\n')
	f.close()
	# format 2: each line is a selected marker with marker index, genomic-region, type, alphas, betas
	# output full details of selected markers into a file
	out_file = out_file_prefix
	with open(out_file, 'w') as f:
		f.write('marker_index\tchr\tstart\tend\tmarker_type\t'+tumor_type+'\t'+normal_type+'\n')
		for (ind, type) in selected:
			if alphas[tumor_type][ind-1]!='NA' and betas[tumor_type][ind-1]!='NA' \
				and alphas[normal_type][ind-1]!='NA' and betas[normal_type][ind-1]!='NA': 
				f.write(str(ind) + '\t' + all_markers[ind]['chr'] + '\t' + \
					all_markers[ind]['start'] + '\t' + all_markers[ind]['end'] + '\t' + \
					type + '\t' + alphas[tumor_type][ind-1] + ':' + betas[tumor_type][ind-1] + '\t' + \
					alphas[normal_type][ind-1] + ':' + betas[normal_type][ind-1] + '\n')
			else:
				sys.stderr.write("marker "+str(ind)+"has NA in alphas and betas\n")

	f.close()
	# format 3: each line is a selected marker with marker index and type
	# We omit this format. it is only for checking
	# out_file = out_file_prefix
	# with open(out_file,'w') as f:
		# for (ind, type) in selected:
			# if alphas[tumor_type][ind-1]!='NA' and betas[tumor_type][ind-1]!='NA' \
				# and alphas[normal_type][ind-1]!='NA' and betas[normal_type][ind-1]!='NA': 
				# f.write(str(ind)+"\t"+type+"\n")
	# f.close()


elif method.startswith('mean-diff'):
	print("mean difference")
	params = method.split(',') # for example, mean-diff,0.3
	print(method)
	diff_min = float(params[1])
	# step 1: modified from the code sel_markers_paired.py for CancerDetector_NAR paper
	selected_hyper  = list(range(n_marker))
	selected_hypo  = list(range(n_marker))

	# step 2: markers that the mean difference between tumor and normal plasma >= cutoff
	NAN = float('nan')
	d_hyper = []
	for i in selected_hyper:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hyper.append( np.nanmean(dataTumor[:,i]) - np.nanmean(dataPlasma[:,i]) )
		else:
			d_hyper.append( NAN )
	final_selected_hyper = np.where( np.array(d_hyper) >= diff_min )
	final_selected_hyper = list( np.array(selected_hyper)[final_selected_hyper[0]] )
	d_hypo = []
	for i in selected_hypo:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hypo.append( np.nanmean(dataPlasma[:,i]) - np.nanmean(dataTumor[:,i]) )
		else:
			d_hypo.append( NAN )
	final_selected_hypo = np.where( np.array(d_hypo) >= diff_min )
	final_selected_hypo = list( np.array(selected_hypo)[final_selected_hypo[0]] )
	print('#final_hyper:', len(final_selected_hyper))
	print('#final_hypo:', len(final_selected_hypo))

	# step 3: output two formats
	final_selected_hyper = list( np.array(final_selected_hyper) +1 ) # convert index from 0-base to 1-base
	final_selected_hypo = list( np.array(final_selected_hypo) +1 )
	selected = [(i, 'hyper') for i in final_selected_hyper] + [(i, 'hypo') for i in final_selected_hypo]
	selected.sort(key=lambda tup: tup[0])
	normal_type = normalType
	# format 1: tumor_type list_of_selected_marker_indexes
	out_file = out_file_prefix + '.sk'
	with open(out_file,'w') as f:
		f.write(tumor_type)
		for (ind, type) in selected:
			if alphas[tumor_type][ind-1]!='NA' and betas[tumor_type][ind-1]!='NA' \
				and alphas[normal_type][ind-1]!='NA' and betas[normal_type][ind-1]!='NA':
				f.write('\t'+str(ind))
		f.write('\n')
	f.close()
	# format 2: each line is a selected marker with marker index, genomic-region, type, alphas, betas
	# output full details of selected markers into a file
	out_file = out_file_prefix
	with open(out_file, 'w') as f:
		f.write('marker_index\tchr\tstart\tend\tmarker_type\t'+tumor_type+'\t'+normal_type+'\n')
		for (ind, type) in selected:
			if alphas[tumor_type][ind-1]!='NA' and betas[tumor_type][ind-1]!='NA' \
				and alphas[normal_type][ind-1]!='NA' and betas[normal_type][ind-1]!='NA':
				f.write(str(ind) + '\t' + all_markers[ind]['chr'] + '\t' + \
					all_markers[ind]['start'] + '\t' + all_markers[ind]['end'] + '\t' + \
					type + '\t' + alphas[tumor_type][ind-1] + ':' + betas[tumor_type][ind-1] + '\t' + \
					alphas[normal_type][ind-1] + ':' + betas[normal_type][ind-1] + '\n')
			else:
				sys.stderr.write("marker "+str(ind)+"has NA in alphas and betas\n")

	f.close()

