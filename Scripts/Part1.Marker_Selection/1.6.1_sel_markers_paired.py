"""Select markers for the paired method using differential criteria."""
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
file = sys.argv[3] #train_1
file_matched = sys.argv[4]
file_marker_def = sys.argv[5]
file_alphas = sys.argv[6]
file_betas = sys.argv[7]
out_file_prefix = sys.argv[8]

normalType = 'plasma_background'
tumorType = tumor_type + '_tumor'
matchedNormaltumorType = tumor_type + '_normal'

dataPlasma = readData(file, normalType)
dataTumor = readData(file, tumorType)
dataMatched = readData(file_matched, matchedNormaltumorType)
print(dataMatched.shape)
(n_matched, _) = dataMatched.shape
print(n_matched)
(nPlasma, _) = dataPlasma.shape
print(nPlasma)


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

if method.startswith('paired-t-test'):
# paired t-test for each marker, and adjust p-value by multiple testing
# however, this does not work, because t-test considers variance
	params = method.split(',') # for example, paired-t-test,0.000005,bonferroni or paired-t-test,0.000005,fdr_bh
	p_cutoff = float(params[1])
	mul_test = params[2]
	tScores = []
	pValues = []
	for i in range(n_marker):
		(t, p) = stats.ttest_rel(dataTumor[:,i], dataMatched[:,i])
		tScores.append(t)
		pValues.append(p)
	tScores = np.array(tScores)
	pValues = np.array(pValues)
	(reject, qValues, _, _) = mt(pValues, alpha=p_cutoff, method=mul_test)
	print(qValues)
	print(reject)
	print(np.sum(reject))
	selected = np.where(reject)
	print("t(%d):"%int(selected[0][0]), dataTumor[:,selected[0][0]])
	print("n(%d):"%int(selected[0][0]), dataMatched[:,selected[0][0]])
	print("p(%d):"%int(selected[0][0]), float(pValues[selected[0][0]]))
	print("q(%d):"%int(selected[0][0]), float(qValues[selected[0][0]]))
elif method.startswith('freq-diff'):
	# We use this method
	params = method.split(',') # for example, freq-diff,0.3
	print(method)
	# step 1: markers that frequently differentiate tumor and matched normal tissue in >= half of the population
	freq_min = n_matched / 2
	diff_min = float(params[1])
	d_hyper = ( dataTumor - dataMatched ) >= diff_min
	d_hypo  = ( dataMatched - dataTumor ) >= diff_min
	#a matrix whose shape as same as dataTumor(dataMatched) contain False or True (0/1)
	freq_hyper = np.sum(d_hyper, axis=0)  # sum columns
	freq_hypo  = np.sum(d_hypo, axis=0)  # sum columns
	selected_hyper = np.where( freq_hyper >= freq_min )
	#filter the cluster which different at least over half samples
	selected_hyper = list(selected_hyper[0])
	selected_hypo  = np.where( freq_hypo >= freq_min )
	selected_hypo  = list(selected_hypo[0])
	print('#hyper:', len(selected_hyper))
	print('#hypo:', len(selected_hypo))

	# test print
	# print "t(%d):"%int(selected_hyper[0]), dataTumor[:,selected_hyper[0]]
	# print "n(%d):"%int(selected_hyper[0]), dataMatched[:,selected_hyper[0]]
	# print "d(%d):"%int(selected_hyper[0]), d_hyper[:, selected_hyper[0]]
	# print "q(%d):"%int(selected_hyper[0]), freq_hyper[selected_hyper[0]]
	# print "t(%d):"%int(selected_hypo[0]), dataTumor[:,selected_hypo[0]]
	# print "n(%d):"%int(selected_hypo[0]), dataMatched[:,selected_hypo[0]]
	# print "d(%d):"%int(selected_hypo[0]), d_hypo[:, selected_hypo[0]]
	# print "q(%d):"%int(selected_hypo[0]), freq_hypo[selected_hypo[0]]
	# end of test print

	# step 2: markers that the median difference between tumor and normal plasma >= cutoff
	NAN = float('nan')
	diff_min2 = diff_min
	nTumor = n_matched
	d_hyper = []
	for i in selected_hyper:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hyper.append( np.nanmedian(dataTumor[:,i]) - np.nanmedian(dataPlasma[:,i]) )
		else:
			d_hyper.append( NAN )
	final_selected_hyper = np.where( np.array(d_hyper) >= diff_min2 )
	final_selected_hyper = list( np.array(selected_hyper)[final_selected_hyper[0]] )
	d_hypo = []
	for i in selected_hypo:
		if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
			d_hypo.append( np.nanmedian(dataPlasma[:,i]) - np.nanmedian(dataTumor[:,i]) )
		else:
			d_hypo.append( NAN )
	final_selected_hypo = np.where( np.array(d_hypo) >= diff_min2 )
	final_selected_hypo = list( np.array(selected_hypo)[final_selected_hypo[0]] )
	print('#final_hyper:', len(final_selected_hyper))
	print('#final_hypo:', len(final_selected_hypo))
	
	# step 3: output two formats
	final_selected_hyper = list( np.array(final_selected_hyper) +1 ) # convert index from 0-base to 1-base
	final_selected_hypo = list( np.array(final_selected_hypo) +1 )
	selected = [(i, 'hyper') for i in final_selected_hyper] + [(i, 'hypo') for i in final_selected_hypo]
	selected.sort(key=lambda tup: tup[0])
	normal_type=normalType
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
elif method.startswith('fold-change'):
        # We first use freq-diff > 0.1, then fold-change > cutoff
        params = method.split(',') # for example, fold-change,10,0.1
        print(method)
        # step 1: markers that frequently differentiate tumor and matched normal tissue in >= half of the population
        freq_min = n_matched / 2
        fc_min = float(params[1])
        diff_min = float(params[2])
        d_hyper = ( dataTumor - dataMatched ) >= diff_min
        d_hypo  = ( dataMatched - dataTumor ) >= diff_min
        freq_hyper = np.sum(d_hyper, axis=0)  # sum columns
        freq_hypo  = np.sum(d_hypo, axis=0)  # sum columns
        selected_hyper = np.where( freq_hyper >= freq_min )
        selected_hyper = list(selected_hyper[0])
        selected_hypo  = np.where( freq_hypo >= freq_min )
        selected_hypo  = list(selected_hypo[0])
        print('#hyper:', len(selected_hyper))
        print('#hypo:', len(selected_hypo))
      
        # test print
        # print "t(%d):"%int(selected_hyper[0]), dataTumor[:,selected_hyper[0]]
        # print "n(%d):"%int(selected_hyper[0]), dataMatched[:,selected_hyper[0]]
        # print "d(%d):"%int(selected_hyper[0]), d_hyper[:, selected_hyper[0]]
        # print "q(%d):"%int(selected_hyper[0]), freq_hyper[selected_hyper[0]]
        # print "t(%d):"%int(selected_hypo[0]), dataTumor[:,selected_hypo[0]]
        # print "n(%d):"%int(selected_hypo[0]), dataMatched[:,selected_hypo[0]]
        # print "d(%d):"%int(selected_hypo[0]), d_hypo[:, selected_hypo[0]]
        # print "q(%d):"%int(selected_hypo[0]), freq_hypo[selected_hypo[0]]
        # end of test print

        # step 2: markers with median difference btw tumor and normal plasma in >= diff_min, exclude the NA's
        NAN = float('nan')
        diff_min2 = diff_min
        nTumor = n_matched
        d_hyper = []
        for i in selected_hyper:
                if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
                        d_hyper.append( np.nanmedian(dataTumor[:,i]) - np.nanmedian(dataPlasma[:,i]) )
                else:
                        d_hyper.append( NAN )
        final_selected_hyper = np.where( np.array(d_hyper) >= diff_min2 )
        final_selected_hyper = list( np.array(selected_hyper)[final_selected_hyper[0]] )
        d_hypo = []
        for i in selected_hypo:
                if ( np.count_nonzero(~np.isnan(dataTumor[:,i])) >= (nTumor/2) ) and ( np.count_nonzero(~np.isnan(dataPlasma[:,i])) >= (nPlasma/2) ):
                        d_hypo.append( np.nanmedian(dataPlasma[:,i]) - np.nanmedian(dataTumor[:,i]) )
                else:
                        d_hypo.append( NAN )
        final_selected_hypo = np.where( np.array(d_hypo) >= diff_min2 )
        final_selected_hypo = list( np.array(selected_hypo)[final_selected_hypo[0]] )
        print('#final_hyper:', len(final_selected_hyper))
        print('#final_hypo:', len(final_selected_hypo))


	# step 3: filter the markers by folder change
        fc_hyper = []
        normal_type=normalType
        for i in final_selected_hyper:
                mr_tumor = float(alphas[tumor_type][i]) / (float(alphas[tumor_type][i])+float(betas[tumor_type][i]))
                mr_normal = float(alphas[normal_type][i]) / (float(alphas[normal_type][i])+float(betas[normal_type][i]))
                fc = mr_tumor / mr_normal
                fc_hyper.append(fc)
        final_final_selected_hyper = np.where( np.array(fc_hyper) >= fc_min )
        final_final_selected_hyper = list( np.array(final_selected_hyper)[final_final_selected_hyper[0]] )
        fc_hypo = []
        for i in final_selected_hypo:
                mr_tumor = float(alphas[tumor_type][i]) / (float(alphas[tumor_type][i])+float(betas[tumor_type][i]))
                mr_normal = float(alphas[normal_type][i]) / (float(alphas[normal_type][i])+float(betas[normal_type][i]))
                fc = (1-mr_tumor) / (1-mr_normal)
                fc_hypo.append(fc)
        final_final_selected_hypo = np.where( np.array(fc_hypo) >= fc_min )
        final_final_selected_hypo = list( np.array(final_selected_hypo)[final_final_selected_hypo[0]] )

        print('#final_final_hyper:', len(final_final_selected_hyper))
        print('#final_final_hypo:', len(final_final_selected_hypo))


        # step 4: output two formats
        final_selected_hyper = list( np.array(final_final_selected_hyper) +1 ) # convert index from 0-base to 1-base
        final_selected_hypo = list( np.array(final_final_selected_hypo) +1 )
        selected = [(i, 'hyper') for i in final_selected_hyper] + [(i, 'hypo') for i in final_selected_hypo]
        #selected = [(i, 'hypo') for i in final_selected_hypo]
        selected.sort(key=lambda tup: tup[0])
        normal_type=normalType
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
	
