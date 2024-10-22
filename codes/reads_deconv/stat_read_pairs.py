# stat_read_pairs.py
#
# Output (1) average overlap rate, (2) overlap rate distribution
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#
# An example of the file of 9 columns. Additional columns may be added, such as 'read_name' or 'RG'.
#
# marker1_index	marker2_index	chr	start_pos	end_pos	strand	num_CpG_sites	CpG_sites	sites_methylation_status
# 0	2	chr1	566950	567024	-	1	567014	1
# 0	2	chr1	566952	567026	-	1	567014	1
# 3	-	chr1	569458	569530	-	1	569507	0
# 3	-	chr1	569537	569611	-	4	569546,569550,569558,569603	1111
# 3	-	chr1	569538	569612	-	4	569546,569550,569558,569603	0011
# ...
#
def stat_single_end_reads( reads_file ):
	r_size_list = []
	if reads_file=='stdin':
		f = sys.stdin
	else:
		f = open(reads_file)
	header_items = f.readline().rstrip().split('\t') # the first line is header
	for line in f:
		items =  line.rstrip().split('\t')
		if len(items)<9:
			sys.stderr.write('Error: only <9 columns!\nExit.\n')
			sys.exit(0)
		r_size = int(items[4]) - int(items[3])
		r_size_list.append( r_size )
	if reads_file!='stdin':
		f.close()
	return r_size_list

#
# An example of the file: the 10th column must be 'endpos_of_read1,startpos_of_read2', otherwise we cannot calculate fragment size and overlap size
#
# marker1_index	marker2_index	chr	start_pos	end_pos	strand	num_CpG_sites	CpG_sites	sites_methylation_status	endpos_of_read1,startpos_of_read2	read_name	RG
# 2	-	chr1	566995	567158	+-	5	567014,567062,567090,567113,567122	11111	567115,567049	ST-E00114:459:HLLN7CCXX:4:1109:20841:50850_1:N:0:CTTGTAAT	N2L
# 2	-	chr1	567034	567181	+-	5	567062,567090,567113,567122,567167	11111	567154,567072	ST-E00114:459:HLLN7CCXX:4:2214:2077:69713_1:N:0:NTTGTAAT	N2L
# 0	4	chr1	713756	713876	-	2	713823,713825	00	NA	ST-E00114:459:HLLN7CCXX:4:2207:4675:47527_1:N:0:CTTGTAAT	N2L
# ...
#
# We do not count those single forward (or reverse) reads, r1 (or r2).
#
def stat_paired_end_reads( reads_file ):
	fragment_size_list = []
	r1_size_list = []
	r2_size_list = []
	overlap_size_list = []
	overlap_rate_list = [] # overlap_size/fragment_size
	if reads_file=='stdin':
		f = sys.stdin
	else:
		f = open(reads_file)
	header_items = f.readline().rstrip().split('\t') # the first line is header
	for line in f:
		items =  line.rstrip().split('\t')
		if len(items)<=9:
			sys.stderr.write('Error: only <=9 columns, no r1_end_position and r2_start_position (in the 10th column) are provided; so fragment size cannot be calculated!\nExit.\n')
			sys.exit(0)
		if header_items[9]!='endpos_of_read1,startpos_of_read2':
			sys.stderr.write('Error: the 10th column is not for r1_end_position and r2_start_position; so fragment size cannot be calculated!\nExit.\n')
			sys.exit(0)
		r1_start = int(items[3])
		r2_end = int(items[4])
		strand = items[5]
		sites_int = np.array( map(int, items[7].split(',')) )
		additional_pos = items[9]
		# (m1_idx, m2_idx, chr, r1_start, r2_end, strand, n_sites, sites, methy_status, additional_pos, read_name, rg) = line.rstrip().split('\t')
		if strand=='+' or strand=='-':
			# forward or reverse single-end read
			continue
		else:
			# read pair - fragment
			(r1_end, r2_start) = additional_pos.split(',')
			r1_end = int(r1_end)
			r2_start = int(r2_start)
			r1_size = r1_end - r1_start
			r2_size = r2_end - r2_start
			fragment_size = r2_end - r1_start
			overlap_size = r1_end - r2_start 
			if (overlap_size<0):
				overlap_size = 0
			r1_size_list.append( r1_size )
			r2_size_list.append( r2_size )
			fragment_size_list.append( fragment_size )
			overlap_size_list.append( overlap_size )
			overlap_rate_list.append( float(overlap_size)/fragment_size )
	if reads_file!='stdin':
		f.close()
	return (r1_size_list, r2_size_list, fragment_size_list, overlap_size_list, overlap_rate_list)

def plt_hist(fig, ax, vals, bins, title, xlabel, ylabel, median_=None):
	bins = np.array(bins)
	width = bins[1] - bins[0]
	center = (bins[:-1] + bins[1:]) / 2
	ax.set_title(title, fontsize=18)
	ax.set_xlabel(xlabel, fontsize=14)
	ax.set_ylabel(ylabel, fontsize=14)
	ax.bar(center, vals, align='center', width=width)
	if median_:
		ax.axvline(x=median_, linewidth=1, color='r')
		# ax.annotate('median size=%.6g'%median_, xy=(median_, 0.75*np.max(vals)), textcoords='offset points', arrowprops=dict(arrowstyle="->", facecolor='red'))

readsFile = sys.argv[1]
readType = sys.argv[2] # 'single-end' or 'paired-end'
figureFile = None
if len(sys.argv)==4:
	figureFile = sys.argv[3] 
num_bins = 20

if readType=='single-end':
	r_size_list = stat_single_end_reads( readsFile )
	r_size_list = np.array(r_size_list)
	print '----------- summary -----------------'
	read_num = '#read: %d'%len(r_size_list)
	r_size = 'r_size: %.6g +- %.6g'%(np.mean(r_size_list), np.std(r_size_list))
	print read_num
	print r_size
	print ''

	mpl.rc('xtick', labelsize=14)
	mpl.rc('ytick', labelsize=14)
	fig, axarr = plt.subplots(1, 2, figsize=(12,6))
	fig.tight_layout()
	plt.subplots_adjust(left=.08, top=.95, bottom=.05, wspace=.15,  hspace=.3)
	print '------------ histogram --------------'
	bins = range(max(r_size_list)+1)
	h, _ = np.histogram(r_size_list, bins)
	plt_hist(fig, axarr[0], h, bins, 'read length distribution', 'read length', 'frequency', np.median(r_size_list))
	print 'read_size_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''

	font = {'family': 'serif', \
		'color':  'darkred', \
		'weight': 'normal', \
		'size': 14, \
        }
	axarr[1].axis('off')
	axarr[1].text(0.05, 0.65,  read_num,  fontdict=font)
	axarr[1].text(0.05, 0.45,  r_size.replace('+-','$\pm$'), fontdict=font)
	if figureFile:
		plt.savefig(figureFile, dpi=300)
	else:
		plt.show()

elif readType=='paired-end':
	(r1_size_list, r2_size_list, fragment_size_list, overlap_size_list, overlap_rate_list) = stat_paired_end_reads( readsFile )
	r1_size_list = np.array(r1_size_list)
	r2_size_list = np.array(r2_size_list)
	fragment_size_list = np.array(fragment_size_list)
	overlap_size_list = np.array(overlap_size_list)
	overlap_rate_list = np.array(overlap_rate_list)
	print '----------- summary -----------------'
	fragment_num = '#fragment (of markers): %d'%len(fragment_size_list)
	fragment_size = 'avg_fragment_size: %.6g +- %.6g'%(np.mean(fragment_size_list), np.std(fragment_size_list))
	fragment_size_median = 'median_fragment_size: %.6g'%(np.median(fragment_size_list))
	r1_size = 'r1_size: %.6g +- %.6g'%(np.mean(r1_size_list), np.std(r1_size_list))
	r2_size = 'r2_size: %.6g +- %.6g'%(np.mean(r2_size_list), np.std(r2_size_list))
	overlap_size = 'overlap_size: %.6g +- %.6g'%(np.mean(overlap_size_list), np.std(overlap_size_list))
	overlap_rate = 'overlap_rate: %.6g +- %.6g'%(np.mean(overlap_rate_list), np.std(overlap_rate_list))
	print fragment_num
	print fragment_size
	print fragment_size_median
	print r1_size
	print r2_size
	print overlap_size
	print overlap_rate
	print ''

	mpl.rc('xtick', labelsize=14)
	mpl.rc('ytick', labelsize=14)
	fig, axarr = plt.subplots(3, 2, figsize=(12,12))
	fig.tight_layout()
	plt.subplots_adjust(left=.08, top=.95, bottom=.05, wspace=.15,  hspace=.3)
	print '------------ histogram --------------'
	bins = range(max(fragment_size_list)+1)
	h, _ = np.histogram(fragment_size_list, bins)
	plt_hist(fig, axarr[0,0], h, bins, 'fragment size distribution', 'fragment size', 'frequency', np.median(fragment_size_list))
	print 'fragment_size_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''
	bins = range(max(r1_size_list)+1)
	h, _ = np.histogram(r1_size_list, bins)
	plt_hist(fig, axarr[0,1], h, bins, 'forward-strand read size distribution', 'read length', 'frequency')
	print 'r1_size_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''
	bins = range(max(r2_size_list)+1)
	h, _ = np.histogram(r2_size_list, bins)
	plt_hist(fig, axarr[1,0], h, bins, 'reverse-strand read length distribution', 'read length', 'frequency')
	print 'r2_size_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''
	bins = range(max(overlap_size_list)+1)
	h, _ = np.histogram(overlap_size_list, bins)
	plt_hist(fig, axarr[1,1], h, bins, 'overlap size distribution', 'overlap size', 'frequency')
	print 'overlap_size_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''
	bins = np.array(range(num_bins+1))/float(num_bins)
	h, _ = np.histogram(overlap_rate_list, bins)
	plt_hist(fig, axarr[2,0], h, bins, 'overlap rate distribution', 'overlap rate', 'frequency')
	print 'overlap_rate_hist: %s'%', '.join(map(str, h))
	print 'bins: %s'%', '.join(map(str, bins))
	print ''

	font = {'family': 'serif', \
		'color':  'darkred', \
		'weight': 'normal', \
		'size': 14, \
        }
	axarr[2,1].axis('off')
	axarr[2,1].text(0.1, 0.95,  fragment_num,  fontdict=font)
	axarr[2,1].text(0.1, 0.8, fragment_size.replace('+-','$\pm$'), fontdict=font)
	axarr[2,1].text(0.1, 0.65, fragment_size_median.replace('+-','$\pm$'), fontdict=font)
	axarr[2,1].text(0.1, 0.5,  r1_size.replace('+-','$\pm$'),       fontdict=font)
	axarr[2,1].text(0.1, 0.35, r2_size.replace('+-','$\pm$'),       fontdict=font)
	axarr[2,1].text(0.1, 0.2,  overlap_size.replace('+-','$\pm$'),  fontdict=font)
	axarr[2,1].text(0.1, 0.05, overlap_rate.replace('+-','$\pm$'),  fontdict=font)
	if figureFile:
		plt.savefig(figureFile, dpi=300)
	else:
		plt.show()
else:
	sys.stderr.write( 'Error: the 2nd argument readType (%s) is wrong!\nExit.\n'%readType )

