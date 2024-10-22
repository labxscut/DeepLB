# trim_read_pairs.py
#
# trim the ending bases of forward read (r1) and the beginning bases of reverse read (r2) for each fragment, so that each read's length is exactly the same as the parameter provided in the command argument (length_kept_per_read).
# print the trimmed reads to stdout
#
import sys
import numpy as np

#
# An example of the file for paired-end reads
# marker1_index	marker2_index	chr	start_pos	end_pos	strand	num_CpG_sites	CpG_sites	sites_methylation_status	endpos_of_read1,startpos_of_read2	read_name	RG
# 2	-	chr1	566995	567158	+-	5	567014,567062,567090,567113,567122	11111	567115,567049	ST-E00114:459:HLLN7CCXX:4:1109:20841:50850_1:N:0:CTTGTAAT	N2L
# 2	-	chr1	567034	567181	+-	5	567062,567090,567113,567122,567167	11111	567154,567072	ST-E00114:459:HLLN7CCXX:4:2214:2077:69713_1:N:0:NTTGTAAT	N2L
# 0	4	chr1	713756	713876	-	2	713823,713825	00	NA	ST-E00114:459:HLLN7CCXX:4:2207:4675:47527_1:N:0:CTTGTAAT	N2L
# ...
#
def trim_paired_reads(reads_file, length_kept_per_read):
	methy_bins = {}
	if reads_file=='stdin':
		f = sys.stdin
	else:
		f = open(reads_file)
	header_line = f.readline().rstrip() # the first line is header
	sys.stdout.write( header_line + '\n' )
	columns_names = header_line.rstrip().split('\t')
	n_columns = len(columns_names)
	n_cytosines = 0
	n_cytosines_cut = 0
	n_bases_covered = 0
	n_bases_covered_cut = 0
	marker_prev = -1
	if n_columns!=11 and n_columns!=12:
		sys.stderr.write('Error: #columns should be 11 or 12!\nExit.\n')
		sys.exit(0)
	if columns_names[9] != 'endpos_of_read1,startpos_of_read2':
		sys.stderr.write('Error: the 10th columns is not endpos_of_read1,startpos_of_read2 for paired reads.\nExit.\n')
		sys.exit(0)
	for line in f:
		if n_columns==11:
			(m1_idx, m2_idx, chr, r1_start, r2_end, strand, n_sites, sites, methy_status, additional_pos, read_name) = line.replace('\n','').split('\t')
		if n_columns==12:
			(m1_idx, m2_idx, chr, r1_start, r2_end, strand, n_sites, sites, methy_status, additional_pos, read_name, rg) = line.replace('\n','').split('\t')
		marker_index = None
		if m2_idx=='-' or m2_idx=='0':
			marker_index = int(m1_idx)
		if m1_idx=='-' or m1_idx=='0':
			marker_index = int(m2_idx)
		if not marker_index:
			sys.stderr.write('Error: no marker is assined to the read %s\nExit.\n'%sread_name)
			sys.exit(0)
		if marker_index!=marker_prev:
			# this is the first read in the new marker
			marker_prev = marker_index
			methy_bins[marker_index] = {'methylation_rate':0, 'total_count':0, 'methylation_count':0, 'unmethylation_count':0, '#reads':0, '#fragments':0}
		sites_int = np.array( map(int,sites.split(',')) )
		n_cytosines += len(sites_int)
		if strand=='+':
			# forward single-end read
			r1_start = int(r1_start)
			r1_end = int(r2_end)
			r1_cut_end   = min( r1_end, r1_start+length_kept_per_read )
			good_sites_idx = np.where(sites_int<r1_cut_end)[0]
			good_methy_status = [methy_status[i] for i in good_sites_idx]
			good_sites = [sites_int[i] for i in good_sites_idx]
			n_sites_cut = len(good_sites)
			if n_sites_cut>0:
				if n_columns==11:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r1_start, r1_cut_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), 'NA', read_name) )
				if n_columns==12:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r1_start, r1_cut_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), 'NA', read_name, rg) )
			n_cytosines_cut += len(sites_int) - len(good_sites)
			n_bases_covered += r1_end - r1_start
			n_bases_covered_cut += r1_end - r1_cut_end
			methy_bins[marker_index]['total_count'] += n_sites_cut
			methy_bins[marker_index]['methylation_count'] += good_methy_status.count('1')
			methy_bins[marker_index]['unmethylation_count'] += good_methy_status.count('0')
			methy_bins[marker_index]['#reads'] += 1
			methy_bins[marker_index]['#fragments'] += 1
		elif strand=='-':
			# reverse single-end read
			r2_start = int(r1_start)
			r2_end = int(r2_end)
			r2_cut_start = max( r2_start, r2_end-length_kept_per_read )
			good_sites_idx = np.where(sites_int>=r2_cut_start)[0]
			good_methy_status = [methy_status[i] for i in good_sites_idx]
			good_sites = [sites_int[i] for i in good_sites_idx]
			n_sites_cut = len(good_sites)
			if n_sites_cut>0:
				if n_columns==11:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r2_cut_start, r2_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), 'NA', read_name) )
				if n_columns==12:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r2_cut_start, r2_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), 'NA', read_name, rg) )
			n_cytosines_cut += len(sites_int) - len(good_sites)
			n_bases_covered += r2_end - r2_start
			n_bases_covered_cut += r2_cut_start - r2_start
			methy_bins[marker_index]['total_count'] += n_sites_cut
			methy_bins[marker_index]['methylation_count'] += good_methy_status.count('1')
			methy_bins[marker_index]['unmethylation_count'] += good_methy_status.count('0')
			methy_bins[marker_index]['#reads'] += 1
			methy_bins[marker_index]['#fragments'] += 1
		else:
			# read pair - fragment
			(r1_end, r2_start) = additional_pos.split(',')
			r1_end = int(r1_end)
			r2_start = int(r2_start)
			r1_start = int(r1_start)
			r2_end = int(r2_end)
			r1_size = r1_end - r1_start
			r2_size = r2_end - r2_start
			r1_cut_end   = min( r1_end, r1_start+length_kept_per_read )
			r2_cut_start = max( r2_start, r2_end-length_kept_per_read )
			fragment_size = r2_end - r1_start
			good_sites_idx = np.unique( np.concatenate( (np.where(sites_int<r1_cut_end)[0], np.where(sites_int>=r2_cut_start)[0]), axis=0) )
			good_methy_status = [methy_status[i] for i in good_sites_idx]
			good_sites = [sites_int[i] for i in good_sites_idx]
			n_sites_cut = len(good_sites)
			if n_sites_cut>0:
				if n_columns==11:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r1_start, r2_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), str(r1_cut_end)+','+str(r2_cut_start), read_name) )
				if n_columns==12:
					sys.stdout.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(m1_idx, m2_idx, chr, r1_start, r2_end, strand, n_sites_cut, ','.join(map(str,good_sites)), ''.join(good_methy_status), str(r1_cut_end)+','+str(r2_cut_start), read_name, rg) )
			n_cytosines_cut += len(sites_int) - len(good_sites)
			if r2_start > r1_end: # no overlap between a pair of old reads
				n_bases_covered += (r1_end-r1_start) + (r2_end-r2_start)
				n_bases_covered_cut += (r1_end-r1_cut_end) + (r2_cut_start-r2_start)
			else: # there are overlaps between a pair of old reads
				n_bases_covered += r2_end - r1_start
				if r2_cut_start > r1_cut_end: # no overlap bewteen a pair of new reads
					n_bases_covered_cut += r2_cut_start-r1_cut_end
				else: # there are overlaps between a pair of new reads
					n_bases_covered_cut += 0 # no change for n_bases_covered_cut, since no bases are cut after trimming
			methy_bins[marker_index]['total_count'] += n_sites_cut
			methy_bins[marker_index]['methylation_count'] += good_methy_status.count('1')
			methy_bins[marker_index]['unmethylation_count'] += good_methy_status.count('0')
			methy_bins[marker_index]['#reads'] += 2
			methy_bins[marker_index]['#fragments'] += 1
	if reads_file!='stdin':
		f.close()
	for marker_index in methy_bins:
		if methy_bins[marker_index]['total_count']==0:
			methy_bins[marker_index]['methylation_rate'] = None
			methy_bins[marker_index]['total_count'] = None
			methy_bins[marker_index]['methylation_count'] = None
			methy_bins[marker_index]['unmethylation_count'] = None
			methy_bins[marker_index]['#reads'] = None
			methy_bins[marker_index]['#fragments'] = None
		else:
			methy_bins[marker_index]['methylation_rate'] = methy_bins[marker_index]['methylation_count'] / float(methy_bins[marker_index]['total_count'])
	return (n_cytosines, n_cytosines_cut, n_bases_covered, n_bases_covered_cut, methy_bins)

#
# An example file of "*.methy_bins" is below
# marker_index	methylation_rate	total_count	methylation_count	unmethylation_count	#reads	#fragments
# 1	0	0	0	0	0	0
# 2	0.277777777777778	36	10	26	11	7
# 3	0.888888888888889	9	8	1	5	2
#
def write_methy_bins_file(file, methy_bins, total_num_markers):
	with open(file,'w') as f:
		f.write('marker_index\tmethylation_rate\ttotal_count\tmethylation_count\tunmethylation_count\t#reads\t#fragments\n')
		for i in range(total_num_markers):
			marker_index = i+1 # marker index is 1-base
			if marker_index in methy_bins:
				if methy_bins[marker_index]['total_count']:
					f.write('%d\t%.15g\t%d\t%d\t%d\t%d\t%d\n'%(marker_index, \
						methy_bins[marker_index]['methylation_rate'], \
						methy_bins[marker_index]['total_count'], \
						methy_bins[marker_index]['methylation_count'], \
						methy_bins[marker_index]['unmethylation_count'], \
						methy_bins[marker_index]['#reads'], \
						methy_bins[marker_index]['#fragments']))
				else:
					f.write('%d\t0\t0\t0\t0\t0\t0\n'%marker_index)
			else:
				f.write('%d\t0\t0\t0\t0\t0\t0\n'%marker_index)

readsFile = sys.argv[1]
length_kept_per_read = int(sys.argv[2])
total_num_markers = int(sys.argv[3])
methy_bins_file = sys.argv[4]
logFile = sys.argv[5]

n_cytosines, n_cytosines_cut, n_bases_covered, n_bases_covered_cut, methy_bins = trim_paired_reads( readsFile, length_kept_per_read )
write_methy_bins_file(methy_bins_file, methy_bins, total_num_markers)

with open(logFile, 'w') as f:
	f.write('Old reads file: %s\n'%readsFile)
	f.write('Number_of_cytosines: %d\n'%n_cytosines)
	f.write('read_length_kept: %d\n'%length_kept_per_read)
	f.write('Number_of_bases(non-overlap): %d\n'%n_bases_covered)
	f.write('Number_of_cytosines: %d\n'%n_cytosines)
	f.write('New reads file: stdout\n')
	f.write('Number_of_bases_cut(non-overlap): %d\n'%n_bases_covered_cut)
	f.write('Number_of_cytosines_cut: %d\n'%n_cytosines_cut)
	f.write('Fraction_of_bases_cut(non-overlap): %.5g\n'%(n_bases_covered_cut/float(n_bases_covered)))
	f.write('Fraction_of_cytosines_cut: %.5g\n'%(n_cytosines_cut/float(n_cytosines)))

