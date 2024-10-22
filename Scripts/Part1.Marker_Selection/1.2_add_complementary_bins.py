# Modify from Wenyuan's version
# py2 to py3, change print to print()
import sys, argparse, os

# File format is as below
#
# chrom	size
# chr1	249250621
# chr2	243199373
# chr3	198022430
# chr4	191154276
#
def read_chrom_sizes_file(file):
	print("read", file)
	chr2size = {}
	n_line = 0
	with open(file) as f:
		for line in f:
			if n_line==0: # skip the header line
				n_line += 1
				continue
			else:
				items = line.rstrip().split('\t')
				chr2size[items[0]] = int(items[1])
				n_line += 1
	print("  ", len(chr2size), "chromosomes")
	return chr2size

# File format of (Wenyuan defined biomarker paired-value features file, which include paired values for each tissue/cancer type) is as below
# Column 1: marker_index (1-base)
# Column 2: chr
# Column 3: start_coordinate (1-base)
# Column 4: end_coordinate (1-base)
# Column 5: marker_value
# Column 6+: tissue/cancer types, each column is a paired values for a tissue/cancer type. The paired values are separated by ","
#
# marker_index	chr	start	end	marker_value	background_normal_pnas2013	kirc_tumor
# 13	chr1	954618	957121	0.4830000000000001	335.68562692221246,2279.7068695500316	607709.111065718,8522078.726067608
# 21	chr1	1065863	1068091	0.41595000000000004	6219.118948814673,1188.5918976913854	38484.38905220434,17080.123767764366
# 22	chr1	1072196	1074335	0.4212499999999999	22554.477075674793,20807.32981052174	2435.522566726656,1616.5051198928027
# 23	chr1	1098000	1100559	0.47675000000000006	14365.577873353086,4674.759939353441	2249.4131243131133,1733.4210284819076
# ...
# 11526	chr22	44350859	44351687	0.39075000000000004	2025.8369941783471,1275529.0241993738	535935.1276418628,1.1925933626528783E7
# 11527	chr22	44568336	44568915	0.41380000000000006	765.8499869567707,1860.8681700329187	1361.083539897045,837.5956210975257
# 11540	chr22	46067073	46068544	0.36839999999999995	5816.301165832917,83596.47617210586	23662.881885284492,141733.55210015707
# 11544	chr22	46465159	46469093	0.5078499999999999	6090.436486501831,175755.39338149034	5896.097907874386,78683.46774263895
# 11555	chr22	50451053	50454138	0.43284999999999996	4222.501038162458,1733.3516650616589	5359.6442211807225,6424.6368375856755
#
# We suppose its chromosome coordinates are 1-base
#
class Biomarker:
	# "start_column_of_genome_location" is the column of chr, then followed by two columns with start and end coordinates (1-base)
	# column index is 1-base
	def __init__(self,str_line, start_column_of_genome_location=1):
		items = str_line.split('\t')
		try:
			self.chr = items[start_column_of_genome_location-1]
			self.start_coord = int(items[start_column_of_genome_location])
			self.end_coord = int(items[start_column_of_genome_location+1])
			self.line = str_line
		except Exception:
			sys.stderr.write("Error: %s\n"%str_line.rstrip())
	def get_chr(self):
		return self.chr
	def get_coords(self):
		return (self.start_coord, self.end_coord)
	def get_line(self):
		return self.line

# "start_column_of_genome_location" is the column of chr, then followed by two columns with start and end coordinates (1-base)
# column index is 1-base
def read_wenyuan_defined_biomarker_feature_file(file, start_column_of_genome_location):
	print("read", file)
	chr2biomarkers = {} # map: chr -> a list of tuples, where a tuple is (start_coord, end_coord, line)
	n_line = 0
	with open(file) as f:
		for line in f:
			if n_line==0: # skip the header line
				header_line = line.rstrip()
				items = line.rstrip().split('\t')
				# "start_column_of_genome_location" refers to the column of "chr", then
				# followed by 2 columns "start_coordinate" and "end_coordinate"
				num_of_addtional_columns = len(items)-(start_column_of_genome_location+2)
				n_line += 1
				continue
			else:
				biomarker = Biomarker(line.rstrip(), start_column_of_genome_location)
				chr = biomarker.get_chr()
				(start_coord, end_coord) = biomarker.get_coords()
				if not (chr in chr2biomarkers): # this is a new chr
					chr2biomarkers[chr] = []
				markers = chr2biomarkers[chr]
				markers.append( (start_coord, end_coord, biomarker.line) )
				n_line += 1
	print("  ", (n_line-1), "biomarkers")
	return (chr2biomarkers, header_line, num_of_addtional_columns)

def merge_two_biomarkers(chr2biomarkers_I, chr2biomarkers_II):
	chr2biomarkers = chr2biomarkers_I
	n_markers = 0
	for chr in chr2biomarkers_II:
		if not (chr in chr2biomarkers): # this is a new chr
			chr2biomarkers[chr] = chr2biomarkers_II[chr]	
		else:
			markers = chr2biomarkers[chr]
			markers += chr2biomarkers_II[chr]
		chr2biomarkers[chr].sort()
		n_markers += len(chr2biomarkers[chr])
	print("  ", n_markers, "merged biomarkers")
	return chr2biomarkers

# an empty or complementary bin is as below:
# Column 1: "0"
# Column 2: chr
# Column 3: start_coordinate (1-base)
# Column 4: end_coordinate (1-base)
# Column 5+: these additional columns are "-", indicating empty/blank columns 
#
# for example:
# 0	chr1	954618	957121	-	-	-	-	-	-
#
def make_a_bin(chr, start_coord, end_coord, num_of_addtional_columns):
	columns = ["0"] + [chr] + [str(start_coord)] + [str(end_coord)] + ["-"]*num_of_addtional_columns
	line = "\t".join(columns)
	return (start_coord, end_coord, line)

def complement_biomarkers(chr2biomarkers, chr2size, num_of_addtional_columns):
	chr2complementary_bins = {}
	prev_chr = ""
	for chr in chr2biomarkers:
		if not (chr in chr2complementary_bins): # this is a new chr
			chr2complementary_bins[chr] = []
			# complement the last bin of previous chr
			if prev_chr: # prev_chr is not an empty string: this is not the 1st chr
				if prev_end_coord < chr2size[prev_chr]:
					# bin's range is [start_coord, end_coord)
					compl_bin_start_coord = prev_end_coord
					compl_bin_end_coord = chr2size[prev_chr] + 1
					chr2complementary_bins[prev_chr].append( make_a_bin(prev_chr, compl_bin_start_coord, compl_bin_end_coord, num_of_addtional_columns) )
			prev_chr = chr
			
		markers = chr2biomarkers[chr]
		i_marker = 0
		prev_end_coord = 1 # We suppose coordinates start from 1
		for (start_coord, end_coord, line) in markers:
			i_marker += 1
			if start_coord > prev_end_coord:
				compl_bin_start_coord = prev_end_coord
				compl_bin_end_coord = start_coord
			else: 
				sys.stderr.write("Warn: marker (%s:%d-%d) has overlap with its previous marker\n"%(chr,start_coord,end_coord))
			if end_coord <= prev_end_coord:
				sys.stderr.write("Warn: marker (%s:%d-%d) is a subset of previous marker\n"%(chr,start_coord,end_coord))
			else:	
				prev_end_coord = end_coord
			chr2complementary_bins[chr].append( make_a_bin(chr, compl_bin_start_coord, compl_bin_end_coord, num_of_addtional_columns) )

	# this is for the last complementary bin of the last chr
	if prev_end_coord<chr2size[prev_chr]:
		# bin's range is [start_coord, end_coord)
		compl_bin_start_coord = prev_end_coord
		compl_bin_end_coord = chr2size[prev_chr] + 1
		chr2complementary_bins[prev_chr].append( make_a_bin(chr, compl_bin_start_coord, compl_bin_end_coord, num_of_addtional_columns) )

	return chr2complementary_bins

def make_complete_chromosomes():
	all_chrs=[]
	for i in range(1,23):
		all_chrs.append( "chr" + str(i) )
	all_chrs.append("chrX")
	all_chrs.append("chrY")
	return all_chrs

def write_markers(chr2biomarkers, header_line, output_file):
	all_chrs = make_complete_chromosomes()
	with open(output_file,'w') as f:
		f.write(header_line+'\n')
		marker_count = 0
		for chr in all_chrs:
			if chr in chr2biomarkers:
				markers = chr2biomarkers[chr]
				for (start_coord, end_coord, line) in markers:
					f.write( line + '\n' )

parser = argparse.ArgumentParser(formatter_class= \
                    argparse.ArgumentDefaultsHelpFormatter, \
                    description="add complementary/blank bins to existing biomarker features; so that all bins and features must cover every coordinate of the genome.")
parser.add_argument("--chromosome_sizes_file", required=True)
parser.add_argument("--biomarker_features_file", required=True, help="biomarker feature file are genomic regions with a set of associated values")
parser.add_argument("--output_file", required=True)
args = parser.parse_args()

chr2size = read_chrom_sizes_file(args.chromosome_sizes_file)
(chr2biomarkers, header_line, num_of_addtional_columns) = read_wenyuan_defined_biomarker_feature_file(args.biomarker_features_file, 2)
chr2complementary_bins = complement_biomarkers(chr2biomarkers, chr2size, num_of_addtional_columns)
chr2all = merge_two_biomarkers(chr2biomarkers, chr2complementary_bins)
write_markers(chr2biomarkers, header_line, args.output_file)

print("Done.")
