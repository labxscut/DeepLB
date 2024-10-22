// reads_binning.cpp
// Given a reads file (SAM), print those reads into biomarker bins which are defined by genome bins annotation files where only biomarker bins have non-zero indexes.
//
// 2016/01/08
// Author: Wenyuan Li
//
// Updates:
//   2016/02/15
//      1. Generate too many warning messages. Make them clean and short
//      2. "segment fault" when using option "-B or --output_bins_methy_counts"
//
//   2017/05/03
//      1. Major change: add the function that can pair the forward and reverse reads to output fragments. The format of output fragments is defined. In "sam.h" and "sam.cpp", two new classes are defined (ReadProcessed and ReadPair). Both classes are used only for processing paired-end sequencing data, when filter_strand="both" and is_paired_end_reads=true.
//      2. "DEBUG" and "DEBUG2" are defined for future debug use in "sam.cpp"
//
//   2017/05/18
//      1. Add a strand filter (variable "filter_strand") to select reads with different types: "both", "forward-only", "reverse-only". The function "binning_reads_of_sam_file_of_specific_strand_with_bins_methycounts()" is added to "sam.cpp"
//      question: the function "binning_reads_of_sam_file()" does not consider removing overlaps of forward and matched reverse reads (paired-end sequencing data), does not consider the strand filter either. This will be resolved in future updates. Please always use the option "-B <bin-level methylation counts file>"
//
//	2017/05/30
//		1. Add an option to output SAM record for each single-end read to a separately supplementary file. For paired-end reads, output SAM records to two separately supplementary files (read1 and read2).
//
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <climits>
#include <iomanip> // std::setprecision
#include <iostream> // std::cout, std::fixed
#include <fstream>
#include <sstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "my_getopt-1.5/getopt.h"
#include "sam.h"
#include "utils.h"
#include "data_types.h"

using namespace std;
using namespace boost;

// type definitions
typedef map<string,map<int,double> > POSITIONS_VALUES; // chr -> map(position -> value)

// Global variables
string bins_annot_file="";
string reads_file="stdin"; // default input file is "standard input" (stdin)
string reads_file_format="SAM"; // default is "SAM" file format
int min_cpg_sites=3; // default minimum CpG sites of reads we will use is 3.
string output_file_bins_methy_counts="";
string output_file_of_sam_reads="";
string strand = "both"; // default is "both" (using all reads of both forward and reverse strands); the other values are "forward-only" (using only forward strand) and "reverse-only" (using only reverse strand).
bool is_paired_end_reads = true; // if filter_strand="forward-only" or "reverse-only", is_paired_end_reads is is of no use (so it is disabled).
bool print_RG_tag=false;
bool do_not_print_reads=false;
int print_level=1; // default is 1, the other optins are 2, 3, 11
bool correct_coordiate_of_reverse_reads = false; // If true, correcting the methylated C in the reverse strand, this "start" cannot be used for the leftmost coordinate. This is only for Dennis Lo's Bismark bam file

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	{ "bins_annot_file",                required_argument, 0,   'b' }, /*       1 */
	{ "reads_file",                     required_argument, 0,   'r' }, /*       2 */
	{ "reads_file_format",              required_argument, 0,   'f' }, /*       3 */
	{ "filter_min_cpg_sites_of_reads",  required_argument, 0,   'c' }, /*       4 */
	{ "output_bins_methy_counts",       required_argument, 0,   'B' }, /*       5 */
	{ "output_SAM_reads",               required_argument, 0,   'R' }, /*       5 */
	{ "strand",                         required_argument, 0,   'S' }, /*       5 */
	{ "single_end_reads",  				no_argument,       0,   's' }, /*       0 */
	{ "print_RG_tag",     				no_argument,       0,   'g' }, /*       0 */
	{ "print_level",                    required_argument, 0,   'P' }, /*       5 */
	{ "do_not_print_reads",     		no_argument,       0,   'p' }, /*       0 */
	{ "correct_reverse_coordinate",     required_argument, 0,   'C' }, /*       5 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hb:B:c:f:gpP:r:R:sS:C"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Binning the reads data and print these reads to standard output" << endl << endl;
	cerr << "USAGE: reads_binning [options] -b [genome bin file] -r [reads file] -f [reads file format]" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -f [reads file format]:        default value is SAM." << endl;
	cerr << "     -c [min number of CpG sites]:  default value is 3. This is for filtering reads covering small number of CpG sites." << endl;
	cerr << "     -S :  output reads of 'both' (default), 'forward-only', and 'reverse-only' strand." << endl;
	cerr << "     -s :  input reads are single end reads." << endl;
	cerr << "     -g :  print RG tag of the reads." << endl;
	cerr << "     -P :  print level 1 (default, for likelihood calc), 2, 3, 11." << endl;
	cerr << "     -p :  do not print the reads. This option disables the option -P" << endl;
	cerr << "     -B [output file of (un)methylation values of bins]:  Optional. If not provided, then won't output this file" << endl;
	cerr << "     -R [output file of SAM reads]:  Optional. If not provided, then won't output this file. For single-end reads, file extension is '*.sam', for paired-end reads, file extensions are '*.sam2', '*.sam2'" << endl;
	cerr << "     -C :  Correct the coordinates of CpG sites in reverse reads only for Dennis Lo's Bismark SAM file of single-end reads (not for paired-end reads): minus 1 for all sites." << endl;
	cerr << endl;
    cerr << "version 1.1, August 9 2016" << endl;
    cerr << "version 1.2, May 3 2017" << endl;
    cerr << "version 1.3, May 30 2017" << endl;
    cerr << "By Wenyuan Li" << endl;
	cerr << endl;
}

void parse_command_line(int argc, char * argv[])
{
	int opt; /* during argument parsing, opt contains the return value from getopt() */
	int longind = 0; /* long option list index */
	/* parse all options from the command line */
	while ((opt=getopt_long(argc, argv, shortopts, longopts, &longind)) != -1) {
		switch (opt) {
			case 'h': /* --help */
				print_usage(); /* 'parms_setting' is a global variable */
				exit(EXIT_FAILURE);
				break;
			case 'b': /* --bins_annot_file */
				bins_annot_file = optarg;
				break;
			case 'r': /* --reads_file */
				reads_file = optarg;
				break;
			case 'f': /* --reads_file_format */
				reads_file_format = optarg;
				break;
			case 'c': /* --filter_min_cpg_sies_of_reads */
				min_cpg_sites = (int)atoi(optarg);
				break;
			case 'B': /* --output_bins_methy_counts */
				output_file_bins_methy_counts = optarg;
				break;
			case 'R': /* --output_file_prefix_of_sam_reads */
				output_file_of_sam_reads = optarg;
				break;
			case 'S': /* --strand */
				strand = optarg;
				break;
			case 's': /* --single_end_reads */
				is_paired_end_reads = false;
				break;
			case 'g': /* --print_RG_tag */
				print_RG_tag = true;
				break;
			case 'P': /* --print_level */
				print_level = (int)atoi(optarg);
				break;
			case 'p': /* --do_not_print_reads */
				do_not_print_reads = true;
				break;
			case 'C': /* --correct_reverse_coordinate */
				correct_coordiate_of_reverse_reads = true;
				break;
			default: /* something unexpected has happened */
				cerr << endl << "Error: wrong parameters!" << endl << endl;
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
	if (bins_annot_file.empty()) {
		cerr << endl << "Error: bins annotation file is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (reads_file.empty()) {
		cerr << endl << "Error: reads file is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);

	//cerr << "genome bin annotation file: " << bins_annot_file << endl;
	//cerr << "reads file: " << reads_file << endl;
	//cerr << "reads file format: " << reads_file_format << endl;
	//cerr << "min number of CpG sites: " << min_cpg_sites << endl;
	//if (!output_bins_methy_counts.empty()) {
		//cerr << "output_bins_methy_counts: " << output_bins_methy_counts << endl;
	//}

	Bins_end_coord bins_end_coord; // This is used for finding which bin a genomic position belongs to.
	Bins_index bins_index; // marker index (start from 1) and index of complementary bins (usually 0)
	Bins_info bins_info; // the line of each bin, representing all the information of a bin
	read_bins_annot_file(bins_annot_file, bins_end_coord, bins_index, bins_info, true);
	//print_bins( cerr, bins_end_coord, bins_index, bins_info);

	//read_sam_file(reads_file);

	// print_flag (print to stdout):
	//   0: print nothing
	//   1: print only valid CpG sites of a read, which fall into the marker bin.
	//   2: print all CpG sites of a read
	if (do_not_print_reads) print_level=0;
	if (output_file_bins_methy_counts.empty()) {
		binning_reads_of_sam_file(reads_file, bins_end_coord, bins_index, min_cpg_sites, print_level, print_RG_tag, correct_coordiate_of_reverse_reads);
	} else {
		Bins2Values bins2methycounts; // two counts in each bin: methylation counts and unmethylation counts
		vector<int> markers_index; // store only indexes of markers
		int num_marker_bins = get_num_of_non_void_bins(bins_index, markers_index);
		if (num_marker_bins>0) {
			vector<string> columns_names;
			cerr << "#marker: " << num_marker_bins << endl;
			cerr << "read type: ";
			if (is_paired_end_reads)
				cerr << "paired-end" << endl;
			else
				cerr << "single-end" << endl;
			cerr << "output reads with strand type: " << strand << endl;
			cerr << "Note: if strand type is 'forward-only' or 'reverse-only', then read type is disabled." << endl << flush;
			if (strand=="forward-only" || strand=="reverse-only") {
				create_Bins2Values(markers_index, 3, 0.0, bins2methycounts); // 3 elements are used for count
				binning_reads_of_sam_file_of_specific_strand_with_bins_methycounts(reads_file, bins_end_coord, bins_index, min_cpg_sites, strand, print_level, print_RG_tag, correct_coordiate_of_reverse_reads, bins2methycounts, output_file_of_sam_reads);
				columns_names.push_back("marker_index");
				columns_names.push_back("methylation_rate");
				columns_names.push_back("total_count");
				columns_names.push_back("methylation_count");
				columns_names.push_back("unmethylation_count");
				columns_names.push_back("#reads");
				write_Bins2Values(bins2methycounts, columns_names, output_file_bins_methy_counts, true);
				double nRead=0;
				Bins2Values::iterator it;
				for (it=bins2methycounts.begin(); it!=bins2methycounts.end(); ++it) {
					vector<double> & values = it->second;
					nRead += values[2];
				}
				cerr << std::fixed;
				cerr << "#read=" << std::setprecision(0) << nRead << endl;
				exit(EXIT_SUCCESS);
			}
			if (is_paired_end_reads) { // paired-end reads
				create_Bins2Values(markers_index, 4, 0.0, bins2methycounts); // 4 elements are used for count
				binning_paired_end_reads_of_sam_file_with_bins_methycounts(reads_file, bins_end_coord, bins_index, min_cpg_sites, print_level, print_RG_tag, bins2methycounts, output_file_of_sam_reads);
				columns_names.push_back("marker_index");
				columns_names.push_back("methylation_rate");
				columns_names.push_back("total_count");
				columns_names.push_back("methylation_count");
				columns_names.push_back("unmethylation_count");
				columns_names.push_back("#reads");
				columns_names.push_back("#fragments");
				write_Bins2Values(bins2methycounts, columns_names, output_file_bins_methy_counts, true);
				double nRead=0, nFragment=0;
				Bins2Values::iterator it;
				for (it=bins2methycounts.begin(); it!=bins2methycounts.end(); ++it) {
					vector<double> & values = it->second;
					nRead += values[2];
					nFragment += values[3];
				}
				cerr << std::fixed;
				cerr << "#read=" << std::setprecision(0) << nRead << endl;
				cerr << "#fragment=" << std::setprecision(0) << nFragment << endl;
			} else { // single-end reads
				create_Bins2Values(markers_index, 3, 0.0, bins2methycounts); // 3 elements are used for count
				binning_reads_of_sam_file_with_bins_methycounts(reads_file, bins_end_coord, bins_index, min_cpg_sites, print_level, print_RG_tag, correct_coordiate_of_reverse_reads, bins2methycounts, output_file_of_sam_reads);
				columns_names.push_back("marker_index");
				columns_names.push_back("methylation_rate");
				columns_names.push_back("total_count");
				columns_names.push_back("methylation_count");
				columns_names.push_back("unmethylation_count");
				columns_names.push_back("#reads");
				write_Bins2Values(bins2methycounts, columns_names, output_file_bins_methy_counts, true);
				double nRead=0;
				Bins2Values::iterator it;
				for (it=bins2methycounts.begin(); it!=bins2methycounts.end(); ++it) {
					vector<double> & values = it->second;
					nRead += values[2];
				}
				cerr << std::fixed;
				cerr << "#read=" << std::setprecision(0) << nRead << endl;
			}
		}
	}
}

