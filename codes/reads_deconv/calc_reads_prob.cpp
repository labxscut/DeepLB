// calc_reads_prob.cpp
// Given a reads-binning file (including methylation data) and a set of tissue files (each file is a tissue), we calculate the probability of each read belonging to the tissue.
//
// tissue file could be
//   (1) position-specific
//   (2) bin-specific with each bin a single value
//   (3) bin-specific with each bin a distribution
//
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "my_getopt-1.5/getopt.h"
#include "utils.h"

using namespace std;
using namespace boost;

// Global variables
string reads_binning_file="";
string tissues_list_file="";
string tissue_data_type="bin_single_value"; // default is "bin_single_value"

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	{ "input_reads_binning_file",       required_argument, 0,   'r' }, /*       1 */
	{ "input_file_of_tissues_list",     required_argument, 0,   't' }, /*       1 */
	{ "tissue_data_type",               required_argument, 0,   'T' }, /*       3 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hr:t:T:"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Calculate the probability of each read that is originated from each tissue. The output is a probability matrix (lines are reads and columns are tissues) with header line denoting each tissue's name." << endl << endl;
	cerr << "USAGE: calc_reads_prob [options] -r [input reads binning file] -t [input file of tissues list]" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -T [tissue data type]:  tissue data type: bin_single_value (default, the tissue file is values_of_bins) or bin_paired_values." << endl;
	cerr << endl;
    cerr << "version 1.0, January 21, 2016" << endl;
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
			case 'r': /* --reads_binning_file */
				reads_binning_file = optarg;
				break;
			case 't': /* --input_file_of_tissues_list */
				tissues_list_file = optarg;
				break;
			case 'T': /* --tissue_file_type */
				tissue_data_type = optarg;
				break;
			default: /* something unexpected has happened */
				cerr << endl << "Error: wrong parameters!" << endl << endl;
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
	if (reads_binning_file.empty()) {
		cerr << endl << "Error: input reads binning file is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (tissues_list_file.empty()) {
		cerr << endl << "Error: input file of tissues list is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (tissue_data_type.compare("bin_single_value")!=0 && tissue_data_type.compare("bin_paired_values")!=0) {
		cerr << endl << "Error: tissue data type should take one of three values: bin_single_value or bin_paired_values!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);

	//cerr << "input reads binning file: " << reads_binning_file << endl;
	//cerr << "file of tissues list: " << tissues_list_file << endl;
	//cerr << "tissue data type: " << tissue_data_type << endl;

	if (tissue_data_type.compare("bin_single_value")==0) {
		Bins2Values bins2values;
		vector<string> tissue_names;
		//
		// value_file_of_bins format: a header line and the following lines are marker bins
		//Column 1: marker_index, 1-based index. Only marker bins are included, those complementary bins do not appear in this file.
		//Column 2: chr
		//Column 3: start coordinate of bin (1-base)
		//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
		//Column 5: marker_type. "I" is Marker_type_I, "II" is Marker_type_II, "-" is the complementary bin, only facilitating the searching.
		//Column 6+: values for this marker (each column is a tissue)
		//
		// for example:
		//marker_index	chr	start	end	marker_type	Liver	T-cells	B-cells	Neutrophils	Placenta
		//1	chr1	855266	855766	II	0.472	0.848	0.607	0.298	0.144
		//2	chr1	969796	970296	II	0.659	0.458	0.352	0.332	0.331
		unsigned int value_column_start_index=6; // column index is 1-base. This is the first value's column index of the file "value_file_of_bins"
		read_value_file_of_bins(tissues_list_file, value_column_start_index, bins2values, tissue_names);
		calc_read_probability_by_bins2values(reads_binning_file, bins2values, tissue_names);
	}
	if (tissue_data_type.compare("bin_paired_values")==0) {
		Bins2PairedValues bins2pairedvalues;
		vector<string> tissue_names;
		// paired_values_file_of_bins format: a header line and the following lines are marker bins
		//Column 1: marker_index, 1-based index. Only marker bins are included, those complementary bins do not appear in this file.
		//Column 2: chr
		//Column 3: start coordinate of bin (1-base)
		//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
		//Column 5: marker_type. "-" is the complementary bin, only facilitating the searching.
		//Column 6+: paired values for this marker (each column is a tissue). Each pair contains two values, delimited by ":"

		//For example:
		//marker_index	chr	start	end	marker_type	normal_plsa	tumor_brca
		//1	chr1	855266	855766	brca	47.2:28.9	8.48:10.7
		//2	chr1	969796	970296	brca	6.59:10.9	5.58:8.9
		//3	chr1	1099044	1099544	brca	7.22:20.8	8.79:3.2
		//4	chr1	1109315	1109815	brca	30.8:20.7	25.8:78.9
		//...
		//5816	chr22	50962109	50962609	I	91.2:10.7	20.7:60.7
		//5817	chr22	50987071	50987571	II	1.04:10.2	15.7:1.8
		//5818	chr22	51016754	51017254	II	32.1:10.7	10.8:2.9
		//5819	chr22	51041995	51042495	II	1.29:1.9	5.7:8.9
		unsigned int value_column_start_index=6; // column index is 1-base. This is the first value's column index of the file "value_file_of_bins"
		read_paired_values_file_of_bins(tissues_list_file, value_column_start_index, bins2pairedvalues, tissue_names);
		calc_read_probability_by_bins2pairedvalues(reads_binning_file, bins2pairedvalues, tissue_names);
	}
}

