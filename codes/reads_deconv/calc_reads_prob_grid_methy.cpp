// calc_reads_prob_grid_methy.cpp
// Given a reads-binning file (including methylation data) and number of bins for [0, 1/n, 2/n, ..., 1], we calculate the probability of each read belonging to each bin of methylation level.
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
int num_methy_level_bins=10;

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	{ "input_reads_binning_file",       required_argument, 0,   'r' }, /*       1 */
	{ "num_methy_level_bins",     		required_argument, 0,   'n' }, /*       1 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hr:n:z"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Calculate the probability of each read that belongs to a bin of methylation level. The output is a probability matrix (lines are reads and columns are methylation level bins) and the header line denoting each tissue's name." << endl << endl;
	cerr << "USAGE: calc_reads_prob_grid_methy [options] -r [input reads binning file]" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -n [number of methylation level bins]: e.g., [0, 1/n, 2/n, ..., 1] and the default is 10" << endl;
	cerr << endl;
    cerr << "version 1.0, September 20, 2016" << endl;
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
			case 'n': /* --num_methy_level_bins */
				num_methy_level_bins = (int)atoi(optarg);
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
	if (num_methy_level_bins<=1) {
		cerr << endl << "Error: number of methylation level bins should be >=2!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);

	//cerr << "input reads binning file: " << reads_binning_file << endl;
	//cerr << "number of methylation level bins: " << num_methy_level_bins << endl;

	calc_read_probability_by_methy_level_bins(reads_binning_file, num_methy_level_bins);
}

