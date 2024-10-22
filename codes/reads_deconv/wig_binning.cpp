// wig_binning.cpp
// Given a wig file consisting of position-value pairs, print these into biomarker bins which are defined by genome bins annotation files where only biomarker bins have non-zero indexes.
//
// 2015/12/29
// Author: Wenyuan Li
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
#include "data_types.h"

using namespace std;
using namespace boost;

// type definitions
typedef map<string,map<int,double> > POSITIONS_VALUES; // chr -> map(position -> value)

// Global variables

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "h"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Binning the wig data and print them to standard output" << endl << endl;
	cerr << "USAGE: wig_binning <option> [genome bin file] [wig file]" << endl;
	cerr << endl;
    cerr << "version 1.0, December 29 2015" << endl;
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
			default: /* something unexpected has happened */
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
}


int main(int argc, char* argv[]) {
	string wig_file, bin_annot_file;
	if ((optind+2)>argc) {
		cerr << endl << "Error: #required paramters <> 3!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	parse_command_line(argc, argv);
	bin_annot_file = argv[optind]; optind++;
	wig_file = argv[optind]; optind++;

	//cerr << "genome bin annotation file: " << bin_annot_file << endl;
	//cerr << "wig file: " << wig_file << endl;

	POSITIONS_VALUES wig_data;
	//cerr << "reading wig file ... '" << wig_file << "'" << endl << flush;
	read_wig_file(wig_file, wig_data);
	
	Bins_end_coord bins_end_coord; // This is used for finding which bin a genomic position belongs to.
	Bins_index bins_index; // marker index (start from 1) and index of complementary bins (usually 0)
	Bins_info bins_info; // the line of each bin, representing all the information of a bin
	read_bins_annot_file(bin_annot_file, bins_end_coord, bins_index, bins_info, true);
	//print_bins( cerr, bins_end_coord, bins_index, bins_info);

	//cerr << "binning wig data ..." << endl << flush;
	Wig_bins_data wig_bins_data;
	binning_wig_data(wig_data, bins_end_coord, bins_index, wig_bins_data);

	//cerr << "printing ... " << endl << flush;
	print_wig_bins_data(wig_bins_data, bins_index, bins_info);
	//print_wig_data(wig_data);
}

