// em.cpp
// Given a reads-binning file (including methylation data) and a reads probability file (each line is a read, each column is a tissue type), perform EM algorithm to get model parameters

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
#include "matrix.h"

using namespace std;
using namespace boost;

#include "template_utils.cpp"

// Global variables
string reads_binning_file=""; // this is only used by semi-supervised or unsupervised EM algorithms
string reads_prob_file="";
string em_algorithm_type="supervise"; // default is "supervised manner"
int em_max_iterations=50; // default is 50 iterations

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	{ "reads_binning_file",             required_argument, 0,   'r' }, /*       2 */
	{ "em_algorithm",                   required_argument, 0,   'a' }, /*       3 */
	{ "em_max_iterations",              required_argument, 0,   'x' }, /*       3 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hr:a:x:"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Perform EM algorithm to obtain model parameters, using reads binning file and reads probability file" << endl << endl;
	cerr << "USAGE: em [options] {reads probability file}" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -a [string]:  EM algorithm type: supervise (default), semi-supervise, or unsupervise." << endl;
	cerr << "     -r [FILE]:  reads binning file must be provided when EM type is semi-supervise or unsupervise." << endl;
	cerr << "     -x [positive integer]:  EM algorithm maximum iteration number: 50 (default)." << endl;
	cerr << endl;
    cerr << "version 1.0, January 28, 2016" << endl;
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
			case 'a': /* --em_algorithm */
				em_algorithm_type = optarg;
				break;
			case 'x': /* --em_max_iterations */
				em_max_iterations = (int)atoi(optarg);
				break;
			default: /* something unexpected has happened */
				cerr << endl << "Error: wrong parameters!" << endl << endl;
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
	if (em_algorithm_type.compare("supervise")!=0 && em_algorithm_type.compare("semi-supervise")!=0 && em_algorithm_type.compare("unsupervise")!=0) {
		cerr << endl << "Error: EM algorithm types should take one of three values: supervise, semi-supervise and unsupervise!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (em_algorithm_type.compare("semi-supervise")==0 || em_algorithm_type.compare("unsupervise")==0) {
		if (reads_binning_file.empty()) {
			cerr << endl << "Error: reads binning file must be provided, when EM type=" << em_algorithm_type << "!" << endl << endl;
			print_usage();
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char* argv[]) {
	if ((optind+1)>argc) {
		cerr << "Error: #required paramters <> 1!" << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	parse_command_line(argc, argv);
	reads_prob_file = argv[optind]; optind++;

	cerr << "reads probability file: " << reads_prob_file << endl;
	cerr << "EM algorithm type: " << em_algorithm_type << endl;
	if (em_algorithm_type.compare("semi-supervise")==0 || em_algorithm_type.compare("unsupervise")==0)
		cerr << "reads binning file: " << reads_binning_file << endl;
	cerr << "EM max iterations: " << em_max_iterations << endl;

	//cout << "debug before constructor Matrix_Double" << endl << flush;
	Matrix_Double p(reads_prob_file, true, 1); // has a header line and a header column
	//cerr << "debug after constructor Matrix_Double" << endl << flush;
	unsigned int n_tissue = p.get_column_num();
	vector<double> theta;
	//cout << "debug before em_type.compare" << endl << flush;
	if (em_algorithm_type.compare("supervise")==0) {
		theta.resize(n_tissue, 0); // assign (n_tissue) space and initialize to 0s.
		em_supervise(p, em_max_iterations, theta);
	}
	if (em_algorithm_type.compare("semi-supervise")==0) {
		theta.resize(n_tissue+1, 0); // assign (n_tissue+1) space and initialize to 0s.
		vector<int> uniq_bin_indexes;
		//cout << "debug: before get unique bin indexes " << endl << flush;
		p.get_unique_row_labels(uniq_bin_indexes);
		//cout << "debug: after get unique bin indexes " << endl << flush;
		//print_vec(cout, uniq_bin_indexes, "\n");
		int nbin = uniq_bin_indexes.size(); // M, number of bins that all N reads belong to.
		//cout << "debug: nbin=" << nbin << endl << flush;
		vector<double> m;
		m.resize(nbin, 0); // assign nbin space and initialize to 0s.
		// read reads_binning_file to get Rm and Rl
		vector<int> Rm, Rl;
		get_reads_methy_data_from_reads_binning_file(reads_binning_file, Rl, Rm);
		//print_vec(cout, Rl, "\n");
		em_semisupervise(p, Rm, Rl, em_max_iterations, theta, m);
	}
	print_vec(cout, theta, "\n");
}

