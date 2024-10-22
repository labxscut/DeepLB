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
	{ "reads_prob_file",                required_argument, 0,   'p' }, /*       2 */
	{ "em_algorithm",                   required_argument, 0,   'a' }, /*       3 */
	{ "em_max_iterations",              required_argument, 0,   'x' }, /*       3 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hp:a:x:"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Perform EM algorithm to obtain model parameters, using reads binning file and reads probability file" << endl << endl;
	cerr << "USAGE: em [options] -p [reads probability file]" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -a [string]:  EM algorithm type: supervise (default), semi-supervise, or unsupervise." << endl;
	cerr << "     -x [positive integer]:  EM algorithm maximum iteration number: 50 (default)." << endl;
	cerr << endl;
    cerr << "version 1.0, January 20, 2016" << endl;
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
			case 'p': /* --reads_prob_file */
				reads_prob_file = optarg;
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
	if (reads_prob_file.empty()) {
		cerr << endl << "Error: reads probability file is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (em_algorithm_type.compare("supervise")!=0 && em_algorithm_type.compare("semi-supervise")!=0 && em_algorithm_type.compare("unsupervise")!=0) {
		cerr << endl << "Error: EM algorithm types should take one of three values: supervise, semi-supervise and unsupervise!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);

	cerr << "reads probability file: " << reads_prob_file << endl;
	cerr << "EM algorithm type: " << em_algorithm_type << endl;
	cerr << "EM max iterations: " << em_max_iterations << endl;

	Matrix_Double p(reads_prob_file, true, 1); // has a header line and a header column
	unsigned int n_tissue = p.get_column_num();
	vector<double> theta;
	theta.resize(n_tissue, 0);
	em_supervise(p, em_max_iterations, theta);
	print_vec(cout, theta);
	cout << endl;
}
