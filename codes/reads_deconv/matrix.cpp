#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace boost;

#include "template_utils.cpp"

Matrix_Double::Matrix_Double(const unsigned int nrow_, const unsigned int ncol_, const double init_value)
{
	nrow = nrow_;
	ncol = ncol_;
	empty = false;
	mat.reserve(nrow_);
	for (int i=0; i<nrow; i++) {
		vector<double> temp;
		temp.reserve(ncol);
		for (int j=0; j<ncol; j++)
			temp.push_back(init_value);
		mat.push_back( temp );
	}
	empty = mat.empty();
}

// a line represents a row of the matrix. Each element in this row is delimited by a TAB.
void Matrix_Double::process_one_line_of_mat_file(string & line, int num_header_column)
{
	vector<string> items;
	split(items, line, is_any_of("\t"));
	vector<double> v;
	for (int i=0; i<items.size(); i++) {
		if (i==0) row_labels.push_back( atoi(items[0].c_str()) ); // the first column is the label of this row.
		if (i<num_header_column) continue; // skip the first num_header_column header columns
		v.push_back( atof(items[i].c_str()) );
	}
	mat.push_back( v );
}

Matrix_Double::Matrix_Double(const string & mat_file,
	bool header_line=false, int num_header_column=1)
{
	unsigned long i_line=0;
	// input is a plain text file
	istream * in=&cin; // default is stdin
	ifstream fin;
	if ( !(mat_file.compare("stdin")==0)) {
		fin.open(mat_file.c_str());
		if (fin.fail()){
			cerr << "Error: Unable to open " << mat_file << " in Matrix_Double()" << endl;
			exit(EXIT_FAILURE);
		}
		in = &fin;
	}
	string line;
	while (!(*in).eof()) {
		getline((*in), line);
		//cerr << "Line: " << line << endl;
		if (header_line) {
			if (i_line==0) {
				// skip the header line
				i_line++;
				continue;
			}
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		process_one_line_of_mat_file(line, num_header_column);
		i_line++;
	}
	if ( !(mat_file.compare("stdin")==0)) {
		fin.close();
	}
	empty = mat.empty();
	nrow = mat.size();
	if (nrow>=1) ncol = mat[0].size();
	else ncol = 0;
}

bool Matrix_Double::get_element(const int i, const int j, double & v)
{
	if (!empty) {
		if (i<nrow && j<ncol) {
			v = mat[i][j];
			return true;
		} else {
			cerr << "Warning: i=" << i << " and j=" << j << " exceed matrix size!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

bool Matrix_Double::set_element(const int i, const int j, double v)
{
	if (!empty) {
		if (i<nrow && j<ncol) {
			mat[i][j] = v;
			return true;
		} else {
			cerr << "Warning: i=" << i << " and j=" << j << " exceed matrix size!" << endl;
			return false;
		}
	} else {
		return false;
	}
}

bool Matrix_Double::get_column_sum(const int j, double & v)
{
	if (!empty) {
		if (j<ncol) {
			v = 0;
			for (int k=0; k<nrow; k++) {
				v += mat[k][j];
			}
			return true;
		} else {
			cerr << "Warning: j=" << j << " exceeds matrix column size (" << ncol << ")!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

bool Matrix_Double::get_row_sum(const int i, double & v)
{
	if (!empty) {
		if (i<nrow) {
			v = 0;
			for (int k=0; k<ncol; k++) {
				v += mat[i][k];
			}
			return true;
		} else {
			cerr << "Warning: i=" << i << " exceeds matrix row size (" << nrow << ")!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

void Matrix_Double::get_unique_row_labels(vector<int> & uniq_labels)
{
	map<int, int> counts;
	for (int i=0; i<row_labels.size(); i++)
		if (counts.find(row_labels[i]) != counts.end())
			counts[row_labels[i]] += 1; // found this label, then increment count
		else
			counts[row_labels[i]] = 1; // not found this label, then count it once
	map<int,int>::iterator it;
	for (it=counts.begin(); it!=counts.end(); ++it)
		uniq_labels.push_back(it->first);
}


ostream & operator<<(ostream & os, Matrix_Double & mat)
{
	double v;
	if (!mat.isempty()) {
		for (int i=0; i<mat.get_row_num(); i++) {
			int j;
			for (j=0; j<mat.get_column_num()-1; j++) {
				mat.get_element(i,j,v);
				os << v << "\t";
			}
			mat.get_element(i,j,v);
			os << v << endl;
		}
	}
	return os;
}

double objective_em_supervise(Matrix_Double & p, vector<double> & theta)
{
	unsigned int ncol = p.get_column_num();
	unsigned int nrow = p.get_row_num();
	double v;
	double obj = 0;
	for (int i=0; i<nrow; i++) {
		double sum = 0;
		for (int j=0; j<ncol; j++) {
			p.get_element(i,j,v);
			sum += theta[j]*v;
		}
		obj += log(sum);
	}
	return obj;
}
//
// output: theta (model parameters), a vector with "ncol" elements. This vector is already allocated with space of ncol units
void em_supervise(Matrix_Double & p, int max_iter, vector<double> & theta)
{
	cout.precision(15);
	cerr.precision(15);
	unsigned int ncol = p.get_column_num();
	unsigned int nrow = p.get_row_num();
	if (theta.size() != ncol) {
		cerr << "Warning (em_supervise): theta size (" << theta.size() << ") is inequal to the column number of matrix 'p' (" << p.get_column_num() << ")!" << endl;
		return;
	}
	// initialize model parameters as uniform distribution
	// alternatively, model parameters can be random numbers by satisfying the crition
	//    (1) \sum_{i=1}^{ncol}{theta_i}=1
	for (int j=0; j<ncol; j++) {
		theta[j] = 1/(double)ncol;
	}
	// create and initialize q with the same size of p and with all elements initialized as 0
	Matrix_Double q(nrow, ncol, 0);
	cerr << "iter 0\t" << objective_em_supervise(p, theta) << endl;
	for (int iter=0; iter<max_iter; iter++) {
		// E-step: estimate q
		double v1, v2;
		for (int i=0; i<nrow; i++) {
			double sum = 0;
			for (int j=0; j<ncol; j++) {
				p.get_element(i,j,v1);
				v2 = theta[j]*v1;
				q.set_element( i, j, v2 );
				sum += v2;
			}
			for (int j=0; j<ncol; j++) {
				q.get_element(i,j,v2);
				v2 /= sum;
				q.set_element( i, j, v2 );
			}
		}
		// M-step: estimate theta
		for (int j=0; j<ncol; j++) {
			double sum=0;
			for (int i=0; i<nrow; i++) {
				q.get_element(i,j,v2);
				sum += v2;
			}
			theta[j] = sum / nrow;
		}
		cerr << "iter " << (iter+1) << "\t" << objective_em_supervise(p, theta) << endl;
	}
}

//
// There are T known tissues and 1 unknown tissue (described by a double vector "m")
//
// input:
//   p is a matrix of N X T, where N is number of reads and T is number of known tissue
//     p.row_label is an int vector of N X 1, each element is the marker bin index (1-base) of a read.
//   Rm is a vector of N X 1, each element is number of valid methylated CpG sites in a read
//   Rl is a vector of N X 1, each element is number of all valid CpG sites in a read
//
// output:
//   theta is a vector of (T+1) X 1, where T is number of known tissue. This vector is already allocated with space of (ncol+1) units.
//   m is a vector of M X 1, where M is number of marker bins. This vector is already allocated with space of nbin units.
//
void em_semisupervise(Matrix_Double & p, vector<int> & Rm, vector<int> & Rl,
	int max_iter, vector<double> & theta, vector<double> & m)
{
	cout.precision(15);
	cerr.precision(15);
	unsigned int ncol = p.get_column_num(); // T, number of known tissues
	unsigned int nrow = p.get_row_num(); // N, number of reads
	int nbin = m.size();

	// initialize model parameters theta as uniform distribution, and m as 0.5
	// alternatively, theta and m can be random numbers, by satisfying:
	//    (1) \sum_{i=1}^{ncol+1}{theta_i}=1
	//    (2) 0 <= m_k <=1 for all k=1,2,...,#bin
	for (int j=0; j<ncol+1; j++) {
		theta[j] = 1/(double)(ncol+1);
	}
	for (int k=0; k<nbin; k++) {
		m[k] = 0.5;
	}
	// create and initialize q with the same size of p and with all elements initialized as 0
	Matrix_Double q(nrow, ncol, 0);
	// create and initialize q_unknown (for the unknown tissue) with nrow and with all elements being 0
	vector<double> q_unknown;
	for (int i=0; i<nrow; i++)
		q_unknown.push_back(0);

	// EM algorithm
	for (int iter=0; iter<max_iter; iter++) {
		// E-step: estimate q (for T known tissues) and q_unknown (for one unknown tissue)
		double v1, v2;
		int bin_index;
		for (int i=0; i<nrow; i++) {
			// process the first T known tissues
			double sum = 0;
			for (int j=0; j<ncol; j++) {
				p.get_element(i,j,v1);
				v2 = theta[j]*v1;
				q.set_element( i, j, v2 );
				sum += v2;
			}
			// process the last unknown tissue
			p.get_row_label(i, bin_index); // get bin_index (1-base) of read i
			double v2 = pow(m[bin_index-1],Rm[i]) * pow(1-m[bin_index-1],Rl[i]-Rm[i]);
			q_unknown[i] = theta[ncol]*v1;
			sum += q_unknown[i];

			// update q and q_unknown
			for (int j=0; j<ncol; j++) {
				q.get_element(i,j,v2);
				v2 /= sum;
				q.set_element( i, j, v2 );
			}
			q_unknown[i] /= sum;
		}

		// M-step 1: estimate m of each marker bin
		double sum1=0, sum2=0;
		int curr_bin_index=-1, prev_bin_index=-1; // all bin index is 1-base
		for (int i=0; i<nrow; i++) {
			// for each row (or read)
			p.get_row_label(i, curr_bin_index); // get the bin_index (1-base) of read i
			if (curr_bin_index != prev_bin_index) {
				// this is the 1st read of a new bin
				// we need to summarize m value of previous bin
				if (prev_bin_index!=-1) { // current bin is not the 1st bin
					m[prev_bin_index-1] = (sum2!=0 ? sum1/sum2 : 0);
					sum1 = 0;
					sum2 = 0;
				}
				prev_bin_index = curr_bin_index;
			}
			sum1 += q_unknown[i]*Rm[i];
			sum2 += q_unknown[i]*Rl[i];
		}
		m[curr_bin_index-1] = (sum2!=0 ? sum1/sum2 : 0); // estimate m of the last bin

		cout << "iter=" << iter << ", " << "m=" << endl; // for debug
		print_vec(cout, m, "\n"); // for debug
		cout << endl; // for debug

		// M-step 2: estimate theta, which has ncol+1 values
		double sum;
		for (int j=0; j<ncol; j++) {
			sum=0;
			for (int i=0; i<nrow; i++) {
				q.get_element(i,j,v2);
				sum += v2;
			}
			theta[j] = sum / nrow;
		}
		sum=0;
		for (int i=0; i<nrow; i++)
			sum += q_unknown[i];
		theta[ncol] = sum / nrow;

		cout << "theta="; // for debug
		print_vec(cout, theta, ", "); // for debug
		cout << endl << endl; // for debug
	}
}

void em_unsupervise(Matrix_Double & p)
{
}

