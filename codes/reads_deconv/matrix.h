#ifndef MATRIX_H
#define MATRIX_H

#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>

using namespace std;

class Matrix_Double {
private:
	vector<vector<double> > mat; // a double matrix of size nrow X ncol
	vector<int> row_labels; // a int vector of size nrow X 1. The first column (integer) is assumed to be the label of a row. 
	unsigned int nrow, ncol;
	bool empty;
	void process_one_line_of_mat_file(string & line, int num_header_column);
	// putting constructors as private can make a class non-copyable
public:
	Matrix_Double() {nrow=0; ncol=0; empty=true;};
	Matrix_Double(const unsigned int nrow_, const unsigned int ncol_, const double init_value);
	Matrix_Double(const string & mat_file, bool header_line, int num_header_column);
	unsigned int get_row_num() {return nrow;};
	unsigned int get_column_num() {return ncol;};
	bool get_column_sum(const int j, double & v);
	bool get_row_sum(const int i, double & v);
	bool get_element(const int i, const int j, double & v);
	bool set_element(const int i, const int j, double v);
	bool isempty() {return empty;};
	bool get_row_label(const int row_index, int & row_label) { // row_index is 0-base
		if (!empty) { row_label = row_labels[row_index]; return true; }
		else return false;};
	vector<int> & get_row_labels() {return row_labels;};
	void get_unique_row_labels(vector<int> & uniq_labels);
};

ostream & operator<<(ostream & os, Matrix_Double & mat);
void em_supervise(Matrix_Double & p, int max_iter, vector<double> & theta);
void em_semisupervise(Matrix_Double & p, vector<int> & Rm, vector<int> & Rl,
	int max_iter, vector<double> & theta, vector<double> & m);
void em_unsupervise(Matrix_Double & p);

#endif
