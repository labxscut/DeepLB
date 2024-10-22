#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "data_types.h"

using namespace std;

void make_complete_chromosomes(vector<string> & all_chrs);

void read_list_strings_from_file(string input_file, vector<string>& strs);
void read_two_columns_of_list_strings_from_file(string input_file, vector<string>& strs1, vector<string>& strs2);
void read_two_columns_of_list_strings_from_file(string input_file, vector<string>& strs1, vector<string>& strs2,
	map<string, vector<string> > & map_str1tostrs2, string delimit);

void read_wig_file(string wig_file, map<string,map<int,double> > & data);
void print_wig_data(map<string,map<int,double> > & data);
void binning_wig_data(map<string,map<int,double> > & data,
	Bins_end_coord & bins_end_coord, Bins_index & bins_index,
	Wig_bins_data & wig_bins_data);
void print_wig_bins_data(Wig_bins_data & wig_bins_data, Bins_index & bins_index,
	Bins_info & bins_info);

void read_value_file_of_bins(string value_file_of_bins,
	unsigned int value_column_start_index, Bins2Values & bins2value,
	vector<string> & value_names);
void read_paired_values_file_of_bins(string value_file_of_bins,
	unsigned int value_column_start_index, Bins2PairedValues & bins2pairedvalues,
	vector<string> & value_names);
void read_position_value_pairs_file_of_bins(string position_value_pairs_file_of_bins,
	Bins2PositionValuePairs bins2positionvaluepairs);

void calc_read_probability_by_bins2values(string reads_binning_file,
	Bins2Values & bins2values, vector<string> & value_names);
void calc_read_probability_by_methy_level_bins(string reads_binning_file, int num_bins_of_methylation_level);
void calc_read_probability_by_bins2pairedvalues(string reads_binning_file,
	Bins2PairedValues & bins2pairedvalues, vector<string> & value_names);
void get_reads_methy_data_from_reads_binning_file(string reads_binning_file,
	vector<int> & num_cpg_sites, vector<int> & num_methy_cpg_sites);

void print_vec_of_uint(ostream& of, vector<unsigned int> & v);
void print_vec_of_ulong(ostream& of, vector<unsigned long> & v);
void print_map_of_strings(ostream& of, map<string, vector<string> > & map_str1tostrs2);
void print_str_vectors(ofstream& o, vector<string> & s);
void strings2floats(string str, vector<float> & vec, string delimit=",");
float min(float a, float b);
float max(float a, float b);
float meta_p(vector<float> & probs);

#endif
