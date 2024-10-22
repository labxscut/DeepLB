#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/assign.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_exp.h>
#include "data_types.h"
#include "utils.h"

using namespace std;
using namespace boost;

void make_complete_chromosomes(vector<string> & all_chrs){
	all_chrs.clear();
	for (int i=1; i<=22; i++) {
		ostringstream convert;
		convert << i;
		all_chrs.push_back("chr" + convert.str());
	}
	all_chrs.push_back("chrX");
	all_chrs.push_back("chrY");
}

// File format (each line is a string)
void read_list_strings_from_file(string input_file, vector<string>& strs) {
	ifstream fin;
	fin.open(input_file.c_str());
	if (fin.fail()){
        cerr << "Error: Unable to open " << input_file << " in read_list_strings_from_file()" << endl;
		cerr << "Exit." << endl;
        exit(EXIT_FAILURE);
    }
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		strs.push_back(line);
	}
	fin.close();
}

// File format (each line has two strings, TAB delimited):
// Column 1: string1
// Column 2: string2
//
// We obtain the unique values of Column 1, and build a map between unique value of Column 1 and their corresponding values in Column 2
void read_two_columns_of_list_strings_from_file(string input_file, vector<string>& strs1, vector<string>& strs2,
	map<string, vector<string> > & map_str1tostrs2, string delimit="\t") {
	cerr << "Load file '" << input_file << "'" << endl;
	ifstream fin;
	fin.open(input_file.c_str());
	if (fin.fail()){
        cerr << "Error: Unable to open " << input_file << " in read_two_columns_of_list_strings_from_file()" << endl;
		cerr << "Exit." << endl;
        exit(EXIT_FAILURE);
    }
	string line;
	unsigned long i=0;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		i++;
		vector<string> items;
		split(items, line, is_any_of(delimit));
		if (items.size()<2) {
        	cerr << "Error (read_two_columns_of_list_strings_from_file):"<< endl;
			cerr << "   File: " << input_file << endl;
			cerr << "   Line " << i << ": There are less than 2 columns!" << endl;
			cerr << "Exit." << endl;
	        exit(EXIT_FAILURE);
		}
		strs1.push_back(items[0]);
		strs2.push_back(items[1]);

		// build a map between strs1 and strs2
		// string in strs1 is a key, and strings of strs2 in corresponding lines are values of this key
		if ( map_str1tostrs2.find(items[0]) == map_str1tostrs2.end() ) {
			// this is a new key
			vector<string> values;
			values.push_back(items[1]);
			map_str1tostrs2.insert(make_pair(items[0], values));
		} else {
			// this key already exists in the map
			vector<string> & values = map_str1tostrs2[items[0]];
			values.push_back(items[1]);
		}
	}
	fin.close();
	// print a summary of this map
	map<string, vector<string> >::iterator it;
	int n=0;
	for (it=map_str1tostrs2.begin(); it!=map_str1tostrs2.end(); it++) {
		cout << it->first << "\t" << it->second.size() << endl;
		n += it->second.size();
	}
	cerr << "Total: " << n << " elements" << endl;
}

// wig File format
// http://www.ensembl.org/info/website/upload/wig.html
// http://genome.ucsc.edu/goldenpath/help/wiggle.html
// wig file format: Wiggle element data values can be integer or real, positive or negative. Chromosome positions are 1-relative, i.e. the first base is 1. Only positions specified have data; unspecified positions will be empty.
// Wiggle format is line-oriented. For wiggle custom tracks, the first line must be a track definition line (i.e., track type=wiggle_0), which designates the track as a wiggle track and adds a number of options for controlling the default display.
//
// track type=wiggle_0 name="UCSD.Adipose_Tissue.Bisulfite-Seq.STL003:methRatio" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean
// variableStep  chrom=chrN  [span=windowSize]
// chromStartA  dataValueA
// chromStartB  dataValueB
// ... etc ...  ... etc ...
//
// For example
// track type=wiggle_0 name="UCSD.Adipose_Tissue.Bisulfite-Seq.STL003:methRatio" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean
// variableStep chrom=chr1
// 10469	0.75
// 10470	0.75
// 10471	0.833333333333333
// 10472	0.833333333333333
// 10484	0.928571428571429
// 10485	0.928571428571429
// We obtain the unique values of Column 1, and build a map between unique value of Column 1 and their corresponding values in Column 2
//
// output variable: map<string,map<int,double> >  data, it is "chr -> map(position -> value)"
// The input line is not the first header line
void process_one_line_of_wig_file(string & line, string & chr, map<string,map<int,double> > & data) {
	if (!std::isdigit(line[0])) {
		// this is the line "variableStep chrom=chr2" or "fixedStep chrom=chr3 start=400601 step=100"
		// this indicates this is a new chr
		string delimit=" "; // a space
		vector<string> items;
		split(items, line, is_any_of(delimit));
		for (int i=0; i<items.size(); i++) {
			//cerr << "Item: " << items[i] << endl << flush;
			if (items[i].substr(0,6)=="chrom=") {
				chr = items[i].substr(6,items[i].size()-6);
				//cerr << chr << endl << flush;
				break;
			}
		}
		map<int,double> position2value;
		data[chr] = position2value; // create an empty map
	} else {
		// this is the line "chromStart  dataValue", which is delimited by a TAB
		string delimit="\t"; // a tab
		vector<string> items;
		split(items, line, is_any_of(delimit));
		//cerr << "Item 1: " << items[0] << endl << flush;
		//cerr << "Item 2: " << items[1] << endl << flush;
		int position = (int)atoi(items[0].c_str());
		double value = atof(items[1].c_str());
		//cerr << position << "\t" << value << endl << flush;
		map<int,double> & position2value = data[chr];
		position2value[position] = value;
	}
}

void read_wig_file(string wig_file, map<string,map<int,double> > & data) {
	//cerr << "reading '" << wig_file << "'" << endl;
	unsigned long i_line=0;
	string chr="";
	// input is a plain text file
	ifstream fin;
	fin.open(wig_file.c_str());
	if (fin.fail()){
		cerr << "Error: Unable to open " << wig_file << " in read_wig_file()" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of wig file
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		process_one_line_of_wig_file(line, chr, data);
		i_line++;
	}
	fin.close();
}

//void print_one_line_of_data(string chr, map<int,double> & position2value, stringstream & out) {
	//map<int,double>::iterator it2;
	//for (it2=position2value.begin(); it2!=position2value.end(); ++it2) {
		//out << chr << "\t" << it2->first << "\t" << it2->second << endl;
	//}
//}

// "data" is wig data
// print wig data to another format:
// Column 1: chr
// Column 2: position
// Column 3: value
void print_wig_data(map<string,map<int,double> > & data) {
	cout.precision(15);
	map<string,map<int,double> >::iterator it1;
	for (it1=data.begin(); it1!=data.end(); ++it1) {
		string chr = it1->first;
		map<int,double> & position2value = it1->second;
		map<int,double>::iterator it2;
		for (it2=position2value.begin(); it2!=position2value.end(); ++it2) {
			cout << chr << "\t" << it2->first << "\t" << it2->second << endl;
		}
	}
}

// add a position,value to wig bins data
void add_a_position_value_to_wig_bins_data(Wig_bins_data & wig_bins_data,
	int bin_internal_index, string chr, int position, double value)
{
	if ( wig_bins_data.find(chr) == wig_bins_data.end() ) {
		// this is a new chr
		map<int,map<int,double> > info_of_a_chr_;
		wig_bins_data[chr] = info_of_a_chr_;
	}
	map<int,map<int,double> > & info_of_a_chr = wig_bins_data[chr];
	if ( info_of_a_chr.find(bin_internal_index) == info_of_a_chr.end() ) {
		// this is a new bin
		map<int,double> info_of_a_bin_;
		info_of_a_chr[bin_internal_index] = info_of_a_bin_;
	}
	map<int,double> & info_of_a_bin = info_of_a_chr[bin_internal_index];
	info_of_a_bin[position] = value;
}

// "data" is wig data
// binning all wig data to "wig_bins_data"
void binning_wig_data(map<string,map<int,double> > & data,
	Bins_end_coord & bins_end_coord, Bins_index & bins_index,
	Wig_bins_data & wig_bins_data)
{
	map<string,map<int,double> >::iterator it1;
	for (it1=data.begin(); it1!=data.end(); ++it1) {
		string chr = it1->first;
		map<int,double> & position2value = it1->second;
		map<int,double>::iterator it2;
		for (it2=position2value.begin(); it2!=position2value.end(); ++it2) {
			int position = it2->first;
			double value = it2->second;
			int bin_internal_index = find_bin_of_position(bins_end_coord, chr, position);
			if (bin_internal_index==-1) {
				// Error: it is wierd this position doesn't belong to any bin. This is impossible. Must debug data or codes.
				// This can hapen, if chr of this wig data doesn't exist in genome bins annotation "bins_end_coord".
				//cerr << "Warn: " << chr << ":" << position << ":" << value << " doesn't exist in bin!" << endl; 
				continue;
			}
			int marker_index = bins_index[chr][bin_internal_index];
			if (marker_index!=0) {
				// this bin is a marker, not a complementary bin
				// record this position and its value
				add_a_position_value_to_wig_bins_data(wig_bins_data, bin_internal_index, chr, position, value);
			}
		}
	}
}

// "wig_of_a_bin" is a map of "position -> value".
// We print this map to one line: a list of position/value pairs, separated by a space; each pair is separated by ":"
void print_wig_of_a_bin(map<int,double> & wig_of_a_bin, ostream & out) {
	map<int,double>::iterator it;
	size_t size_ = wig_of_a_bin.size();
	out << size_ << "\t";
	size_t i = 1;
	for (it=wig_of_a_bin.begin(); it!=wig_of_a_bin.end(); ++it) {
		if (i<size_)
			out << it->first << ":" << it->second << " ";
		else
			out << it->first << ":" << it->second << endl;
		++i;
	}

}

// "data" is wig data
// print wig data to marker bins:
//
//Each line is a marker bin (or feature). All columns are delimited by TAB. There is one header line.
//Column 1: marker_index, 1-based index. Those bins that are not markers, are indexed as 0.
//Column 2: chr
//Column 3: start coordinate of bin (1-base)
//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
//Column 5: marker_type. "I" is Marker_type_I, "II" is Marker_type_II, "-" is the complementary bin, only facilitating the searching.
//Column 6: a list of position/value pairs, separated by a space; each pair is separated by ":"
//
// The following is an example file
void print_wig_bins_data(Wig_bins_data & wig_bins_data, Bins_index & bins_index,
	Bins_info & bins_info)
{
	vector<string> all_chrs;
	make_complete_chromosomes( all_chrs );
	cout.precision(15);
	vector<string>::iterator chr_it;
	for (chr_it=all_chrs.begin(); chr_it!=all_chrs.end(); ++chr_it) {
		// print wig data in the right order of chromosomes, which is as the order of chromosomes in "all_chrs".
		string chr = *chr_it;
		vector<int> & markers_index_of_a_chr = bins_index[chr];
		vector<string> & markers_info_of_a_chr = bins_info[chr];
		if (wig_bins_data.find(chr) == wig_bins_data.end()) {
			// this chr doesn't exist in wig_bins_data
			// We still output marker bin's info, but no wig data
			for (int bin_internal_index=0; bin_internal_index<markers_index_of_a_chr.size(); ++bin_internal_index) {
				if (markers_index_of_a_chr[bin_internal_index]!=0) {
					// this bin is a marker; complementary bin's index is zero
					cout << markers_info_of_a_chr[bin_internal_index] << "\t"; // print the first five columns
					cout << "0" << "\t" << "-" << endl;
				}
			}
		} else {
			map<int,map<int,double> > & wig_data_of_a_chr = wig_bins_data[chr];
			for (int bin_internal_index=0; bin_internal_index<markers_index_of_a_chr.size(); ++bin_internal_index) {
				if (markers_index_of_a_chr[bin_internal_index]!=0) {
					// this bin is a marker; complementary bin's index is zero
					cout << markers_info_of_a_chr[bin_internal_index] << "\t"; // print the first five columns
					if (wig_data_of_a_chr.find(bin_internal_index) == wig_data_of_a_chr.end()) {
						// this bin doesn't exist in wig_bins_data
						cout << "0" << "\t" << "-" << endl;
					} else {
						map<int,double> & wig_of_a_bin = wig_data_of_a_chr[bin_internal_index];
						print_wig_of_a_bin(wig_of_a_bin, cout); // print number of position/value pairs and the list of these pairs.
					}
				}
			}
		}
	}
}


//file of a single value (of each tissue) for each marker (generated by program "build_features_bins.py")
//This file should have the same number lines (bins) as the number of markers (excluding complementary bins) in bins annotation file.
//Each line is a marker bin (or feature). All columns are delimited by TAB. There is one header line.
//The first 5 columns are the same as File 1.
//Column 1: marker_index, 1-based index. Only marker bins are included, those complementary bins do not appear in this file.
//Column 2: chr
//Column 3: start coordinate of bin (1-base)
//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
//Column 5: marker_type. "I" is Marker_type_I, "II" is Marker_type_II, "-" is the complementary bin, only facilitating the searching.
//Column 6+: values for this marker (each column is a tissue)

//For example:
//marker_index	chr	start	end	marker_type	Liver	T-cells	B-cells	Neutrophils	Placenta
//1	chr1	855266	855766	II	0.472	0.848	0.607	0.298	0.144
//2	chr1	969796	970296	II	0.659	0.458	0.352	0.332	0.331
//3	chr1	1099044	1099544	II	0.722	0.879	0.928	0.636	0.134
//4	chr1	1109315	1109815	II	0.308	0.25	0.179	0.24	0.152
//...
//5816	chr22	50962109	50962609	I	0.912	0.894	0.935	0.952	0.888
//5817	chr22	50987071	50987571	II	0.104	0.037	0.038	0.621	0.029
//5818	chr22	51016754	51017254	II	0.321	0.487	0.63	0.298	0.153
//5819	chr22	51041995	51042495	II	0.129	0.116	0.051	0.429	0.051
//5820	chr22	51136171	51136671	II	0.531	0.898	0.926	0.337	0.1
//
// input file is TAB-delimited plain text. The file has a header line
// value_column_start_index: the column index of the first value (index is 1-base). In this format, it should be 6.
// bins2values: a map for bin_index -> a vector of value. bin_index is always 1-base. In detail, it is "map<unsigned int, vector<double> >"
void read_value_file_of_bins(string value_file_of_bins,
	unsigned int value_column_start_index, Bins2Values & bins2values,
	vector<string> & value_names)
{
	//cerr << "reading '" << value_file_of_bins << "'" << endl;
	unsigned long i_line=0;
	// input is a plain text file
	ifstream fin;
	fin.open(value_file_of_bins.c_str());
	if (fin.fail()){
		cerr << "Error: Unable to open " << value_file_of_bins << " in read_value_file_of_bins()" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of value_file_of_bins
			vector<string> items;
			split(items, line, is_any_of("\t"));
			for (int i=value_column_start_index-1; i<items.size(); i++) {
				value_names.push_back( items[i] );
			}
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		vector<string> items;
		split(items, line, is_any_of("\t"));
		if (value_column_start_index>items.size()) {
			cerr << "Error(read_value_file_of_bins): the value column starts from Column " << value_column_start_index << " that is > total number of columns (" << items.size() << ") in Line " << i_line+1 << "!" << endl;
			exit(EXIT_FAILURE);
		}
		int marker_index = (int)atoi(items[0].c_str());
		vector<double> values;
		for (int i=value_column_start_index-1; i<items.size(); i++) {
			values.push_back( (double)atof(items[i].c_str()) );
		}
		bins2values[marker_index] = values;
		i_line++;
	}
	fin.close();
}


//file of a single value (of each tissue) for each marker (generated by program "build_features_bins.py")
//This file should have the same number lines (bins) as the number of markers (excluding complementary bins) in bins annotation file.
//Each line is a marker bin (or feature). All columns are delimited by TAB. There is one header line.
//The first 5 columns are the same as File 1.
//Column 1: marker_index, 1-based index. Only marker bins are included, those complementary bins do not appear in this file.
//Column 2: chr
//Column 3: start coordinate of bin (1-base)
//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
//Column 5: marker_type. "-" is the complementary bin, only facilitating the searching.
//Column 6+: paired values for this marker (each column is a tissue). Each pair contains two values, delimited by ":"

//For example:
//marker_index	chr	start	end	marker_type	normal_plsa	tumor_brca
//1	chr1	855266	855766	brca	47.2,28.9	8.48,10.7
//2	chr1	969796	970296	brca	6.59,10.9	5.58,8.9
//3	chr1	1099044	1099544	brca	7.22,20.8	8.79,3.2
//4	chr1	1109315	1109815	brca	30.8,20.7	25.8,78.9
//...
//5816	chr22	50962109	50962609	I	91.2,10.7	20.7,60.7
//5817	chr22	50987071	50987571	II	1.04,10.2	15.7,1.8
//5818	chr22	51016754	51017254	II	32.1,10.7	10.8,2.9
//5819	chr22	51041995	51042495	II	1.29,1.9	5.7,8.9
//5820	chr22	51136171	51136671	II	53.1,10.7	3.0,7.1
//
// input file is TAB-delimited plain text. The file has a header line
// value_column_start_index: the column index of the first value (index is 1-base). In this format, it should be 6.
// bins2pairedvalues: a map for bin_index -> a pair of two vectors of values. bin_index is always 1-base. In detail, it is "map<int, pair<vector<double>,vector<double> > >"
void read_paired_values_file_of_bins(string value_file_of_bins,
	unsigned int value_column_start_index, Bins2PairedValues & bins2pairedvalues,
	vector<string> & value_names)
{
	//cerr << "reading '" << value_file_of_bins << "'" << endl;
	unsigned long i_line=0;
	// input is a plain text file
	ifstream fin;
	fin.open(value_file_of_bins.c_str());
	if (fin.fail()){
		cerr << "Error: Unable to open " << value_file_of_bins << " in read_paired_values_file_of_bins()" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of value_file_of_bins
			vector<string> items;
			split(items, line, is_any_of("\t"));
			for (int i=value_column_start_index-1; i<items.size(); i++) {
				value_names.push_back( items[i] );
			}
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		vector<string> items;
		split(items, line, is_any_of("\t"));
		if (value_column_start_index>items.size()) {
			cerr << "Error(read_paired_values_file_of_bins): the value column starts from Column " << value_column_start_index << " that is > total number of columns (" << items.size() << ") in Line " << i_line+1 << "!" << endl;
			exit(EXIT_FAILURE);
		}
		int marker_index = (int)atoi(items[0].c_str());
		vector<double> values1, values2;
		for (int i=value_column_start_index-1; i<items.size(); i++) {
			vector<string> subitems;
			split(subitems, items[i], is_any_of(",:"));
			if (subitems.size()!=2) {
				cerr << "Error(read_paired_values_file_of_bins): the value column " << (i+1) << " has no paired values (only " << subitems.size() << " values) in Line " << i_line+1 << "!" << endl;
				exit(EXIT_FAILURE);
			}
			values1.push_back( (double)atof(subitems[0].c_str()) );
			values2.push_back( (double)atof(subitems[1].c_str()) );
		}
		bins2pairedvalues[marker_index] = make_pair(values1, values2);
		i_line++;
	}
	fin.close();
}

//file of position-specific value of a tissue for each marker (generated by program "wig_binning")
//This file should have the same number lines (bins) as the number of markers (excluding complementary bins) in bins annotation file.
//Each line is a marker bin (or feature). All columns are delimited by TAB. There is one header line.
//The first 5 columns are the same as File 1.
//Column 1: marker_index, 1-based index. Those bins that are not markers, are indexed as 0.
//Column 2: chr
//Column 3: start coordinate of bin (1-base)
//Column 4: end coordinate of bin (1-base). The range of the bin is [start, end)
//Column 5: marker_type. "I" is Marker_type_I, "II" is Marker_type_II, "-" is the complementary bin, only facilitating the searching.
//Column 6: number of CpG sites in this marker region
//Column 7: a list of position/value pairs, separated by a space; each pair is separated by ":"

//For example:
//1	chr1	855266	855766	II	44	855298:0.8125 855299:0.8125 855404:0.933333333333333 855405:0.933333333333333 855409:0.8 855410:0.8 855425:0.75 855426:0.75 855428:0.818181818181818 855429:0.818181818181818 855437:1 855438:1 855445:1 855446:1 855472:0.833333333333333 855473:0.833333333333333 855478:0.857142857142857 855479:0.857142857142857 855486:0.909090909090909 855487:0.909090909090909 855544:0.266666666666667 855545:0.266666666666667 855549:0.133333333333333 855550:0.133333333333333 855567:0.571428571428571 855568:0.571428571428571 855609:0.1 855610:0.1 855643:0.5 855644:0.5 855649:0.777777777777778 855650:0.777777777777778 855654:0.454545454545455 855655:0.454545454545455 855657:0.333333333333333 855658:0.333333333333333 855664:0.545454545454545 855665:0.545454545454545 855692:0.0666666666666667 855693:0.0666666666666667 855708:0.357142857142857 855709:0.357142857142857 855741:0.363636363636364 855742:0.363636363636364
//2	chr1	969796	970296	II	64	969825:0.2 969826:0.2 969828:0.333333333333333 969829:0.333333333333333 969852:0.2 969853:0.2 969949:0.5 969950:0.5 969955:0.5 969956:0.5 969966:0.6 969967:0.6 969973:0.6 969974:0.6 969980:0.6 969981:0.6 969996:0.333333333333333 969997:0.333333333333333 970019:0.5 970020:0.5 970043:0.833333333333333 970044:0.833333333333333 970048:0.5 970049:0.5 970052:0.6 970053:0.6 970058:0.8 970059:0.8 970060:0.8 970061:0.8 970093:0.4 970094:0.4 970095:0.2 970096:0.2 970099:0.6 970100:0.6 970104:0.666666666666667 970105:0.666666666666667 970107:0.833333333333333 970108:0.833333333333333 970129:0.2 970130:0.2 970158:0.4 970159:0.4 970160:0.2 970161:0.2 970163:0.4 970164:0.4 970169:0.8 970170:0.8 970171:0.5 970172:0.5 970180:0.75 970181:0.75 970191:0.8 970192:0.8 970214:0 970215:0 970223:0.285714285714286 970224:0.285714285714286 970237:0.625 970238:0.625 970294:0.6 970295:0.6
//
// input file is plain text
// Bins2PositionValuePairs: a map for bin_index -> a pair of (position vector, value vector). bin_index is always 1-base. In detail, it is "map<unsigned int, pair<vector<GENOME_POSITION>,vector<double> > >".
void read_position_value_pairs_file_of_bins(string position_value_pairs_file_of_bins,
	Bins2PositionValuePairs bins2positionvaluepairs)
{
	//cerr << "reading '" << position_value_pairs_file_of_bins << "'" << endl;
}

//
//input file of indexing reads by bins (or features). (generated by program "reads_binning")
//Each line is a read. All columns are delimited by TAB. There is one header line.
//Column 1: marker1_index, 1-based index. "0" indicates the bin "marker1_index" is not a marker.
//Column 2: marker2_index, 1-based index. "-" indicates the read completely falls into "marker1"; "0" indicates the bin "marker2_index" is not a marker.
//Column 3: chr
//Column 4: start coordinate of read (1-base)
//Column 5: end coordinate of read (1-base). The range of the read is [start, end)
//Column 6: strand (+ or -)
//Column 7: number of CpG sites
//Column 8: list of CpG coordinates observed in this read (delimited by a comma). For example, 10469,10471,10484,10489,10493,10497,10525,10542
//Column 9: a binary vector of methylation status for all CpG sites in this read (no delimit). This vector should have the same size as list size in Column 6. For example, 00111100
//
//For example:
//1	-	chr1	855441	855515	+	4	855445,855472,855478,855486	1111
//1	-	chr1	855442	855516	+	4	855445,855472,855478,855486	1111
//1	-	chr1	855518	855577	-	3	855544,855549,855567	111
//1	-	chr1	855536	855608	-	3	855544,855549,855567	110
//0	2	chr1	969782	969856	+	3	969825,969828,969852	101
//2	-	chr1	969824	969897	-	4	969825,969828,969852,969883	1110
//2	-	chr1	969844	969918	-	4	969852,969883,969898,969914	1000
//...
//
// cpg_sites_column_index is 1-base
// methy_status_column_index is 1-base
void process_one_line_of_reads_binning_file(string & line, int cpg_sites_column_index, int methy_status_column_index, int & marker_index, vector<GENOME_POSITION> & cpg_sites, vector<int> & methy_status)
{
	vector<string> items;
	split(items, line, is_any_of("\t"));
	marker_index = (unsigned int)atoi(items[0].c_str());
	if (marker_index==0) {
		// indicates this read covers two bins: a complementary bin (marker_index=0) and a marker bin (marker_index is non-zero).
		// so we should take the second item as the marker_index
		// for example:
		// 0	2	chr1	969782	969856	+	3	969825,969828,969852	101
		marker_index = (unsigned int)atoi(items[1].c_str());
	}
	vector<string> subitems;
	split(subitems, items[cpg_sites_column_index-1], is_any_of(","));
	for (int i=0; i<subitems.size(); i++)
		cpg_sites.push_back((GENOME_POSITION)atoi(subitems[i].c_str()));
	string & str_methy_status = items[methy_status_column_index-1];
	for (int i=0; i<str_methy_status.size(); i++)
		methy_status.push_back(str_methy_status.at(i)-'0'); // convert each char '0' or '1' to 0 or 1, respectively.
}

double calc_one_read_prob_by_single_value(double value, vector<int> & methy_status)
{
	if (!methy_status.empty()) {
		double p=1.0;
		for (int i=0; i<methy_status.size(); i++) {
			// for each CpG site
			if (methy_status[i]==1) {
				// methylation status is 1 (methylated)
				p *= value;
			} else {
				// methylation status is 0 (unmethylated)
				p *= 1-value;
			}
		}
		return p;
	} else {
		return -1;
	}
}

//
//input file of indexing reads by bins (or features). (generated by program "reads_binning")
//Each line is a read. All columns are delimited by TAB. There is one header line.
//
// bins2values: a map for bin_index -> a vector of value. bin_index is always 1-base. In detail, it is "map<unsigned int, vector<double> >"
//
// Function: calculate the probability of each read that belongs to a tissue
//
// Output: file of probability value for each read
//
//Each line is a read. All columns are delimited by TAB. There is one header line.
//Column 1: marker index (1-base)
//Columns 2+: each column is a probability value of a specific tissue. For example, given 15 tissues, there will be 15 additional columns.
void calc_read_probability_by_bins2values(string reads_binning_file, 
	Bins2Values & bins2values, vector<string> & value_names)
{
	int column_index_of_cpg_sites = 8; // 1-base index
	int column_index_of_methy_status = 1 + column_index_of_cpg_sites; // 0-base index
	cout.precision(15);
	//cerr << "reading '" << reads_binning_file << "'" << endl;
	// print out the header line
	cout << "marker_index";
	for (int i=0; i<value_names.size(); i++)
		cout << "\t" << value_names[i];
	cout << "\n";
	// process each read (i.e., line) of the input file, then print out the calculated probabilities for each read.
	unsigned long i_line=0; // line count of the input file
	// input is a plain text file
	istream * in=&cin; // default is stdin
	ifstream fin;
	if ( !(reads_binning_file.compare("stdin")==0)) {
		fin.open(reads_binning_file.c_str());
		if (fin.fail()){
			cerr << "Error: Unable to open " << reads_binning_file << " in calc_read_probability_by_bins2values()" << endl;
			exit(EXIT_FAILURE);
		}
		in = &fin;
	}
	string line;
	while (!(*in).eof()) {
		getline((*in), line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of wig file
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		int marker_index;
		vector<GENOME_POSITION> cpg_sites;
		vector<int> methy_status;
		process_one_line_of_reads_binning_file(line, column_index_of_cpg_sites, column_index_of_methy_status, marker_index, cpg_sites, methy_status);
		if (!(bins2values.find(marker_index) == bins2values.end())) {
			//cerr << "marker " << marker_index << " does not exist in markers list" << endl;
		//} else {
			cout << marker_index; 
			vector<double> & values = bins2values[marker_index];
			for (int t=0; t<values.size(); t++) {
				double p = calc_one_read_prob_by_single_value(values[t], methy_status);
				cout << "\t" << p; 
			}
			cout << "\n";
		}
		i_line++;
	}
	if ( !(reads_binning_file.compare("stdin")==0) ) {
		fin.close();
	}
}

//
//input file of indexing reads by bins (or features). (generated by program "reads_binning")
//Each line is a read. All columns are delimited by TAB. There is one header line.
//
// bins2values: a map for bin_index -> a vector of value. bin_index is always 1-base. In detail, it is "map<unsigned int, vector<double> >"
//
// Function: calculate the probability of each read that belongs to a tissue
//
// Output: file of probability value for each read
//
//Each line is a read. All columns are delimited by TAB. There is one header line.
//Column 1: marker index (1-base)
//Columns 2+: each column is a probability value of a specific tissue. For example, given 15 tissues, there will be 15 additional columns.
void calc_read_probability_by_methy_level_bins(string reads_binning_file, int num_bins_of_methylation_level)
{
	vector<double> values;
	vector<string> value_names;
	double mvalue_step = 1.0/(double)num_bins_of_methylation_level;
	char str[64];
	for (int i=0; i<num_bins_of_methylation_level; i++) {
		double v = i*mvalue_step + mvalue_step/2;
		values.push_back( v );
		sprintf(str,"%.6g", i*mvalue_step);
		string start = string( str );
		sprintf(str,"%.6g", (i+1)*mvalue_step);
		string end = string( str );
		sprintf(str,"%d", i+1 );
		string idx = string( str );
		sprintf(str,"%.6g", v );
		string v_str = string( str );
		value_names.push_back( idx + ":(" + start + "," + end + "):" + v_str );
	}

	int column_index_of_cpg_sites = 8; // 1-base index
	int column_index_of_methy_status = 1 + column_index_of_cpg_sites; // 0-base index
	cout.precision(15);
	//cerr << "reading '" << reads_binning_file << "'" << endl;
	// print out the header line
	cout << "marker_index";
	for (int i=0; i<value_names.size(); i++)
		cout << "\t" << value_names[i];
	cout << "\n";

	// process each read (i.e., line) of the input file, then print out the calculated probabilities for each read.
	unsigned long i_line=0; // line count of the input file

	// input is a plain text file
	ifstream fin;
	fin.open(reads_binning_file.c_str());
	if (fin.fail()){
		cerr << "Error: Unable to open " << reads_binning_file << " in calc_read_probability_by_methy_level_bins()" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of wig file
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		int marker_index;
		vector<GENOME_POSITION> cpg_sites;
		vector<int> methy_status;
		process_one_line_of_reads_binning_file(line, column_index_of_cpg_sites, column_index_of_methy_status, marker_index, cpg_sites, methy_status);
		cout << marker_index; 
		for (int t=0; t<values.size(); t++) {
			double p = calc_one_read_prob_by_single_value(values[t], methy_status);
			cout << "\t" << p; 
		}
			cout << "\n";
			i_line++;
	}
	fin.close();
}


// Function: calculate the probability of each read that belongs to a tissue
//    the formula is B(x+alpha, 1-x+beta)/B(alpha,beta)
//    where B() is beta function, x is methylation status of a CpG site (0 or 1), alpha and beta are a pair of values for describing a Beta-distribution provided by the input parameter "bins2pairedvalues".
//
// background_B: a precomputed B(alpha,beta)
//
double calc_one_read_prob_by_a_paired_value(double alpha, double beta, double background_B, vector<int> & methy_status)
{
	if (!methy_status.empty()) {
		double p=1.0;
		for (int i=0; i<methy_status.size(); i++) {
			// for each CpG site
			//cerr << "gammaln_of_all: " << ( gsl_sf_lngamma(methy_status[i]+alpha) + gsl_sf_lngamma(1-methy_status[i]+beta) - gsl_sf_lngamma(alpha+beta+1) - background_B ) << endl << flush;
			p *= gsl_sf_exp( gsl_sf_lngamma(methy_status[i]+alpha) + gsl_sf_lngamma(1-methy_status[i]+beta) - gsl_sf_lngamma(alpha+beta+1) - background_B );
		}
		return p;
	} else {
		return -1;
	}
}
//
//input file of indexing reads by bins (or features). (generated by program "reads_binning")
//Each line is a read. All columns are delimited by TAB. There is one header line.
//
// bins2pairedvalues: a map for bin_index -> a pair of two vectors of value. bin_index is always 1-base. In detail, it is "map<unsigned int, pair<vector<double>,vector<double> > >"
//
// Function: calculate the probability of each read that belongs to a tissue
//    the formula is B(x+alpha, 1-x+beta)/B(alpha,beta)
//    where B() is beta function, x is methylation status of a CpG site (0 or 1), alpha and beta are a pair of values for describing a Beta-distribution provided by the input parameter "bins2pairedvalues".
//
// Output: file of probability value for each read
//
//Each line is a read. All columns are delimited by TAB. There is one header line.
//Column 1: marker index (1-base)
//Columns 2+: each column is a probability value of a specific tissue. For example, given 15 tissues, there will be 15 additional columns.
void calc_read_probability_by_bins2pairedvalues(string reads_binning_file,
	Bins2PairedValues & bins2pairedvalues, vector<string> & value_names)
{
	int column_index_of_cpg_sites = 8; // 1-base index
	int column_index_of_methy_status = 1 + column_index_of_cpg_sites; // 0-base index
	// precompute the B(alpha,beta) for each bin
	map<int,vector<double> > B;
	Bins2PairedValues::iterator it;
	for (it=bins2pairedvalues.begin(); it!=bins2pairedvalues.end(); ++it) {
		int bin_index = it->first;
		vector<double> & alphas = (it->second).first;
		vector<double> & betas = (it->second).second;
		vector<double> values;
		int n = alphas.size();
		for (int t=0; t<n; t++) {
			//values.push_back( gsl_sf_beta(alphas[t], betas[t]) );
			values.push_back( gsl_sf_lngamma(alphas[t]) + gsl_sf_lngamma(betas[t]) - gsl_sf_lngamma(alphas[t]+betas[t]) );
			//cerr << "a: " << alphas[t] << ", b: " << betas[t] << "gammaln(z)+gammaln(w)-gammaln(z+w): " << gsl_sf_lngamma(alphas[t]) + gsl_sf_lngamma(betas[t]) - gsl_sf_lngamma(alphas[t]+betas[t]) << endl;
		}
		B[bin_index] = values;
	}
	cout.precision(15);
	//cerr << "reading '" << reads_binning_file << "'" << endl;
	// print out the header line
	cout << "marker_index";
	for (int i=0; i<value_names.size(); i++)
		cout << "\t" << value_names[i];
	cout << "\n";
	// process each read (i.e., line) of the input file, then print out the calculated probabilities for each read.
	unsigned long i_line=0; // line count of the input file
	// input is a plain text file
	istream * in=&cin; // default is stdin
	ifstream fin;
	if ( !(reads_binning_file.compare("stdin")==0)) {
		fin.open(reads_binning_file.c_str());
		if (fin.fail()){
			cerr << "Error: Unable to open " << reads_binning_file << " in calc_read_probability_by_bins2values()" << endl;
			exit(EXIT_FAILURE);
		}
		in = &fin;
	}
	string line;
	while (!(*in).eof()) {
		getline((*in), line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of wig file
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		int marker_index;
		vector<GENOME_POSITION> cpg_sites;
		vector<int> methy_status;
		process_one_line_of_reads_binning_file(line, column_index_of_cpg_sites, column_index_of_methy_status, marker_index, cpg_sites, methy_status);
		if (!(bins2pairedvalues.find(marker_index) == bins2pairedvalues.end())) {
			//cerr << "marker " << marker_index << " does not exist in markers list" << endl;
		//} else {
			cout << marker_index; 
			vector<double> & alphas = bins2pairedvalues[marker_index].first;
			vector<double> & betas = bins2pairedvalues[marker_index].second;
			for (int t=0; t<alphas.size(); t++) {
				double p = calc_one_read_prob_by_a_paired_value(alphas[t], betas[t], B[marker_index][t], methy_status);
				cout << "\t" << p; 
			}
			cout << "\n";
		}
		i_line++;
	}
	if ( !(reads_binning_file.compare("stdin")==0) ) {
		fin.close();
	}
}
// Function: calculate the probability of each read that belongs to a tissue
//
// Output: two vectors
//    (1) number of methylated sites of each read 
//    (2) number of all sites of each read
//
void get_reads_methy_data_from_reads_binning_file(string reads_binning_file,
	vector<int> & num_cpg_sites, vector<int> & num_methy_cpg_sites)
{
	int column_index_of_cpg_sites = 8; // 1-base index
	int column_index_of_methy_status = 1 + column_index_of_cpg_sites; // 0-base index
	cout.precision(15);
	//cerr << "reading '" << reads_binning_file << "'" << endl;
	// print out the header line
	cout << "marker_index" << "\t" << "num_CpG_sites" << "\t" << "num_methy_CpG_sites" << endl;
	// process each read (i.e., line) of the input file, then print out the numbers of CpG sites and methylated CpG sites for each read.
	unsigned long i_line=0; // line count of the input file
	// input is a plain text file
	ifstream fin;
	fin.open(reads_binning_file.c_str());
	if (fin.fail()){
		cerr << "Error: Unable to open " << reads_binning_file << " in get_reads_methy_data_from_reads_binning_file()" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (!fin.eof()) {
		getline(fin, line);
		//cerr << "Line: " << line << endl;
		if (i_line==0) {
			// skip the header line of wig file
			i_line++;
			continue;
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		int marker_index;
		vector<GENOME_POSITION> cpg_sites;
		vector<int> methy_status;
		process_one_line_of_reads_binning_file(line, column_index_of_cpg_sites, column_index_of_methy_status, marker_index, cpg_sites, methy_status);
		int num_methy_sites=0;
		for (int t=0; t<methy_status.size(); t++)
			if (methy_status[t]==1)
				num_methy_sites++;
		num_methy_cpg_sites.push_back(num_methy_sites);
		num_cpg_sites.push_back(methy_status.size());
		//cout << marker_index << "\t" << methy_status.size() << "\t" << num_methy_sites << endl; 
		i_line++;
	}
	fin.close();
}

// print vector of uint to the format, e.g., 10469,10471,10484,10489,10493,10497,10525,10542
void print_vec_of_uint(ostream& of, vector<unsigned int> & v) {
	switch (v.size()) {
		case 0:
			return;
		case 1:
			of << v[0];
			break;
		default: // v.size()>=2
			int i;
			for (i=0; i<v.size()-1; i++) of << v[i] << ",";
			of << v[i];
	}
}

// print vector of ulong to the format, e.g., 10469,10471,10484,10489,10493,10497,10525,10542
void print_vec_of_ulong(ostream& of, vector<unsigned long> & v) {
	switch (v.size()) {
		case 0:
			return;
		case 1:
			of << v[0];
			break;
		default: // v.size()>=2
			int i;
			for (i=0; i<v.size()-1; i++) of << v[i] << ",";
			of << v[i];
	}
}

void print_map_of_strings(ostream& of, map<string, vector<string> > & map_str1tostrs2) {
	map<string, vector<string> >::iterator it;
	for (it=map_str1tostrs2.begin(); it!=map_str1tostrs2.end(); it++) {
		string str1 = it->first;
		vector<string> & strs2 = it->second;
		for (int i=0; i<strs2.size(); i++)
			of << str1 << "\t" << strs2[i] << endl;
	}
}

void print_str_vectors(ofstream& o, vector<string> & s)
{
	vector<string>::iterator it;
	for (it=s.begin(); it!=s.end(); ++it)
		o << *it << endl;
}

void strings2floats(string str, vector<float> & vec, string delimit) {
	vector<string> strs;
	split(strs, str, is_any_of(delimit));
	if (!strs.empty())
		for (int i=0; i<strs.size(); i++)
			vec.push_back((float)atof(strs[i].c_str()));
}

float min(float a, float b) {
	return a < b ? a : b;
}

float max(float a, float b) {
	return a > b ? a : b;
}

float meta_p(vector<float> & probs) {
	if (!probs.empty()) {
		double x = 0;
		int n_effective = 0;
		for (int i=0; i<probs.size(); i++)
			if (probs[i]!=-1) {
				x += -2*log(probs[i]);
				n_effective++;
			}
		int n_degree = 2*n_effective;
		if (n_effective==0) return -1;
		else return gsl_cdf_chisq_Q(x, n_degree);
	} else
		return -1;
}
