"""Extract reads and methylation strings from BAMs for train/test.

Author: Yin Liang
Date: 10/16/2024
"""
import pysam
import numpy as np
from collections import Counter
import time
import argparse
import os
import os.path
import env_module

def reverse_complement(read, methy):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}    
    reverse_read = read[::-1]
    reverse_methy = methy[::-1]
    complement_read = "".join(complement.get(base, base) for base in reverse_read)
    return complement_read, reverse_methy

def process_line(line):
    sig = line.cigarstring
    if "I" in sig or "N" in line.query_sequence:
        return None, None, None, None
    name = line.reference_name
    start = str(line.reference_start)
    end = str(line.reference_end)
    read = line.query_sequence
    flag = line.flag
    methy = line.get_tag("XM")
    return name, start, read, flag, methy

def process_methylation(flag, read, methy):
    if flag == 16:
        read, methy = reverse_complement(read, methy)
    methylation = ""
    for i in range(len(methy)):
        if methy[i] == "Z":
            if i < len(methy) - 1:
                if read[i] == "C" and read[i + 1] == "G":
                    methylation += "1"
                else:
                    methylation += "0"
            elif i == len(methy) - 1:
                if read[i] == "C":
                    methylation += "1"
                else:
                    methylation += "0"
            else:
                methylation += "0"
        else:
            methylation += "0"
            
        if methy[i] != ".":  ##T to C
                read = list(read)
                read[i] = "C" 
                read = "".join(read)

    return read, methylation

def input_sample(i):
    test = sample_list_dir + "test_" + str(i) + ".list"
    train = sample_list_dir + "train_healthy_" + str(i) + ".list"
    test_list = []
    train_list = []  
    with open(test, "r") as samples1:
        for t in samples1:
            test_list.append(t)
    with open(train, "r") as samples2:
        for t in samples2:
            train_list.append(t)
    return test_list,train_list

def process_files(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    total_replacements = 0 

    for filename in os.listdir(input_folder):
        if filename.endswith('.reads'):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)

            with open(input_path, 'r') as input_file, open(output_path, 'w') as output_file:
                replacements_in_file = 0  

                for line in input_file:
                    fields = line.strip().split('\t')
                    if len(fields) == 4 :
                        fields[3] = fields[3].replace('2', '0')
                        updated_line = '\t'.join(fields)
                        if fields[0] == "chrX":
                            continue
                        output_file.write(updated_line + '\n')
                        replacements_in_file += line.count('2')

                total_replacements += replacements_in_file

                print(f"File '{filename}': replaced {replacements_in_file} occurrences of 2")

    print(f"Total replacements: {total_replacements} occurrences of 2")

parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-m', '--marker',  help='Marker type')#"hyper" "hypo"
parser.add_argument('-s', '--threshold', help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group',  help='Divide the numbers')#"allsample" " onlylow" "onlyhigh"
parser.add_argument('-a', '--approach', help='Divide the numbers') # method: "paired" or "tumoronly"; default paired without extra tag
parser.add_argument('-r', '--rep',help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort',  help='Divide the numbers') #list
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 bottom30 bottom40 bottom50
args = parser.parse_args()

purity = args.purity
if purity is None:
    purity = ""
threshold = args.threshold
tumor_type = args.tumor 
marker_type = args.marker
rep =  int(args.rep)
cohort = int(args.cohort)
group = args.group
method = args.approach + "/" + purity + "/"
Group = ["train","test"]
root_dir = env_module.root_out_dir + method + "/" +  tumor_type + "-" + group + "/"
sample_list_dir = root_dir  + env_module.preprocess_folder +"sample_list/"

for i in range(cohort,rep):
    start_time = time.time()
    print("-----Process : " + tumor_type + "  " + marker_type + "  " +str(i) + "_"+ threshold +"-----")
    test,train_healthy = input_sample(i)
    for g in Group : 
        bam_dir = root_dir + env_module.real_bam_dir + g + "/" +marker_type + "/" +  str(i ) + "_" + threshold + "/"
        
        if g == "train":
            samples = train_healthy
            out_dir = root_dir + env_module.reads_to_train_dir + "/healthy/" + marker_type + "/" +  str(i ) + "_" + threshold + "/"   
        else:
            samples = test
            out_dir = root_dir + env_module.reads_to_test_dir + marker_type + "/" +  str(i ) + "_" + threshold + "/"   
        
        os.makedirs(out_dir , exist_ok=True)
        print(bam_dir)
        print(out_dir)
        miss_sample = []

        for s in samples:
            bam_path = bam_dir + s.strip() + "_process.bam"
            reads_file = out_dir  + s.strip() + "_bam.reads"
            if os.path.exists(reads_file) :
                print(f"{ s.strip()} exist")
            else:
                if os.path.exists(bam_path) :
                    with pysam.AlignmentFile(bam_path, "rb") as bam_file, open(reads_file, "w") as output_file:
                        for line in bam_file:
                            result = process_line(line)
                            if result[0] is None:
                                continue
                            name, start, read, flag, methy = result
                            read2, methylation = process_methylation(flag, read, methy)
                            if len(methylation) < 71 or name == "chrX"  or name == "chrY":
                                continue
                            output_file.write(name+'\t'+start+'\t'+read2[5:71]+'\t'+methylation[5:71]+'\n')
                else:
                    print(f"{s.strip()} does not exist. Skipping.")
                    miss_sample.append(s.strip())
    
    print(miss_sample)
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Complete!! Script executed in {execution_time:.2f} seconds")
