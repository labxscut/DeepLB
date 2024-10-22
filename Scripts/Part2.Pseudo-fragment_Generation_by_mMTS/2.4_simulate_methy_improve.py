"""
Author: Yin Liang
Date: 10/16/2024
"""
import random
import os
import pandas as pd
import time
import argparse
import env_module
from collections import defaultdict
import numpy as np

#Reading the methylation file for CpG sites
def read_beta_file(beta_file):
    beta_data = {}
    with open(beta_file, "r") as file:
        for line in file:
            fields = line.strip().split()
            pos = f"{fields[0].split()[0]},{fields[1]},{fields[2]}"
            beta_data[pos] = fields[3]
    return beta_data

# reads the reads file and stores the location and base sequence
def read_reads_file(reads_file):
    reads_data = {}
    with open(reads_file, "r") as file:
        for line in file:
            fields = line.strip().split()
            pos = f"{fields[0]} {fields[1]} {fields[2]}"
            seq = fields[3]
            methy = fields[4]
            reads_data[pos] = [seq, methy]
    return reads_data

def initialize_methylation_dict(reads_data):
    reads_methy_dict = {}
    for index, sequence in reads_data.items():
        methylation_vector = sequence[1]
        seq = sequence[0]
        reads_methy_dict[index] = {seq: methylation_vector}
    return reads_methy_dict

def read_bed_file(bed_file):
    dmr_list = []
    with open(bed_file, "r") as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 3:
                dmr_list.append((parts[0], int(parts[1]), int(parts[2])))
    return dmr_list

def index_beta_data(beta_data):
    beta_index = defaultdict(list)
    for cpg_key, beta_value in beta_data.items():
        if beta_value == "nan":
            continue
        cpg_chr, cpg_start, cpg_end = cpg_key.split(",")
        beta_index[cpg_chr].append((int(cpg_start), int(cpg_end), beta_value))
    return beta_index

def extract_dmr_beta_data(dmr_list, beta_index):
    dmr_beta_data = defaultdict(list)
    for dmr in dmr_list:
        dmr_chr, dmr_start, dmr_end = dmr
        if dmr_chr in beta_index:
            for cpg_start, cpg_end, beta_value in beta_index[dmr_chr]:
                if dmr_start <= cpg_start <= dmr_end:
                    dmr_beta_data[dmr].append((cpg_start, cpg_end, beta_value))
    return dmr_beta_data


def reform_DMR_reads(reads_methy_dict,dmr_start,dmr_end):
    chr_positions  = []
    dna_bases_matrix  = []
    methylation_matrix  = []
    segment_length = dmr_end - dmr_start + 1

    for key, base_methylation in reads_methy_dict.items():
        read_chr, read_start, read_end = key.split()
        start, end = int(read_start), int(read_end)
        chr_positions.append([read_chr, start])
        dna_sequence  = list(base_methylation.keys())[0]
        methylation_sequence = list(base_methylation.values())[0]

        bases_row = ['o'] * segment_length
        methylation_row = [3] * segment_length

        for i, (base, methylation) in enumerate(zip(dna_sequence, methylation_sequence)):
            relative_position = start + i - dmr_start
            if 0 <= relative_position < segment_length:
                bases_row[relative_position] = base
                methylation_row[relative_position] = int(methylation)
        
        dna_bases_matrix.append(bases_row)
        methylation_matrix.append(methylation_row)
    
    chr_positions_df = pd.DataFrame(chr_positions, columns=['Chromosome', 'Start_Position'])
    dna_bases_matrix = np.array(dna_bases_matrix)
    methylation_matrix = np.array(methylation_matrix)
    
    return chr_positions_df, dna_bases_matrix, methylation_matrix

def process_dmr_regions(dmr_list, reads_methy_dict, region_beta):
    pos_mat = []
    dna_base_mat = []
    dna_methy_mat = []
    #region_beta = Region_beta_data

    for dmr in dmr_list:
        dmr_reads_methy_dict = {}
        dmr_chr, dmr_start, dmr_end = dmr
        cpg_list = region_beta[dmr]
        for read_key, read_data in reads_methy_dict.items():
            read_chr, read_start, read_end = read_key.split()
            if read_chr == dmr_chr and dmr_start <= int(read_start) <= dmr_end:
                dmr_reads_methy_dict[read_key] = read_data
                
        pos_info,DNA_info,Methy_info = reform_DMR_reads(dmr_reads_methy_dict,dmr_start,dmr_end)

        cpg_index = []
        for cpg_info in cpg_list:
            cpg_start = cpg_info[0]
            cpg_end = cpg_info[1]
            cpg_beta = float(cpg_info[2])
            #print(cpg_bata)
            if cpg_beta == "nan":
                continue
            else:
                cpg_rela_pos = cpg_start - dmr_start - 1
                cpg_index.append(cpg_rela_pos)


        for index in cpg_index:
            #print(index)
            cover_index = np.where(Methy_info[:, index] != 3)[0]
            #print(f"cover_index :{list(cover_index)}")
            num_methylated_reads = round(
                        float(cpg_beta) * len(cover_index)
                    )
            #print(f"num_methylated_reads: {num_methylated_reads}")
            if num_methylated_reads > 0:
                indices = random.sample(
                    list(cover_index),num_methylated_reads
                )
                num = 0 
                #print(indices)
                for ind in indices:
                    if DNA_info[ind,index] == "C" and DNA_info[ind,index+1] == "G":
                        Methy_info[ind,index] = "1"
                        num = num + 1
                    elif DNA_info[ind,index] == "C" and DNA_info[ind,index+1] == "o":
                        Methy_info[ind,index] = "1"
                        num = num + 1
                    else:
                        print(DNA_info[ind,index])
                        print(DNA_info[ind,index+1])
                #print(num)
        
        pos_mat.append(pos_info)
        dna_base_mat.append(DNA_info)
        dna_methy_mat.append(Methy_info)
    return(pos_mat,dna_base_mat,dna_methy_mat)

def write_output_file(output_file, pos_mat,dna_base_mat,dna_methy_mat):
    with open(output_file, "w") as file:
        for dmr in list(range(len(pos_mat))):
            for ind in list(range(len(pos_mat[dmr]))):
                chr = pos_mat[dmr].iloc[ind,0]
                start = pos_mat[dmr].iloc[ind,1]
                read = dna_base_mat[dmr][ind][dna_base_mat[dmr][ind] != "o"]
                methy = dna_methy_mat[dmr][ind][dna_methy_mat[dmr][ind] != 3]
                read_chr = str(chr)
                read_start = str(start)
                read_seq = ''.join(read)
                df_methy = methy.astype(str)
                read_methy = ''.join(df_methy)
                file.write(f"{read_chr}\t{read_start}\t{read_seq}\t{read_methy}\n")
def parse_args():
    parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
    parser.add_argument('-t', '--tumor',  help='Add the tumor type') #"lihc" "paad" "stad" "brca"
    parser.add_argument('-s', '--threshold', help='Subtract the numbers') #"0.15" "0.25"
    parser.add_argument('-g', '--group', help='Divide the numbers')#
    parser.add_argument('-a', '--approach', help='Divide the numbers') #method :"paied" or "tumoronly"
    parser.add_argument('-r', '--rep',  help='Divide the numbers') #total list num
    parser.add_argument('-l', '--cohort', help='Divide the numbers') #list
    parser.add_argument('-m', '--marker',  help='Marker type')#"hyper" "hypo"
    parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity')#top30 top40 top50 bottom30 bottom40 bottom50
    return parser.parse_args()

TUMOR_PREFIX = {
    "lihc": "LIHC",
    "paad": "PAAD",
    "stad": "STAD",
    "brca": "BRCA"
}

def main():
    args = parse_args()
    purity = args.purity if args.purity else ""
    start_time = time.time()
    threshold = args.threshold
    method = args.approach + "/" + purity + "/"
    tumor_type = args.tumor 
    marker_type = args.marker
    rep =  int(args.rep)
    cohort = int(args.cohort)
    group = args.group
    #simulated sample prefix
    if not tumor_type:
        raise ValueError("Tumor type must be specified.")
    prefix = TUMOR_PREFIX.get(tumor_type)
    root_dir =  env_module.root_out_dir + method + "/" + tumor_type + "-" + group + "/"
    sample_list_dir =  root_dir  + env_module.preprocess_folder +"sample_list/"

    for i in range(cohort,rep):
        sample_list = sample_list_dir + "train_tumor_" + str(i) + ".list"
        samples = pd.read_csv(sample_list ,delimiter='\t', header=None)
        samples.columns = ["sample"]
        sample_read_dir = root_dir + env_module.ref_reads_dir + marker_type + "/" + str(i ) + "_" + threshold + "/"
        bed_dir = root_dir  + env_module.preprocess_folder + "DMR/"
        bed_file = bed_dir + str(i) + "_" + threshold + "_" + marker_type + ".bed"
        beta_dir = root_dir + env_module.ref_cpg_beta_dir + marker_type + "/" + str(i ) + "_" + threshold + "/"
        output_dir = root_dir + env_module.simulated_reads_dir + marker_type + "/" +  str(i ) + "_" + threshold + "/"
        os.makedirs(output_dir , exist_ok=True)
        print(output_dir)
        print("-----Process : " + tumor_type + "  " + marker_type + "  " +str(i) + "_"+ threshold +"-----")
        print(" Complete : " + str(i-1) + "/" + str(rep))

        for s in samples["sample"]:
            sample_start = time.time()
            print("→processed Sample:" + str(s))
            beta_file = beta_dir + prefix +"_"+ str(s)+ ".refcpg.beta"
            reads_file = sample_read_dir  + prefix +"_"+  str(s) + ".ref.reads"
            output_file = output_dir + prefix +"_" +  str(s)+ ".methy.reads"
            beta_data = read_beta_file(beta_file)
            reads_data = read_reads_file(reads_file)
            reads_methy_dict = initialize_methylation_dict(reads_data)
            dmr_list = read_bed_file(bed_file)
            beta_index = index_beta_data(beta_data)
            Region_beta_data = extract_dmr_beta_data(dmr_list, beta_index)
            POS,DNA,METHY = process_dmr_regions(dmr_list, reads_methy_dict, Region_beta_data)
            write_output_file(output_file,POS,DNA,METHY)
            sample_end = time.time()
            sample_times = sample_end - sample_start
            print(f"sample {prefix}_{str(s)} executed in {sample_times:.2f} seconds")
            print("⇒Complete sample: " +str(s) )
        print("==============================")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Script executed in {execution_time:.2f} seconds")


if __name__ == "__main__":
    main()
