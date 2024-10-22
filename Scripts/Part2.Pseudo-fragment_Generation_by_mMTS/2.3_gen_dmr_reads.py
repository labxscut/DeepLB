"""
Author: Yin Liang
Date: 10/16/2024
"""
import random
from Bio import SeqIO
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import env_module
# parameters
parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor',  help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-m', '--marker', help='Marker type')#"hyper" "hypo"
parser.add_argument('-s', '--threshold',  help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group', help='Divide the numbers')#"allsample" " onlylow" "onlyhigh"
parser.add_argument('-a', '--approach',  help='Divide the numbers') #method :"paied" or "tumoronly",默认paired不加任何标签
parser.add_argument('-r', '--rep',  help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort',  help='Divide the numbers') #list
parser.add_argument('-L', '--Length',  help='Divide the numbers') #simulated parameter 3
parser.add_argument('-D', '--Depth',  help='Divide the numbers') #simulated parameter 66
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 bottom30 bottom40 bottom50
args = parser.parse_args()


purity = args.purity
if purity is None:
    purity = ""
threshold = args.threshold
method = args.approach + "/" + purity + "/"
tumor_type = args.tumor 
marker_type = args.marker
rep =  int(args.rep)
cohort = int(args.cohort)
D = int(args.Depth)
L = int(args.Length)
group = args.group
#simulated sample prefix
if tumor_type == "lihc":
    prefix = "LIHC"
elif tumor_type == "paad":
    prefix = "PAAD"
elif tumor_type == "stad":
    prefix = "STAD"
elif tumor_type == "brca":
    prefix = "BRCA"

root_dir =  env_module.root_out_dir  + method + "/" + tumor_type + "-" + group
ref_dir = root_dir  + env_module.preprocess_folder + "DMR_ref_fasta/"+ marker_type + "/" 
sample_list_dir = root_dir  + env_module.preprocess_folder +"sample_list/"

def generate_read(fasta_records, output_file):
    fragments = {}
    for key, record in fasta_records.items():
        header = record.description
        sequence = record.seq
        parts = header.split(":")
        chrom = parts[0]
        pos = parts[1].rstrip("a").split("-")
        start = int(pos[0])
        end = int(pos[1])  #
        R = end - start
        # print(start,end,sequence)
        S = int(D * R / L)
        positions = random.sample(range(R - L - 1), S)
        # 原：positions = [random.randint(0, 500 - L -1 ) for _ in range(S)]
        for position in positions:
            fi_start = start + position
            fi_end = fi_start + L
            fragment_pos = f"{chrom} {fi_start} {fi_end}"
            position = int(position)
            fragment_sequence = sequence[position : position + L + 1]
            cg_positions = []
            for i in range(len(fragment_sequence) - 1):
                if fragment_sequence[i] == "C" and fragment_sequence[i + 1] == "G":
                    cg_positions.append(i)
            fragment_methy = np.zeros(len(fragment_sequence), dtype=int)
            if len(cg_positions) != 0:
                for pos in cg_positions:
                    fragment_methy[pos] = 2
            else:
                pass
            fragments[fragment_pos] = [fragment_sequence[:-1], fragment_methy[:-1]]

    # 将生成的片段写入输出文件
    with open(output_file, "w") as output_file:
        for key, value in fragments.items():
            methy = "".join([str(x) for x in value[1]])
            output_line = f"{key} {str(value[0])} {methy}\n"
            output_file.write(output_line)


for i in tqdm(range(cohort, rep), desc="Processing files"):
    input = ref_dir + str(i) + "_" + threshold + ".ref.fa"
    output_dir = root_dir + env_module.ref_reads_dir + marker_type + "/"+ str(i) + "_" + threshold + "/"
    os.makedirs(output_dir,exist_ok=True)
    fasta_records = SeqIO.to_dict(SeqIO.parse(input, "fasta"))
    sample_list = sample_list_dir + "train_tumor_" + str(i) + ".list"
    samples = pd.read_csv(sample_list ,delimiter='\t', header=None)
    samples.columns = ["sample"]
    print("samples num: "+str( len(samples)))
    for x in samples["sample"]:
        print("→Process Samples: " + str(x))
        out_file = output_dir + prefix +"_" + str(x) + ".ref.reads"
        generate_read(fasta_records, out_file)
