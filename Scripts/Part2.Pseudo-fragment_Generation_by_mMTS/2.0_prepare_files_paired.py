"""
Author: Yin Liang
Date: 10/16/2024
"""

import os
import random
import re
import pandas as pd
from tqdm import tqdm
import argparse
import env_module
#run for lihc paired list1
# parameters


parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold', help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group', help='Divide the numbers')##"TH" "MH" "PH"
parser.add_argument('-a', '--approach', help='Divide the numbers') #method :"paied" or "tumoronly",默认paired不加任何标签
parser.add_argument('-r', '--rep', help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort', help='Divide the numbers') #list
parser.add_argument('-c', '--choose', help='Divide the numbers', default="freq-diff") #freq-diff mean-diff
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 
args = parser.parse_args()

if args.purity is None:
    purity = ""
else:
    purity = args.purity
tumor_type = args.tumor 
group = args.group
method = args.approach 
threshold = args.threshold
rep =  int(args.rep)
cohort = int(args.cohort)
marker_type = ["hypo","hyper"]
find_marker_method = args.choose
#simulated sample prefix
prefix = tumor_type.upper()


## 1、train cohort and test cohort
root_dir = env_module.root_out_dir + method + "/" + tumor_type + "-" + group
all_samples_file = env_module.all_samples_file
all_samples = pd.read_csv(all_samples_file, delimiter='\t', header=None)
all_samples.columns = ["sample","sample_type","type","group"]
meta_data_file = env_module.sample_metadata
meta_data = pd.read_excel(meta_data_file)
list_dir = env_module.sample_list_dir + method + "/" + tumor_type + "/" + group + "/"

meta_data["Group"].unique()
all_healthy_blood = meta_data.loc[(meta_data["Group"] == "Healthy") & (meta_data["Sample type"]== "cfDNA")]
for i in range(cohort,rep):
    train_list = list_dir + "train_" + str(i)  + ".list"
    print(train_list)
    train_samples  = pd.read_csv(train_list,delimiter='\t', header=None)
    train_tumor = train_samples.loc[~train_samples[1].isna(),0]
    train_healthy = train_samples.loc[train_samples[1].isna(),0]
    test_healthy = all_healthy_blood["Run"][~meta_data["Run"].isin(train_healthy)]
    #!!! need to change tumor sample index
    test_tumor = meta_data[meta_data['Group'].str.contains('breast cancer', case=False, na=False)]["Run"] 
    print(test_tumor)
    test_sample =  test_healthy.append(test_tumor).drop_duplicates()

    list_out_dir = root_dir  + env_module.preprocess_folder +"sample_list/"
    os.makedirs(list_out_dir,exist_ok=True)
    # each rep contain 3 sample list
    train_tumor_path = os.path.join(list_out_dir,f"train_tumor_{i}.list")
    train_tumor.to_csv(train_tumor_path, index=False, header=False)

    train_healthy_path = os.path.join(list_out_dir,f"train_healthy_{i}.list")
    train_healthy.to_csv(train_healthy_path,index=False, header=False)


    test_file_path = os.path.join(list_out_dir, f"test_{i}.list")
    test_sample.to_csv(test_file_path, index=False, header=False)
    print(f"Saved test list to: {list_out_dir}")


## 2、 prepare dmr.bed (all_markers,hypo,hyper)
print(f"------->Process {tumor_type} dmr of threshold about {threshold}")
marker_dir =  env_module.marker_dir + method+ "/" + find_marker_method  + "/"+ tumor_type + "/" + group +"/"
dmr_dir = root_dir  + env_module.preprocess_folder + "DMR/"
os.makedirs(dmr_dir,exist_ok=True)
for i in tqdm(range(cohort, rep), desc="Processing files"):
    marker_file = marker_dir + str(i) + "_" + threshold
    hypo_marker = dmr_dir + str(i)+ "_"+ threshold+ "_hypo.bed"
    hyper_marker = dmr_dir + str(i)+ "_"+ threshold+ "_hyper.bed"
    with open(marker_file, "r") as m, open(hypo_marker, "w") as hypo,open(hyper_marker,"w") as hyper :
        next(m)
        for line in tqdm(m, desc=f"Processing file {marker_file}", leave=False):
            info = line.rstrip().split("\t")
            if info[1] == "chrY" or info[1] == "chrX":#排除性染色体
                continue
            else:
                if info[4] == "hyper":
                    hyper.write(f"{info[1]}\t{info[2]}\t{info[3]}\n")
                if info[4] == "hypo":
                    hypo.write(f"{info[1]}\t{info[2]}\t{info[3]}\n")


## 3、 prepare DMR cpg all site
print("----------->prepare DMR cpg all site")
import subprocess
bed_dir = dmr_dir
allcpg_out_dir = root_dir + env_module.preprocess_folder + "DMR_cpgsite/"
os.makedirs(allcpg_out_dir,exist_ok=True)
refcpg = env_module.refcpg
for i in range(cohort,rep):
    for t in marker_type:
        bed_file = bed_dir + str(i) + "_" + threshold + "_" + t +".bed"
        cpg_outfile = allcpg_out_dir + str(i) + "_" + threshold + "_"+ t + "_all.cpgsites"
        command = "bedtools intersect -a " + refcpg + " -b " + bed_file + " > "  + cpg_outfile
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) 

## 4、prepare train tumor aver.beta
# prepare pcc annotation
#(1)paired
print("-------> prepare train tumor aver.beta")
def readPairedData(fileName1):
    pairedSampleList = []
    data1 = {}
    with open(fileName1) as f:
        for line in f:
            line = re.sub("NA", "nan", line)
            info = line.rstrip().split("\t")
            index = info.pop(0)
            data1[index] = info  # np.array(map(float,info))
            pairedSampleList.append([index])
        return (data1, pairedSampleList)
    
allbeta_file = env_module.tumor_beta_dir + tumor_type + "_tumor.tsv"
(allMethyTumor, pairedSampleList) = readPairedData(allbeta_file)
pcc_dir = env_module.pcc
pcc = []
with open(pcc_dir, "r") as p:
    next(p)
    for line in p:
        pcc.append(line)

for i in tqdm(range(cohort, rep), desc="Processing files"):
    train_cohort = list_out_dir + "train_tumor_" + str(i) + ".list"
    tumor_sample = pd.read_csv(train_cohort,delimiter='\t', header=None)
    beta_output_dir = root_dir + env_module.preprocess_folder + "sample_beta/" + "train_" + str(i) + "/"
    os.makedirs(beta_output_dir,exist_ok=True)
    tumor_sample.columns = ["sample"]
    for index in  tqdm(tumor_sample["sample"], desc=f"Processing files for iteration {i}", leave=False):
        beta_file = beta_output_dir + prefix + "_" + str(index) + ".aver.beta"
        betavalue = allMethyTumor[str(index)]
        with open(beta_file, "w") as beta:
            for x in range(0, len(pcc)):
                b = betavalue[x]
                pos = pcc[x].strip("\n").split()
                beta.write(f"{pos[1]}\t{pos[2]}\t{pos[3]}\t{b}\n")
