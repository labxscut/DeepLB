"""
Author: Yin Liang
Date: 10/16/2024
"""

#input 
import os
import pandas as pd
from tqdm import tqdm
import argparse
import env_module

def process_cpg_sites(input_file, beta_file, output_file):
    # 读取 BED 文件，构建区间字典
    bed_intervals = {}
    with open(beta_file, "r") as bed:
        for line in bed:
            fields = line.strip().split()
            chrom, start, end, value = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                float(fields[3]),
            )
            bed_intervals[(chrom, start, end)] = value

    # 处理输入文件中的 CpG 位点
    with open(input_file, "r") as cpg_file, open(output_file, "w") as output:
        for line in cpg_file:
            fields = line.strip().split()
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])

            # 检查 CpG 位点是否在任何区间内

            for interval, value in bed_intervals.items():
                if chrom == interval[0] and start >= interval[1] and end <= interval[2]:
                    fields.append(str(value))
                    output.write(
                        f"{fields[0]}\t{str(fields[1])}\t{ str(fields[2])}\t{str(fields[4])}\n"
                    )

parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold',help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group',  help='Divide the numbers')#"allsample" " onlylow" "onlyhigh"
parser.add_argument('-a', '--approach', help='Divide the numbers') #method :"paied" or "tumoronly",默认paired不加任何标签
parser.add_argument('-r', '--rep', help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort', help='Divide the numbers') #list
parser.add_argument('-m', '--marker', help='Marker type')#"hyper" "hypo"
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 bottom30 bottom40 bottom50
args = parser.parse_args()

purity = args.purity
if purity is None:
    purity = ""
tumor_type = args.tumor 
group = args.group
method = args.approach + "/" + purity + "/"
threshold = args.threshold
rep =  int(args.rep)
cohort = int(args.cohort)
marker_type = args.marker

#simulated sample prefix
if tumor_type == "lihc":
    prefix = "LIHC"
elif tumor_type == "paad":
    prefix = "PAAD"
elif tumor_type == "stad":
    prefix = "STAD"
elif tumor_type == "brca":
    prefix = "BRCA"


parameter = [threshold,tumor_type,rep,marker_type,group,prefix]
print(f"Process parameter: {parameter}")
root_dir =  env_module.root_out_dir + method + "/" + tumor_type + "-" + group
sample_list_dir = root_dir  + env_module.preprocess_folder +"sample_list/"
input_dir = root_dir + env_module.preprocess_folder + "DMR_cpgsite/"
for i in tqdm(range(cohort, rep), desc="Processing files"):
    print(i)
    train_cohort = sample_list_dir + "train_tumor_" + str(i) + ".list"
    samples = pd.read_csv(train_cohort,delimiter='\t', header=None)
    samples.columns = ["sample"]
    print("samples num: "+str( len(samples)))
    beta_dir = root_dir + env_module.preprocess_folder + "sample_beta/" + "train_" + str(i) + "/"
    output_dir =root_dir + env_module.ref_cpg_beta_dir + marker_type + "/"+ str(i) + "_" +threshold + "/"
    os.makedirs(output_dir, exist_ok=True)
    for s in  tqdm(samples["sample"], desc=f"Processing files for iteration {i}", leave=False):
        print("→Process Sample: "+ str(s))
        input_file = input_dir + str(i) + "_" + threshold + "_" + marker_type +"_all.cpgsites"
        output_file = output_dir +  prefix + "_" +  str(s) + ".refcpg.beta"
        if os.path.exists(output_file) and os.stat(output_file).st_size != 0:
                print(f"{s} has been processed!")
        else:
            beta_file = beta_dir + prefix +"_" + str(s)+ ".aver.beta"
            process_cpg_sites(input_file, beta_file, output_file)
