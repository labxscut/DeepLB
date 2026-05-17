"""Select healthy BAM reads overlapping DMRs for training.

Author: Yin Liang
Date: 10/16/2024
"""
import os
import numpy as np
import time
import argparse
import env_module
#parameters
parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-m', '--marker',  help='Marker type')#"hyper" "hypo"
parser.add_argument('-s', '--threshold',  help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group',  help='Divide the numbers')#
parser.add_argument('-a', '--approach',  help='Divide the numbers') # method: "paired" or "tumoronly"; default paired without extra tag
parser.add_argument('-r', '--rep',  help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort',  help='Divide the numbers') #list
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity')#top30 top40 top50 bottom30 bottom40 bottom50
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

root_dir =  env_module.root_out_dir + method + "/" + tumor_type + "-" + group + "/"
sample_list_dir = root_dir  + env_module.preprocess_folder +"sample_list/"

bed_dir = root_dir  + env_module.preprocess_folder + "DMR/"
wgbs_dir = env_module.wgbs_bam_dir

def input_sample(i):
    sample_list = sample_list_dir + "train_healthy_" + str(i) + ".list"
    samples = []
    with open(sample_list, "r") as SAMPLE:
        for s in SAMPLE:
                samples.append(s)
    return samples

for i in range(cohort,rep):
    start_time = time.time()
    print("-----Process : " + tumor_type + "  " + marker_type + "  " +str(i) + "_"+ threshold +"-----")
    bed_name = bed_dir + str(i) + "_" + threshold + "_" + marker_type + ".bed" #or hypo
    out_dir =  root_dir + env_module.real_bam_dir + "train/" + marker_type + "/" +  str(i ) + "_" + threshold + "/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    samples = input_sample(i)
    for s in samples:
        bam_dir = wgbs_dir + "/" + env_module.backgroud_cohort #+ "/normal/"  #plasma train data dir
        sample_file = bam_dir + s.strip() + ".bam"
        out_name = out_dir + s.strip() + "_process.bam"
        os.system('samtools view -b  -h ' + sample_file + ' -L "' + bed_name + '" >' + out_name )
        print("Complete sample: " + s.strip() )
        print("——————————————")
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Script executed in {execution_time:.2f} seconds")







