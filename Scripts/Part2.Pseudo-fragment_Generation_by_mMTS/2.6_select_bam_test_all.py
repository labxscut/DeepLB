"""Select test BAM reads overlapping DMRs for evaluation.

Author: Yin Liang
Date: 10/16/2024
"""
import os
import time
import argparse
import env_module

#parameters
parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor',  help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-m', '--marker',  help='Marker type')#"hyper" "hypo"
parser.add_argument('-s', '--threshold', help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group',help='Divide the numbers')#
parser.add_argument('-a', '--approach',  help='Divide the numbers') # method: "paired" or "tumoronly"; default paired without extra tag
parser.add_argument('-r', '--rep', help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort', help='Divide the numbers') #list
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

root_dir =  env_module.root_out_dir + method + "/" + tumor_type + "-" + group + "/"
sample_list_dir = root_dir  + env_module.preprocess_folder +"sample_list/"

bed_dir = root_dir  + env_module.preprocess_folder + "DMR/"
# test bam data dir:  set the dir according cohort
dir1 = env_module.wgbs_bam_dir + env_module.cohort1 # CRR05*
dir2 = env_module.wgbs_bam_dir + env_module.cohort2 # CRR03*
dir3 = env_module.wgbs_bam_dir + env_module.cohort3 # NC0*
dir4 = env_module.wgbs_bam_dir + env_module.backgroud_cohort #"TBR*/CTR*/HOT*"

#normal_dir = env_module.wgbs_bam_dir + "normal/"
#tumor_dir = env_module.wgbs_bam_dir + "tumor/"

def input_sample(i):
    sample_list = sample_list_dir + "test_" + str(i) + ".list"
    sample_name = []  
    with open(sample_list, "r") as samples:
        for t in samples:
            sample_name.append(t)
    return sample_name

for i in range(cohort,rep):
    start_time = time.time()
    print("-----Process: " + tumor_type  + "-" + group + marker_type + "  " +str(i) + "_"+ threshold +"-----")
    bed_name = bed_dir + str(i) + "_" + threshold + "_" + marker_type + ".bed" #or hypo
    out_dir = root_dir + env_module.real_bam_dir + "/test/" + marker_type + "/" +  str(i ) + "_" + threshold + "/"
    os.makedirs(out_dir , exist_ok=True)
    print(f"out dir : {out_dir}")

    sample_name = input_sample(i)
    print(f"sample num: {len(sample_name)}")
    for s in sample_name:
#        if "HOT" in s: #"tumor sample prex"
#            bam_dir = normal_dir
#        elif "CTR" in s :#"normal sample prex"
#            bam_dir = normal_dir
        if "HOT" in s or "TBR" in s or "CTR" in s or "HBV" in s:
            bam_dir = dir4
        elif "NC" in s or "GC"in s or "PC" in s:
            bam_dir = dir3
        elif "CRR05" in s:
            bam_dir = dir1
        elif "CRR03" in s:
            bam_dir = dir2
        sample_file = bam_dir + s.strip() + ".bam"       
        out_name = out_dir + s.strip() + "_process.bam"
        if os.path.exists(out_name) and os.stat(out_name).st_size != 0:
                print(f"{ s.strip()} has been processed!")
        else:
            if os.path.exists(sample_file):
                print(f"Processing {s.strip()}")
                os.system('samtools view -b  -h ' + sample_file + ' -L "' + bed_name + '" >' + out_name )
            else:
                print(f"{s.strip()} does not exist. Skipping.")

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Script executed in {execution_time:.2f} seconds")




    






