"""
Author: Yin Liang
Date: 10/16/2024
"""
import argparse
import os
import sys
import time
import env_module


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


parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor',  help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold',help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-g', '--group',  help='Divide the numbers')#
parser.add_argument('-a', '--approach',  help='Divide the numbers') #method :"paied" or "tumoronly"
parser.add_argument('-r', '--rep',  help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort',  help='Divide the numbers') #list
parser.add_argument('-m', '--marker', help='Marker type')#"hyper" "hypo"
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 bottom30 bottom40 bottom50
args = parser.parse_args()


purity = args.purity
if purity is None:
    purity = ""
threshold = args.threshold
method = args.approach + "/" + purity + "/"
rep =  int(args.rep)
cohort = int(args.cohort)
tumor_type = args.tumor 
marker_type = args.marker
group = args.group
root_dir = env_module.root_out_dir + method + "/" + tumor_type + "-" + group +"/"


for i in range(cohort,rep):
    start_time = time.time()
    print("-----Process : " + tumor_type + "  " + marker_type + "  " +str(i) + "_"+ threshold +"-----")
    input_folder = root_dir + env_module.simulated_reads_dir + marker_type + "/" +  str(i ) + "_" + threshold + "/"
    output_folder = root_dir +  env_module.reads_to_train_dir + "/tumor/"  + marker_type + "/" +  str(i ) + "_" + threshold + "/"
    process_files(input_folder, output_folder)

end_time = time.time()
execution_time = end_time - start_time
print(f"Script executed in {execution_time:.2f} seconds")
