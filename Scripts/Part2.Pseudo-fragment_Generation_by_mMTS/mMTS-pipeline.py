"""
a standard pipeline for generating Pseudo-fragment by mMTS
Author: Yin Liang
Date: 10/16/2024
"""

import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor
import time
from tqdm import tqdm
import os
import env_module

parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep/skip steps')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold', help='Marker Selection marker') #"0.15" "0.25"
parser.add_argument('-m', '--marker', help='Marker type')#"hyper" "hypo"
parser.add_argument('-g', '--group', help='train healthy sample type')#"TH" "MH" "PH"
parser.add_argument('-a', '--approach', help='train sample method') #method :"paied" or "tumorOnly"
parser.add_argument('-r', '--rep',  help='total list num') #total list num
parser.add_argument('-l', '--cohort', help='list') #list first
parser.add_argument('-L', '--Length', help='simulated parameter 3') #simulated parameter 66
parser.add_argument('-D', '--Depth', help='simulated parameter 66') #simulated parameter 3
parser.add_argument('-e','--values', nargs='+', help='Except some steps',default="") #except some steps
parser.add_argument('-c', '--choose', help='choose diff method to find marker') #freq-diff mean-diff
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity') #top30 top40 top50 
args = parser.parse_args()

threshold = args.threshold
tumor_type = args.tumor
marker_type = args.marker
rep = args.rep
cohort = args.cohort
method = args.approach
group = args.group
D = args.Depth
L = args.Length
exceptstep = args.values
find_marker_method = args.choose
purity = args.purity
if purity is None:
    purity = ""

log_dir = env_module.root_out_dir + method + "/" + purity +  "/" + tumor_type + "-" + group +"/log/" + threshold + "_" + marker_type + "/"
os.makedirs(log_dir,exist_ok = True)
start_time = time.time()

#environment
#step_env = "conda activate base"
#print(step_env)
#try:
#    step_env_result = subprocess.Popen(step_env, shell=True)
#    step_env_result.wait()
#    print("→Environment activate successfully.")
#except subprocess.CalledProcessError as e:
#    print(f"Error executing command {step_env}:", e)


#step1-4 and step5-6 para processing,step7 and step8 is final 
#(0)prepare all needed file
#(1)prepare DMR ref gene
#(2)prepare simulated sample ref cpg beta
#(3)prepare DMR ref reads
#(4)simulated reads
#(5)select train healthy bam
#(6)select test bam
#(7)translate train tumor to reads
#(8)extract train healthy and test bam to reads

#(0)prepare all needed file

log0 = log_dir + "step0.log"
if method == "tumorOnly":
    step0 = f"python 2.0_prepare_files_tumorOnly.py -t {tumor_type} -s {threshold}  -r {rep} -l {cohort} -a {method} -g {group} -c {find_marker_method} -p {purity} > {log0}"
else:
    step0 = f"python 2.0_prepare_files_paired.py -t {tumor_type} -s {threshold}  -r {rep} -l {cohort} -a {method} -g {group} -c {find_marker_method} -p {purity} > {log0}"

    
#创建空的变量
step1 = ""
step2 = ""
step3 = ""
step4 = ""
step5 = ""
step6 = ""
step7 = ""
step8 = ""

steps = [
    ("2.1_gen_DMRref_gene_more1.py","step1"),
    ("2.2_generate_cpg_ref_beta.py", "step2"),
    ("2.3_gen_dmr_reads.py", "step3"),
    ("2.4_simulate_methy_improve.py", "step4"),
    ("2.5_select_healthy_bam_to_train.py", "step5"),
    ("2.6_select_bam_test_all.py", "step6"),
    ("2.7_trans_to_train_read.py", "step7"),
    ("2.8_extract_bam_to_reads.py", "step8"),
    ("2.9_process_MCTAdata.py", "step9")
]


for script, log_prefix in steps:
    log_file = f"{log_dir}{log_prefix}.log"
    if log_prefix in ["step3"]:
        step_command = f"python -u {script} -t {tumor_type} -s {threshold} -r {rep} -m {marker_type} -l {cohort} -a {method} -g {group} -D {D} -L {L} -p {purity} > {log_file}"
    else:
        step_command = f"python -u {script} -t {tumor_type} -s {threshold} -r {rep} -m {marker_type} -l {cohort} -a {method} -g {group} -p {purity} > {log_file}"
    exec(f"{log_prefix} = '{step_command}'") #创建变量

para1 = [step1,step2,step3,step4]
para2 = [step5,step6]
para3 = [step7,step8,step9]


#filter some steps
variables = {
    "step0": step0,
    "step1": step1,
    "step2": step2,
    "step3": step3,
    "step4": step4,
    "step5": step5,
    "step6": step6,
    "step7": step7,
    "step8": step8,
}

print("(0)all steps")
for key, value in variables.items():
    print(f'{key}: {value}')

#创建一个空列表，用于存储需要删除的变量
to_remove = []

#根据变量的值查找变量名称，并将需要删除的变量名称添加到列表中
for para in [para1, para2, para3]:
    for step in para:
        for var_name, var_value in variables.items():
            if var_value == step and var_name in exceptstep:
                to_remove.append(step)
print(f"Steps to remove {to_remove}")
# 从列表中删除需要删除的变量
for remove_step in to_remove:
    for para in [para1, para2, para3]:
        if remove_step in para:
            para.remove(remove_step)
print("############################")
print(f"para1 : {para1}")
print(f"para2 : {para2}")
print(f"para3 : {para3}")
print("         ")
def execute_commands(commands, desc):
    with tqdm(total=len(commands), desc=desc) as pbar:
        i = 0
        for command in commands:
            print("################################ Begin Process ########################################")
            print(f"({i + 1})Process → {command}")
            try:
                command_result = subprocess.Popen(command, shell=True)
                command_result.wait()
                print(f"!! executed successfully.")
                pbar.update(1)  # 更新进度条
            except subprocess.CalledProcessError as e:
                print(f"Error executing {command}: {e}")
                pbar.update(1)  # 更新进度条


#(0)prepare all needed file
try:
    if "step0" in exceptstep:
        print("(1)Skip Step0!!!!")
    else:
        print(f"Process → {step0}")
        step0_result = subprocess.Popen(step0, shell=True)
        step0_result.wait()
        print(f"!!executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error executing {step0}:", e)



# 创建线程池，同时执行para1和para2中的命令
with ThreadPoolExecutor() as executor:
    future1 = executor.submit(execute_commands, para1, "Executing para1")
    future2 = executor.submit(execute_commands, para2, "Executing para2")

    # 等待para1和para2中的命令执行完成
    future1.result()
    future2.result()

# 执行para3中的命令
execute_commands(para3, "Executing para3")

end_time = time.time()
execution_time = end_time - start_time
execution_time_minutes = execution_time/ 60
print(f"Script start at {start_time}, end at {end_time}")
print(f"Script executed in {execution_time_minutes:.2f} minutes")


