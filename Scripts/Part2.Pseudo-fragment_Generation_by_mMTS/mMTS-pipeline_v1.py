"""Run the mMTS pseudo-fragment pipeline with logging and step control.

a standard pipeline for generating Pseudo-fragment by mMTS
Author: Yin Liang
Date: 3/25/2025
"""

import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from tqdm import tqdm
import os
import logging
import env_module

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("pipeline.log"), logging.StreamHandler()]
)
logger = logging.getLogger()

parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep/skip steps')
parser.add_argument('-t', '--tumor', help='Add the tumor type')  # "lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold', help='Marker Selection marker')  # "0.15" "0.25"
parser.add_argument('-m', '--marker', help='Marker type')  # "hyper" "hypo"
parser.add_argument('-g', '--group', help='train healthy sample type')  # "TH" "MH" "PH"
parser.add_argument('-a', '--approach', help='train sample method')  # method :"paied" or "tumorOnly"
parser.add_argument('-r', '--rep', help='total list num')  # total list num
parser.add_argument('-l', '--cohort', help='list')  # list first
parser.add_argument('-L', '--Length', help='simulated parameter 3')  # simulated parameter 66
parser.add_argument('-D', '--Depth', help='simulated parameter 66')  # simulated parameter 3
parser.add_argument('-e', '--values', nargs='+', help='Except some steps', default=[])  # except some steps
parser.add_argument('-c', '--choose', help='choose diff method to find marker')  # freq-diff mean-diff
parser.add_argument('-p', '--purity', nargs='?', default=None, help='for tumor purity')  # top30 top40 top50
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
exceptstep = set(args.values)  # Convert to set for faster lookup
find_marker_method = args.choose
purity = args.purity if args.purity else ""

# Create log directory
log_dir = os.path.join(env_module.root_out_dir, method, purity, f"{tumor_type}-{group}", "log", f"{threshold}_{marker_type}")
os.makedirs(log_dir, exist_ok=True)

start_time = time.time()

# Define step commands
steps = [
    ("2.0_prepare_files_tumorOnly.py" if method == "tumorOnly" else "2.0_prepare_files_paired.py", "step0"),
    ("2.1_gen_DMRref_gene_more1.py", "step1"),
    ("2.2_generate_cpg_ref_beta.py", "step2"),
    ("2.3_gen_dmr_reads.py", "step3"),
    ("2.4_simulate_methy_improve.py", "step4"),
    ("2.5_select_healthy_bam_to_train.py", "step5"),
    ("2.6_select_bam_test_all.py", "step6"),
    ("2.7_trans_to_train_read.py", "step7"),
    ("2.8_extract_bam_to_reads.py", "step8")
]

commands = {}
for script, log_prefix in steps:
    log_file = os.path.join(log_dir, f"{log_prefix}.log")
    if log_prefix == "step3":
        command = f"python -u {script} -t {tumor_type} -s {threshold} -r {rep} -m {marker_type} -l {cohort} -a {method} -g {group} -D {D} -L {L} -p {purity} > {log_file}"
    else:
        command = f"python -u {script} -t {tumor_type} -s {threshold} -r {rep} -m {marker_type} -l {cohort} -a {method} -g {group} -p {purity} > {log_file}"
    commands[log_prefix] = command

# Filter steps to skip
filtered_commands = {k: v for k, v in commands.items() if k not in exceptstep}

# Execute commands
def execute_command(command, desc):
    try:
        logger.info(f"Executing: {desc} → {command}")
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"Successfully executed: {desc}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing {desc}: {e.stderr.decode('utf-8')}")
        raise

def execute_commands(commands, desc):
    with tqdm(total=len(commands), desc=desc) as pbar:
        for command in commands:
            try:
                execute_command(command, desc)
                pbar.update(1)
            except Exception as e:
                logger.error(f"Error in {desc}: {e}")
                pbar.update(1)

# Execute step0
if "step0" not in exceptstep:
    execute_command(commands["step0"], "Step0")
else:
    logger.info("Skipping Step0")

# Run para1 and para2 in parallel
para1 = [commands[k] for k in ["step1", "step2", "step3", "step4"] if k in filtered_commands]
para2 = [commands[k] for k in ["step5", "step6"] if k in filtered_commands]

with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [
        executor.submit(execute_commands, para1, "Executing para1"),
        executor.submit(execute_commands, para2, "Executing para2")
    ]
    for future in as_completed(futures):
        future.result()

# Run para3 sequentially
para3 = [commands[k] for k in ["step7", "step8"] if k in filtered_commands]
execute_commands(para3, "Executing para3")

end_time = time.time()
execution_time = end_time - start_time
execution_time_minutes = execution_time / 60
logger.info(f"Script executed in {execution_time_minutes:.2f} minutes")