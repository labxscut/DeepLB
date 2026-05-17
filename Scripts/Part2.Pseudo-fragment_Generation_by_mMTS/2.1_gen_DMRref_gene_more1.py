"""Extract DMR regions into reference fasta files.

Author: Yin Liang
Date: 10/16/2024
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from tqdm import tqdm
import argparse
import env_module


parser = argparse.ArgumentParser(description='Parameters for threshold/tumor type/marker type/group/rep')
parser.add_argument('-t', '--tumor', help='Add the tumor type') #"lihc" "paad" "stad" "brca"
parser.add_argument('-s', '--threshold',help='Subtract the numbers') #"0.15" "0.25"
parser.add_argument('-m', '--marker', help='Marker type')#"hyper" "hypo"
parser.add_argument('-g', '--group',  help='Divide the numbers')#"TH" "MH" "PH"
parser.add_argument('-a', '--approach',  help='Divide the numbers') # method: "paired" or "tumoronly"; default paired without extra tag
parser.add_argument('-r', '--rep',  help='Divide the numbers') #total list num
parser.add_argument('-l', '--cohort', help='Divide the numbers') #list
parser.add_argument('-p', '--purity',  nargs='?', default=None,help='for tumor purity')
args = parser.parse_args()

purity = args.purity
if purity is None:
    purity = ""
threshold = args.threshold
tumor_type = args.tumor
marker_type = args.marker
rep = int(args.rep)
cohort = int(args.cohort)
method = args.approach + "/" + purity + "/"
group = args.group
root_dir =  env_module.root_out_dir + method + "/" + tumor_type + "-" + group
input_fasta = env_module.ref_fastq
bed_dir = root_dir + env_module.preprocess_folder + "DMR/"
out_dir = root_dir + env_module.preprocess_folder + "DMR_ref_fasta/" + marker_type + "/"
os.makedirs(out_dir, exist_ok=True)

def extract_fasta_regions(input_fasta, input_bed, output_fasta):
    # Read reference genome fasta
    sequences = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

    # Prepare output fasta records
    output_records = []

    with open(input_bed, "r") as bed_file:
        for line in bed_file:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2]) + 1  # Add 1 bp for downstream CG site check
            region_name = f"{chrom}:{start}-{end}"  # Use BED chrom/start/end as record ID

            if chrom in sequences:
                sequence = sequences[chrom][start:end]
                record = SeqRecord(sequence.seq, id=region_name, description="")
                output_records.append(record)

    # Write output fasta
    SeqIO.write(output_records, output_fasta, "fasta")

for i in range(cohort, rep):
    bed_file = bed_dir + str(i) + "_" + threshold + "_" + marker_type + ".bed"
    out_file = out_dir + str(i) + "_" + threshold + ".ref.fa"
    extract_fasta_regions(input_fasta, bed_file, out_file)


