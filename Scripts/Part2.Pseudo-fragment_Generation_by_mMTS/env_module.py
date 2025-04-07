"""
Author: Yin Liang
Date: 20250407
define the file path of proccess
"""
#dir
current_dir = "/home/yinliang/PROJECT/DeepLB/"
root_out_dir = f"{current_dir}/Result/2.simulation_result/"
sample_list_dir = f"{current_dir}/Result/1.4_gen_rand_runs/"
marker_dir = f"{current_dir}/Result/1.6_select_marker/"
tumor_beta_dir = f"{current_dir}/Result/1.3_cal_methylratio_bins/450K/"
wgbs_bam_dir = "/bigData/yinliang/Project/WGBS_data/"
# diff cohort test data dir
backgroud_cohort = "dismir_sort_bam_all/"
cohort1 = "CRA001537/"
cohort2 = "CRA001142/"
cohort3 = "inhouse_sort_bam/"


#f"{current_dir}/Predata/WGBS/"
#file
all_samples_file = f"{current_dir}/Predata/metadata/all_samples_annotation.txt"
sample_metadata = f"{current_dir}/Predata/metadata/all_WGBS_sample_metadata.xlsx"
refcpg = f"{current_dir}/Predata/reference/cpg_system/hg38_CG_noheader_full.CpGsite"
ref_fastq = f"{current_dir}/Predata/reference/GRCh38/GRCh38.fa"
pcc = f"{current_dir}/Result/1.2_define_bins/blocks_D200_S3.txt"
mcta_list = f"{current_dir}/Predata/metadata/MCTA_sample.txt"
#folder name
preprocess_folder = "/2-0_process/"
#including:DMR,sample_list,sample_beta,DMR_cpgsites,DMR_ref_fasta
ref_cpg_beta_dir = "/2-1_ref_cpg_beta/"
#root_dir + env_module.ref_cpg_beta_dir + marker_type + "/" + str(i ) + "_" + threshold + "/"
ref_reads_dir = "/2-2_ref_reads/"
#root_dir + env_module.ref_reads_dir + marker_type + "/" + str(i ) + "_" + threshold + "/"
simulated_reads_dir = "/2-3_simulated_reads/"
#root_dir + env_module.simulated_reads_dir + marker_type + "/" +  str(i ) + "_" + threshold + "/"
real_bam_dir = "/2-5_real_process_bam/"
#root_dir + env_module.real_bam_dir + "train/" + marker_type + "/" +  str(i ) + "_" + threshold + "/"
reads_to_train_dir = "/2-4_reads_to_train/"
reads_to_test_dir = "/2-6_reads_to_test/"
mcta_reads_dir = "/2-7_test_MCTA_reads/"
ajust_reads_dir = "/2-8_ajust_reads/"
ajust_reads_to_train_dir = "/2-9_ajust_to_train/"