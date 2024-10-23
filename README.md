# DeepBL
## Overview
In this study, we introduce our novel deep learning model, DeepLB, developed for early cancer detection through liquid biopsy. This model utilizes pseudo-fragment data generated from 450K methylation array data of tumor tissues, alongside cfDNA whole-genome bisulfite sequencing (WGBS) data from healthy donors, to enhance diagnostic capabilities. DeepLB consists of four key modules: feature selection, methylation Microarray-To-Sequence data converter (mMTS), deep learning model (ResTran), and cancer score estimation. We analyzed cancer-specific features from tumor tissue and healthy plasma WGBS data, generating tumor pseudo-fragment data through the mMTS converter. The ResTran model, which integrates residual networks with Transformer architecture, accurately identifies tumor-derived DNA fragments from cfDNA methylation sequencing data. The workflow in DeepLB is illustrated in the figure below.

![Alt text](png1.png)

## Installation
```R
install.packages('devtools')
devtools::install_github("labxscut/DeepLB")
```
```
../DeepLB/
├── env
├── codes # from CancerLocator and CancerDetector
├── log # for process log
├── Predata # The data should be prepare before Process DeepLB
├── Result 
└── Scripts
    ├── Part1.Marker_Selection
    ├── Part2.Pseudo-fragment_Generation_by_mMTS
    └──  Part3.ResTran_model_training
```

## Prepare data
!!!!! The data should be prepare before Process DeepLB
```
Predata/
├── 450K #Download from GEO and TCGA
│   ├── illuminaMethyl450_hg38_GDC 
│   ├── TCGA_450_probe 
│   ├── TCGA-LIHC.methylation450.tsv #tumor 450K data
│   └── TCGA_Study_Abbreviations.txt
├── metadata
│   ├── all_samples_annotation.txt #tumor tissue and normal plasma annotation
│   ├── sample_list.normal # normal plasma list
│   └── sample_metadata.txt # plasma sample metadata
├── reference
│   ├── cpg_system
│   │   ├── chrom_list.txt
│   │   ├── hg38_CG_noheader_full.CpGsite
│   ├── GRCh38
│   │   ├── Bisulfite_Genome
│   │   ├── GRCh38.fa
│   │   ├── GRCh38.fa.fai
│   │   ├── hg38.chrom.sizes
│   │   ├── index.html
│   │   └── sequence.fa
│   ├── hg38_chromInfo.table
│   ├── hg38.genome.refined.table
│   ├── hg38_repeat
│   └── reformat.sh # for generate hg38.genome.refined.table
├── tumor_purity
│   └── LIHC_purity.csv # generate by Part1
└── WGBS 
    ├── normal
    │   ├── CTR101.bam
    │   └── CTR101.bam.bai
    └── tumor
        ├── HOT151.bam
        └── HOT151.bam.bai
```
## Part1: Marker Selection
The related environment:
```
conda env create -f env/part1_env.yaml
```
In this part we  contain 6 steps:
```
Scripts/
├── Part1.Marker_Selection
│   ├── 1.1_reformat_TCGA.R
│   ├── 1.1_subsample_for_tumorOnly.R
│   ├── 1.2_add_complementary_bins.py
│   ├── 1.2_define_blocks_according_CancerLocator_method.R
│   ├── 1.2_refine_blocks.sh
│   ├── 1.3.1_get_methylation_ratio_blocks_TCGA_paired.R
│   ├── 1.3.2_get_methylation_ratio_blocks_TCGA_tumorOnly.R
│   ├── 1.3.3_extract_reads_from_normal_plasma.sh
│   ├── 1.3.4_get_methylation_ratio_blocks_normal_plasma.R
│   ├── 1.4.1_gen_rand_runs_paired.py
│   ├── 1.4.1_submit_paired_PH.sh
│   ├── 1.4.2_gen_rand_runs_tumorOnly.py
│   ├── 1.4.2_submit_tumorOnly_THorMH.sh
│   ├── 1.5_para_est_mom_mc.R
│   ├── 1.5_submit_est.sh
│   ├── 1.6.1_sel_markers_paired.py
│   ├── 1.6.1_submit_sel_markers_paired.sh
│   ├── 1.6.2_sel_markers_tumorOnly.py
│   ├── 1.6.2_submit_sel_markers_tumorOnly.sh
│   └── 1.6.3_add_complementary_bins.py
Result/
├── 1.1_450K_reformate_data
├── 1.2_define_bins
├── 1.3_cal_methylratio_bins
├── 1.3_extract_reads
├── 1.4_gen_rand_runs
├── 1.5_train_model_for_parameters
├── 1.6_select_marker
```
< Parameters (args) >
---
- root_dir: Directory where DeepLB is located
- tumor: TCGA Abbreviations (LIHC/lihc)
- group : PH TH MH
- subsample : top30 for TH and subsample30 for MH
- annotation_file: tumor tissue and normal plasma annotation
- normal_sample_list: normal plasma list
- marker_selection_threshold: for marker selection ("0.1 0.15 0.2")
- validation: cohorts for cross validation 

---

### 1.1 Prepare 450K data

```shell
#for TH or MH group subsample
#args : root_dir tumor group tumor sample number
Rscript 1.1_subsample_for_tumorOnly.R /home/yinliang/PROJECT/DeepLB LIHC TH 30
Rscript 1.1_subsample_for_tumorOnly.R /home/yinliang/PROJECT/DeepLB LIHC MH 30
#for PH/TH/MH 
#args: root_dir tumor
Rscript 1.1_reformat_TCGA.R /home/yinliang/PROJECT/DeepLB LIHC

```
### 1.2 Identify CpG Clusters
this step include 3 scripts to finish, run the 1.2_refine_blocks.sh to process and remember check the file path
```shell
#args: root_dir tumor
Rscript 1.2_define_blocks_according_CancerLocator_method.R /home/yinliang/PROJECT/DeepLB
bash 1.2_refine_blocks.sh /home/yinliang/PROJECT/DeepLB
```

### 1.3 Quantifying the methylation level of CpG clusters
this step inclue 3 scripts to calculate the CpG Clusters methylation status of 450K and WGBS:
- 1.3.1 for pairwise method, calculate the tumor tissue and ajacent normal tissue
- 1.3.2 for tumorOnly method, only calculate the all tumor tissue
- 1.3.3-1.3.4 for normal plasma sample WGBS data, first extract the part in previous defined CpG Cluster and then calculate the methylation status
```shell
# args：root_dir tumor annotation_file group subsample_num plasma_sample_list
Rscript 1.3.1_get_methylation_ratio_blocks_TCGA_paired.R /home/yinliang/PROJECT/DeepLB lihc all_samples_annotation.txt
Rscript 1.3.2_get_methylation_ratio_blocks_TCGA_tumorOnly.R /home/yinliang/PROJECT/DeepLB lihc TH 30
bash 1.3.3_extract_reads_from_normal_plasma.sh -i /home/yinliang/PROJECT/DeepLB -s sample_list.normal
Rscript 1.3.4_get_methylation_ratio_blocks_normal_plasma.R /home/yinliang/PROJECT/DeepLB sample_list.normal 
```
### 1.4 Split Cohort for cross validation
This step is for generate multi cohort to cross validation
```shell
# args：root_dir tumor annotation_file group subsample plasma_sample_list
bash 1.4.1_submit_paired_PH.sh -i /home/yinliang/PROJECT/DeepLB -g PH -t lihc -s all_samples_annotation.txt
bash 1.4.2_submit_tumorOnly_THorMH.sh -i /home/yinliang/PROJECT/DeepLB -g TH -t lihc -s all_samples_annotation.txt -v 2 -n top30
```
### 1.5 Estimate statistical parameters 
This step is for estimating parameters for marker selection, including three files ,alphas/betas/dist2norm
```shell
# args：root_dir tumor annotation_file group subsample
bash 1.5_submit_est.sh -i /home/yinliang/PROJECT/DeepLB -g PH -t lihc -s paired -v 1 -n ""
bash 1.5_submit_est.sh -i /home/yinliang/PROJECT/DeepLB -g TH -t lihc -s tumorOnly -v 1 -n "top30"
```

### 1.6 Marker Selection
```shell
# args：root_dir tumor annotation_file group subsample_num marker_selection_threshold(-c) and validation(-v)
bash 1.6.1_submit_sel_markers_paired.sh -i /home/yinliang/PROJECT/DeepLB -g PH -t lihc -v 1 -c "0.1 0.15"
bash 1.6.2_submit_sel_markers_tumorOnly.sh -i /home/yinliang/PROJECT/DeepLB -g TH -t lihc -v 1 -c "0.1 0.15" -n "top30"
```

## Part2: Pseudo-fragment Generation by mMTS
![Alt text](png2.png)
```
Scripts/Part2.Pseudo-fragment_Generation_by_mMTS/
├── 2.0_prepare_files_paired.py
├── 2.0_prepare_files_tumorOnly.py
├── 2.1_gen_DMRref_gene_more1.py
├── 2.2_generate_cpg_ref_beta.py
├── 2.3_gen_dmr_reads.py
├── 2.4_simulate_methy_improve.py
├── 2.5_select_healthy_bam_to_train.py
├── 2.6_select_bam_test_all.py
├── 2.7_trans_to_train_read.py
├── 2.8_extract_bam_to_reads.py
├── env_module.py # Must change the direction before process this part!! 
├── mMTS-pipeline.py
Result/2.simulation_result/
├── paired
│   └── lihc-PH
│       ├── 2-0_process
│       ├── 2-1_ref_cpg_beta
│       ├── 2-2_ref_reads
│       ├── 2-3_simulated_reads
│       ├── 2-4_reads_to_train
│       ├── 2-5_real_process_bam
│       ├── 2-6_reads_to_test
│       └── log
└── tumorOnly
```
### This pipeline include 9 steps
- step0 : prepare sample list,marker bed file, markers' CpG methylation status file and CpG site position file
- step1 : obtain the markers reference gene from GRCH38.fa
- step2 : prepare each CpG site's methylation status in markers
- step3 : generate reference DNA reads from markers reference gene
- step4 : simulate methylation status for reference DNA reads(methy.reads)
- step5-6 : obtain the part of markers from plasma WGBS bam file 
- step7 : transfer methy.reads to pseudo-fragments for model training
- step8 : extract DNA reads and methy status from bam file of plasma sample

### Warning!!
Before use mMTS, please check the env_module.py and make sure each file path is right

### Use
```shell
python -u mMTS-pipeline.py -t lihc -s 0.2 -m hypo -g PH -a paired -r 2 -l 1 -L 66 -D 3 -e "" -c freq-diff -p paired
python -u mMTS-pipeline.py -t lihc -s 0.2 -m hypo -g TH -a tumorOnly -r 2 -l 1 -L 66 -D 3 -e "" -c mean-diff -p TumorOnly -p top30
```
### parameters:
- -t : tumor type 
- -s : marker selection threshold
- -g : group(TH MH PH)
- -a : approach(Pairwise or TumorOnly)
- -r : repeat times
- -l : start list
- -c : CpG Clusters methylation status(mean-diff or freq-diff)
- -p : tumor purity (top30 or subsamples)
- -e : except steps
```
Result/3.ResTran_results
├── train_result
│   └── lihc-PH
│       ├── 1_0.4
└── train_result
│       ├── 1_0.4
```

## Part3: Deep Learning (ResTran) model training
```
conda env create -f env/part3_env.yaml
```
```SHELL

bash ResTran.sh -i /home/yinliang/PROJECT/DeepLB -s 1 -t lihc -g PH -m hypo -c 0.4 -v 1
```
Parameters
- i：root dir
- s: sample list
- t :tumor
- g: group
- m: marker type
- c :cutoff/threshold
- v: validation times


# Citation

**If you use this code for your research, please cite paper:**

Kang S, Li Q, Chen Q, Zhou Y, Park S, Lee G, Grimes B, Krysan K, Yu M, Wang W, Alber F, Sun F, Dubinett SM, Li W, Zhou XJ. CancerLocator: non-invasive cancer diagnosis and tissue-of-origin prediction using methylation profiles of cell-free DNA. Genome Biol. 2017 Mar 24;18(1):53. doi: 10.1186/s13059-017-1191-5. PMID: 28335812; PMCID: PMC5364586.

Li W, Li Q, Kang S, Same M, Zhou Y, Sun C, Liu CC, Matsuoka L, Sher L, Wong WH, Alber F, Zhou XJ. CancerDetector: ultrasensitive and non-invasive cancer detection at the resolution of individual reads using cell-free DNA methylation sequencing data. Nucleic Acids Res. 2018 Sep 6;46(15):e89. doi: 10.1093/nar/gky423. PMID: 29897492; PMCID: PMC6125664.

Li J, Wei L, Zhang X, Zhang W, Wang H, Zhong B, Xie Z, Lv H, Wang X. DISMIR: Deep learning-based noninvasive cancer detection by integrating DNA sequence and methylation information of individual cell-free DNA reads. Brief Bioinform. 2021 Nov 5;22(6):bbab250. doi: 10.1093/bib/bbab250. PMID: 34245239; PMCID: PMC8575022.

