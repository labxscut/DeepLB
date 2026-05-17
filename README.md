# DeepLB
## Overview
In this study, we introduce our novel deep learning model, DeepLB, developed for early cancer detection through liquid biopsy. This model utilizes pseudo-fragment data generated from 450K methylation array data of tumor tissues, alongside cfDNA whole-genome bisulfite sequencing (WGBS) data from healthy donors, to enhance diagnostic capabilities. DeepLB consists of four key modules: feature selection, methylation Microarray-To-Sequence data converter (mMTS), deep learning model (ResTran), and cancer score estimation. We analyzed cancer-specific features from tumor tissue and healthy plasma WGBS data, generating tumor pseudo-fragment data through the mMTS converter. The ResTran model, which integrates residual networks with Transformer architecture, accurately identifies tumor-derived DNA fragments from cfDNA methylation sequencing data. The workflow in DeepLB is illustrated in the figure below.

<div align="center">
  <img src="png1.png" alt="DeepLB framework overview" width="80%"/>
  <p><em>Figure 1: Overview of the DeepLB framework.</em></p>
</div>

<div align="center">
  <img src="png2.png" alt="DeepLB computational pipeline" width="80%"/>
  <p><em>Figure 2: The computational pipeline of DeepLB.</em></p>
</div>

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Code organization by language](#code-organization-by-language)
- [Prepare data](#prepare-data)
  - [Warning!!](#warning)
- [Part1: Marker Selection](#part1-marker-selection)
- [Part2: Pseudo-fragment Generation by mMTS](#part2-pseudo-fragment-generation-by-mmts)
- [Part3: Deep Learning (ResTran) model training](#part3-deep-learning-restran-model-training)
- [Use](#use)
- [Citation](#citation)


## Installation

### Prerequisites
- R >= 4.0
- Python >= 3.9

### Install DeepLB: R package + Python wrapper
Installs the R analysis package and an optional Python CLI wrapper to run the pipeline.

```bash
# Install the R package
Rscript install_packages.R

# (Optional) Install the Python wrapper (new, compatibility-first)
# The wrapper simply calls the existing shell pipeline without changing its behavior.
# Public API docs: https://labxscut.github.io/deeplb/

pip install -e .

# From the repo root, you can test it with:
python -m deeplb -- -h

# For an editable install that provides a `deeplb` command:
deeplb -- -h
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
    └── Part3.ResTran_model_training
├── deeplb # compatibility-first Python API/CLI wrapper
├── docs
├── vignettes
```

## Prepare data
!!!!! The data should be prepare before Process DeepLB

[GDC Xena Hub (TCGA-LIHC Methylation Data Entry)](https://xenabrowser.net/datapages/?host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
```
Predata/
├── 450K #Download from GEO and TCGA
│   ├── illuminaMethyl450_hg38_GDC 
│   ├── TCGA_450_probe 
│   ├── TCGA-LIHC.methylation450.tsv #tumor 450K data
│   └── TCGA_Study_Abbreviations.txt
├── metadata
│   ├── all_samples_annotation.txt #tumor tissue and normal plasma annotation
│   ├── background_for_train.txt # normal plasma list
│   └── all_WGBS_sample_metadata.xlsx # plasma sample metadata
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
```

### Warning!!
Before use mMTS, please check the env_module.py and make sure each file path is right

## Code organization by language
DeepLB is a mixed-language pipeline. It is organized by pipeline stage, and each
stage combines R/Python/shell differently.

- **Shell orchestration (`.sh`)**
  - Top-level pipeline entrypoint: `Scripts/DeepLB_pipeline.sh`
  - Stage wrappers/job submit scripts:
    - `Scripts/Part1.Marker_Selection/*.sh`
    - `Scripts/Part3.ResTran_model_training/ResTran.sh`
  - Role: glue code that sequences steps and passes arguments to R/Python tools.

- **R analysis/statistics (`.R`)**
  - Mainly in `Scripts/Part1.Marker_Selection/`
  - Typical tasks: 450K preprocessing, bin-level methylation ratio computation,
    parameter estimation.
  - Examples: `1.1_reformat_TCGA.R`, `1.3.1_get_methylation_ratio_blocks_TCGA_paired.R`,
    `1.5_para_est_mom_mc.R`.

- **Python processing/modeling (`.py`)**
  - Part1 helpers in `Scripts/Part1.Marker_Selection/` (selection/splitting helpers)
  - Part2 mMTS pipeline in `Scripts/Part2.Pseudo-fragment_Generation_by_mMTS/`
  - Part3 ResTran training/inference in `Scripts/Part3.ResTran_model_training/`
  - Role: data transformation, pseudo-fragment generation, deep model train/predict.

- **C/C++ performance utilities**
  - `codes/reads_deconv/` (compiled utilities used by older/low-level steps)
  - Role: read/binning/probability calculations for methylation read processing.

- **Python compatibility API layer (new, optional)**
  - `deeplb/` (`pipeline.py`, `cli.py`)
  - Role: thin Python API/CLI wrappers that call the existing shell pipeline
    without changing core behavior.

In short: **Part1 is R+Python+shell**, **Part2 is mostly Python (+ shell entry)**,
**Part3 is Python model code + shell runner**, and **`codes/` keeps lower-level C/C++ tools**.

## Part1: Marker Selection
- 1.1 Prepare 450K data
- 1.2 Identify CpG Clusters
- 1.3 Quantifying the methylation level of CpG clusters
- 1.4 Split Cohort for cross validation
- 1.5 Estimate statistical parameters
- 1.6 Marker Selection
  
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

## Part2: Pseudo-fragment Generation by mMTS
- step0 : prepare sample list,marker bed file, markers' CpG methylation status file and CpG site position file
- step1 : obtain the markers reference gene from GRCH38.fa
- step2 : prepare each CpG site's methylation status in markers
- step3 : generate reference DNA reads from markers reference gene
- step4 : simulate methylation status for reference DNA reads(methy.reads)
- step5-6 : obtain the part of markers from plasma WGBS bam file 
- step7 : transfer methy.reads to pseudo-fragments for model training
- step8 : extract DNA reads and methy status from bam file of plasma sample
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

## Part3: Deep Learning (ResTran) model training
```
Scripts/Part3.ResTran_model_training/
├── training.py
├── models.py
├── predict_reads_source.py
├── cal_risk.py
├── ResTran.sh
Result/3.ResTran_results
├── train_result
│   └── lihc-PH
│       ├── 1_0.4
└── test_result
  └── lihc-PH
    ├── 1_0.4
``` 

# Use
```
Usage: DeepLB_pipeline.sh [OPTIONS]
Options (full names and abbreviations):
  --root_dir, -r <ROOT_DIR>          : Root directory of the project
  --tumor, -t <TUMOR>               : Tumor sample identifier (e.g., lihc, brca)
  --group, -g <GROUP>              : Group identifier (e.g., PH, TH, MH)
  --subsample, -s <SUBSAMPLE>      : Subsample identifier (e.g., top30, sub30)
  --annotation_file, -a <FILE>     : Annotation file path
  --normal_sample_list, -n <FILE> : Normal sample list file path
  --marker_selection_threshold, -m <THRESHOLD> : Marker selection threshold (e.g., "0.1 0.15 0.2")
  --validation, -v <VALIDATION>     : Validation identifier (1~10)
  --window_size, -w <SIZE>         : Window size for CpG clusters
  --min_probe, -p <COUNT>         : Minimum number of probes
  --module, -u <MODULE>            : Module to run (part1, part2, part3, all)
  --marker_type, -k <TYPE>        : Marker type for pseudo-fragment generation
  --begin_list, -b <FILE>         : Begin list file path
  --coverage, -c <COVERAGE>       : Coverage threshold
  --generation_threshold, -q <THRESHOLD> : Generation threshold (selected from marker selection threshold)
  --fragment_length, -l <LENGTH>  : Fragment length for pseudo-fragment generation
  --meta, -x <FILE>               : Meta file path
  --dry_run, -d                   : Enable dry-run mode
  --help, -h                      : Display this help message
```

Compatibility wrapper usage:
```bash
python -m deeplb -- -r /path/to/DeepLB -t lihc -g TH -u part3 -q 0.1 -k hyper -v 1
python -m deeplb --dry-run -- -r /path/to/DeepLB -t lihc -g TH -u part3 -q 0.1 -k hyper -v 1
```
```
#Examples
# Only part1
bash DeepLB_pipeline.sh -r /path/to/DeepLB -t lihc -g TH -s top30 -a all_samples_annotation.txt -n background_for_train.txt -m "0.1 0.15 0.2" -v 1 -w 100 -p 3 -u part1
# Only part2
bash DeepLB_pipeline.sh -r /path/to/DeepLB -t lihc -g TH -q 0.1 -k "hyper" -v 1 -s top30 -b begin_list.txt -c 3 -l 66 -u part2
# Only part3
bash DeepLB_pipeline.sh -r /path/to/DeepLB -t lihc -g TH -q 0.1 -k "hyper" -v 1 -u part3
# ALL
bash DeepLB_pipeline.sh -r /path/to/DeepLB -t lihc -g TH -s top30 -a all_samples_annotation.txt -n background_for_train.txt -m "0.1 0.15 0.2" -v 1 -w 100 -p 3 -q 0.1 -k marker_type -b begin_list.txt -c 30 -l 66 -u all

```

For detailed function documentation of the scripts inside `Scripts/` (including R, Python, and shell scripts), please refer to:

- **Public API Docs (Python wrapper & pipeline scripts):** [https://labxscut.github.io/deeplb/](https://labxscut.github.io/deeplb/)

This documentation covers the usage of the `deeplb` Python package as well as the individual scripts used in Part1, Part2, and Part3.

# Citation

## Primary citation
If you use DeepLB in your research, please cite the main paper:

> Yin Liang, Zhanyu Liang, Xiaoxin Yao, Yajie Guo, Xufeng Kong, Guangquan Zhang,
> Jia-Bin Wang, Kaida Ning, Yulin Liu, Tong Wang, Qingjiao Li, Li C. Xia.
> **DeepLB: Noninvasive Early Cancer Detection via Deep Learning on
> Pseudo-Fragments with Methylation and Sequence Features.** *In press.*

## Required citations for integrated methods
DeepLB builds upon the following methods. If your analysis uses their specific
modules, please also cite the corresponding works:

| Method | Module / Origin | Citation |
|--------|----------------|----------|
| CancerLocator | CpG cluster definition and marker selection framework | Kang S, Li Q, Chen Q, et al. *Genome Biol.* 2017;18:53. DOI: 10.1186/s13059-017-1191-5 |
| CancerDetector | Probabilistic read-level scoring model | Li W, Li Q, Kang S, et al. *Nucleic Acids Res.* 2018;46:e89. DOI: 10.1093/nar/gky423 |
| DISMIR | Deep learning architecture reference (CNN_LSTM_CNN baseline) | Li J, Wei L, Zhang X, et al. *Brief Bioinform.* 2021;22:bbab250. DOI: 10.1093/bib/bbab250 |

# Acknowledgement

Our model implementation references the open‑source code of DISMIR
(https://github.com/XWangLabTHU/DISMIR). We thank the authors for making their
work publicly available.

