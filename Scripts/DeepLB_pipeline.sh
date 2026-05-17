#!/bin/bash
# Summary: Run DeepLB part1/2/3 pipeline steps with CLI parameters and dry-run support.
# Parameters
DRY_RUN=false
# Part1 parameters
ROOT_DIR="" # /path/to/DeepLB
TUMOR="" # lihc brca stad
GROUP="" # TH MH PH
SUBSAMPLE="" # TH:top30 top50 MH:sub30 sub50
ANNOTATION_FILE="" # all_samples_annotation.txt
NORMAL_SAMPLE_LIST="" # background_for_train.txt
MARKER_SELECTION_THRESHOLD="" # "0.1 0.15 0.2"
VALIDATION="" # validation list 1~10
WINDOW_SIZE="" # 100
MIN_PROBE="" # 3
MODULE="" # part1~3
# Part2 parameters
MARKER_TYPE="" # for generate pseudo fragment data
BEGIN_LIST="" # begin_list
COVERAGE="" # synthesize coverage
GENERATION_THRESHOLD="" # select from MARKER_SELECTION_THRESHOLD
FRAGMENT_LENGTH="" # synthesize fragment length, 66

# Function to display usage information
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options (full names and abbreviations):"
    echo "  --root_dir, -r <ROOT_DIR>          : Root directory of the project"
    echo "  --tumor, -t <TUMOR>               : Tumor sample identifier (e.g., lihc, brca)"
    echo "  --group, -g <GROUP>              : Group identifier (e.g., PH, TH, MH)"
    echo "  --subsample, -s <SUBSAMPLE>      : Subsample identifier (e.g., top30, sub30)"
    echo "  --annotation_file, -a <FILE>     : Annotation file path"
    echo "  --normal_sample_list, -n <FILE> : Normal sample list file path"
    echo "  --marker_selection_threshold, -m <THRESHOLD> : Marker selection threshold (e.g., \"0.1 0.15 0.2\")"
    echo "  --validation, -v <VALIDATION>     : Validation identifier (1~10)"
    echo "  --window_size, -w <SIZE>         : Window size for CpG clusters"
    echo "  --min_probe, -p <COUNT>         : Minimum number of probes"
    echo "  --module, -u <MODULE>            : Module to run (part1, part2, part3, all)"
    echo "  --marker_type, -k <TYPE>        : Marker type for pseudo-fragment generation"
    echo "  --begin_list, -b <FILE>         : Begin list file path"
    echo "  --coverage, -c <COVERAGE>       : Coverage threshold"
    echo "  --generation_threshold, -q <THRESHOLD> : Generation threshold (selected from marker selection threshold)"
    echo "  --fragment_length, -l <LENGTH>  : Fragment length for pseudo-fragment generation"
    echo "  --meta, -x <FILE>               : Meta file path"
    echo "  --dry_run, -d                   : Enable dry-run mode"
    echo "  --help, -h                      : Display this help message"
    exit 1
}

# Parse command-line options
while getopts "a:b:c:dn:g:hk:l:m:n:p:q:r:s:t:u:v:w:x:" opt; do
    case ${opt} in
        a ) ANNOTATION_FILE=$OPTARG ;;
        b ) BEGIN_LIST=$OPTARG ;;
        c ) COVERAGE=$OPTARG ;;
        d|n ) DRY_RUN=true ;;
        g ) GROUP=$OPTARG ;;
        h ) usage ;;
        k ) MARKER_TYPE=$OPTARG ;;
        l ) FRAGMENT_LENGTH=$OPTARG ;;
        m ) MARKER_SELECTION_THRESHOLD=$OPTARG ;;
        n ) NORMAL_SAMPLE_LIST=$OPTARG ;;
        p ) MIN_PROBE=$OPTARG ;;
        q ) GENERATION_THRESHOLD=$OPTARG ;;
        r ) ROOT_DIR=$OPTARG ;;
        s ) SUBSAMPLE=$OPTARG ;;
        t ) TUMOR=$OPTARG ;;
        u ) MODULE=$OPTARG ;;
        v ) VALIDATION=$OPTARG ;;
        w ) WINDOW_SIZE=$OPTARG ;;
        x ) META=$OPTARG ;;
        \? ) echo "Invalid option: -$OPTARG requires an argument" 1>&2
             usage ;;
        : ) echo "Invalid option: -$OPTARG requires an argument" 1>&2
            usage ;;
    esac
done

# Function to display version information
version(){
    echo "DeepLB version: 1.0"
    echo "author: Yinliang ZhanyuLiang"
    echo "updated: 2025-04-07"
}

# Function to execute a command or print it (dry-run mode)
execute_or_print() {
    if $DRY_RUN; then
        echo "DRY-RUN: $1"
        echo "====================================="
    else
        echo "Executing: $1"
        eval $1
        echo "====================================="
    fi
}

# Validate parameters based on the specified module
validate_parameters() {
    case $MODULE in
        part1)
            if [ -z "$ROOT_DIR" ] || [ -z "$TUMOR" ] || [ -z "$GROUP" ] || [ -z "$SUBSAMPLE" ] || [ -z "$ANNOTATION_FILE" ] || [ -z "$NORMAL_SAMPLE_LIST" ] || [ -z "$MARKER_SELECTION_THRESHOLD" ] || [ -z "$VALIDATION" ] || [ -z "$WINDOW_SIZE" ] || [ -z "$MIN_PROBE" ]; then
                echo "Error: Missing required parameters for part1."
                usage
            fi
            ;;
        part2)
            if [ -z "$ROOT_DIR" ] || [ -z "$TUMOR" ] || [ -z "$GROUP" ] || [ -z "$GENERATION_THRESHOLD" ] || [ -z "$MARKER_TYPE" ] || [ -z "$VALIDATION" ] || [ -z "$SUBSAMPLE" ]; then
                echo "Error: Missing required parameters for part2."
                usage
            fi
            ;;
        part3)
            if [ -z "$ROOT_DIR" ] || [ -z "$TUMOR" ] || [ -z "$GROUP" ] || [ -z "$GENERATION_THRESHOLD" ] || [ -z "$MARKER_TYPE" ] || [ -z "$VALIDATION" ]; then
                echo "Error: Missing required parameters for part3."
                usage
            fi
            ;;
        all)
            if [ -z "$ROOT_DIR" ] || [ -z "$TUMOR" ] || [ -z "$GROUP" ] || [ -z "$SUBSAMPLE" ] || [ -z "$ANNOTATION_FILE" ] || [ -z "$NORMAL_SAMPLE_LIST" ] || [ -z "$MARKER_SELECTION_THRESHOLD" ] || [ -z "$VALIDATION" ] || [ -z "$WINDOW_SIZE" ] || [ -z "$MIN_PROBE" ] || [ -z "$GENERATION_THRESHOLD" ] || [ -z "$MARKER_TYPE" ]; then
                echo "Error: Missing required parameters for all parts."
                usage
            fi
            ;;
        *)
            echo "Invalid module specified: $MODULE"
            usage
            ;;
    esac
}

# Validate parameters
validate_parameters

# Part1: Marker Selection
run_part1(){
    cd Part1.Marker_Selection
    echo "Running Part1: Marker Selection"
    tumor=$(echo "$TUMOR" | tr '[:lower:]' '[:upper:]')
    number=$(echo "$SUBSAMPLE" | awk '{gsub(/[^0-9]/,""); print}')
    echo $number
    # 1.1 Prepare 450K data
    echo "Step 1.1: Prepare 450K data"
    echo ">>Running 1.1_reformat_TCGA.R"
    echo "==>output dir:  Result/1.1_450K_reformate_data/ "
    execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.1_reformat_TCGA.R ${ROOT_DIR} ${tumor}"
    execute_or_print "Rscript 1.1_prepare_sample_annotation.R ${ROOT_DIR} ${META} ${tumor} ${NORMAL_SAMPLE_LIST}"

    # 1.1 Subsample for tumorOnly (if SUBSAMPLE is TH or MH)
    if [ "${GROUP}" == "TH" ] || [ "${GROUP}" == "MH" ]; then
        echo "Step 1.1: Subsample for tumorOnly"
        echo ">>running 1.1_subsample_for_tumorOnly.R"
        execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.1_subsample_for_tumorOnly.R ${ROOT_DIR} ${tumor} ${GROUP} ${number} ${ANNOTATION_FILE}"
    else
        echo "Step 1.1: Skipping subsample for paired"
        #touch ${ROOT_DIR}/Predata/tumor_purity/${TUMOR}_${SUBSAMPLE}.txt
    fi
    echo "output dir: Predata/tumor_purity/"

    # 1.2 Identify CpG Clusters
    echo "Step 1.2: Identify CpG Clusters"
    echo ">>running 1.2_define_blocks_according_CancerLocator_method.R"
    execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.2_define_blocks_according_CancerLocator_method.R ${ROOT_DIR} ${WINDOW_SIZE} ${MIN_PROBE}"
    echo ">>running 1.2_refine_blocks.sh"
    execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.2_refine_blocks.sh ${ROOT_DIR} "
    echo "output dir: Result/1.2_define_bins/"

    # 1.3 Quantifying the methylation level of CpG clusters
    echo "Step 1.3: Quantifying the methylation level of CpG clusters"
    echo "processing 450K data of tumor/normal tissue"
    if [ "${GROUP}" == "PH" ] ; then
        echo ">>PH:running 1.3.1_get_methylation_ratio_blocks_TCGA_paired.R"
        execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.3.1_get_methylation_ratio_blocks_TCGA_paired.R ${ROOT_DIR} ${TUMOR} ${ANNOTATION_FILE}"
    else
        echo ">>TH or MH:running 1.3.2_get_methylation_ratio_blocks_TCGA_tumorOnly.R"
        execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.3.2_get_methylation_ratio_blocks_TCGA_tumorOnly.R ${ROOT_DIR} ${TUMOR} ${GROUP} ${number} ${ANNOTATION_FILE}"
    fi
    echo "processing WGBS data of normal plasma"
    echo ">>running 1.3.3_extract_reads_from_normal_plasma.sh"
    execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.3.3_extract_reads_from_normal_plasma.sh -i ${ROOT_DIR} -s ${NORMAL_SAMPLE_LIST}"
    echo ">>running 1.3.4_get_methylation_ratio_blocks_normal_plasma.R"
    execute_or_print "Rscript ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.3.4_get_methylation_ratio_blocks_normal_plasma.R ${ROOT_DIR} ${NORMAL_SAMPLE_LIST}"
    echo "output dir: Result/1.3_cal_methylratio_bins"

    # 1.4 Split Cohort for cross validation
    echo "Step 1.4: Split Cohort for cross validation"
    if [ "${GROUP}" == "PH" ] ; then
        echo ">>running 1.4.1_submit_paired_PH.sh"
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.4.1_submit_paired_PH.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -s ${ANNOTATION_FILE} -v ${VALIDATION}"
    else
        echo ">>running 1.4.2_submit_tumorOnly_THorMH.sh"
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.4.2_submit_tumorOnly_THorMH.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -s ${ANNOTATION_FILE} -v ${VALIDATION} -n ${SUBSAMPLE}"
    fi
    echo "output dir: Result/1.4_gen_rand_runs/"

    # 1.5 Estimate statistical parameters
    echo "Step 1.5: Estimate statistical parameters"
    echo ">>running 1.5_submit_est.sh"
    if [ "${GROUP}" == "PH" ] ; then
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.5_submit_est.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -s paired -v ${VALIDATION} -n \"\""
    else
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.5_submit_est.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -s tumorOnly -v ${VALIDATION} -n \"${SUBSAMPLE}\""
    fi
    echo "output dir: Result/1.5_train_model_for_parameters/"

    # 1.6 Marker Selection
    echo "Step 1.6: Marker Selection"
    if [ "${GROUP}" == "PH" ] ; then
        echo ">>running 1.6.1_submit_sel_markers_paired.sh"
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.6.1_submit_sel_markers_paired.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -v ${VALIDATION} -c \"${MARKER_SELECTION_THRESHOLD}\""
    else
        echo ">>running 1.6.2_submit_sel_markers_tumorOnly.sh"
        execute_or_print "bash ${ROOT_DIR}/Scripts/Part1.Marker_Selection/1.6.2_submit_sel_markers_tumorOnly.sh -i ${ROOT_DIR} -g ${GROUP} -t ${TUMOR} -v ${VALIDATION} -c \"${MARKER_SELECTION_THRESHOLD}\" -n \"${SUBSAMPLE}\""
    fi
    echo "output dir: Result/1.6_select_marker/"
}

# Part2: Pseudo-fragment Generation by mMTS
run_part2(){
    cd "${ROOT_DIR}/Scripts/Part2.Pseudo-fragment_Generation_by_mMTS"
    echo "Running Part2: Pseudo-fragment Generation by mMTS"
    if [ "${GROUP}" == "PH" ] ; then
        CANCER_MARKER_SELECTION="paired"
        FILTER_TYPE="freq-diff"
    else
        CANCER_MARKER_SELECTION="tumorOnly"
        FILTER_TYPE="mean-diff"
    fi
    execute_or_print "python -u ${ROOT_DIR}/Scripts/Part2.Pseudo-fragment_Generation_by_mMTS/mMTS-pipeline.py -t ${TUMOR} -s ${GENERATION_THRESHOLD} -m ${MARKER_TYPE} -g ${GROUP} -a ${CANCER_MARKER_SELECTION} -r ${VALIDATION} -l ${BEGIN_LIST} -L ${FRAGMENT_LENGTH} -D ${COVERAGE} -e \"\" -c ${FILTER_TYPE} -p ${SUBSAMPLE}"
}

# Part3: Deep Learning (ResTran) model training
run_part3(){
    cd "${ROOT_DIR}/Scripts/Part3.ResTran_model_training"
    echo "Running Part3: Deep Learning (ResTran) model training"
    execute_or_print "conda activate part3_env"  # Ensure the correct environment is activated
    execute_or_print "bash ${ROOT_DIR}/Scripts/Part3.ResTran_model_training/ResTran.sh -i ${ROOT_DIR} -s ${VALIDATION} -t ${TUMOR} -g ${GROUP} -m ${MARKER_TYPE} -c ${GENERATION_THRESHOLD} -v ${VALIDATION}"
}

# Main execution based on the specified module
case $MODULE in
    part1)
        version
        run_part1
        ;;
    part2)
        version
        run_part2
        ;;
    part3)
        version
        run_part3
        ;;
    all)
        version
        run_part1
        run_part2
        run_part3
        ;;
    *)
        echo "Invalid module specified: $MODULE"
        usage
        ;;
esac

echo "DeepLB workflow completed successfully!"