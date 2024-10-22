# FOR processing normal WGBS bam file
get_line(){
  sed -n ''$2'p' $1
}

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}

while getopts "i:s:" opt; do
  case $opt in
    i) input_dir="$OPTARG"
       ;;
	s) sample_list="$OPTARG"
	   ;;
    \?) echo "无效的选项: -$OPTARG" >&2
        exit 1
       ;;
  esac
done

# parameters
pe_suffix='.bam'
#se_suffix='_trimmed_bismark_bt2.deduplicated.bam'

min_cpg_sites_for_reads="1"
bam_suffix=$pe_suffix

#public_dir="/bigData/qingjial/project/DATA_liquid_biopsy/"
samples_dir=${input_dir}"/Predata/WGBS/normal"
project_dir=${input_dir}

all_markers_file=${project_dir}"/Result/1.2_define_bins/blocks_D200_S3.all_bins"
samples_list_file=${project_dir}"/Predata/metadata/"${sample_list} #WGBS data list
out_dir=${project_dir}"/Result/1.3_extract_reads/all_reads_${min_cpg_sites_for_reads}sites_with_pos_and_sam"
checkMakeDirectory $out_dir
job_log_dir=${project_dir}"/log/1.3.3/"
checkMakeDirectory $job_log_dir

cmd_dir=${project_dir}"/codes/reads_deconv"

sample_list=`cat $samples_list_file`

for sample in $sample_list; do
	echo $sample
	JOB=${job_log_dir}"a_${sample}_extract_read.pbs"
	out_reads_file="$out_dir/${sample}.reads_bins.gz"
	out_methy_rates_file="$out_dir/${sample}.methy_bins"
	out_sam_reads_file="$out_dir/${sample}.sam"

	# change -s to -p for paired end mode!!
	# we still use -s for inhouse data, treat them as single-end 
	echo "samtools view ${samples_dir}/${sample}${bam_suffix} | ${cmd_dir}/reads_binning -s -b $all_markers_file -c $min_cpg_sites_for_reads -r stdin -B ${out_methy_rates_file} -R $out_sam_reads_file | gzip -f > ${out_reads_file}" > $JOB

	chmod a+x $JOB
	nohup $JOB > ${JOB}.log 2>&1 &
done
#gzip ${out_sam_reads_file}1
#gzip ${out_sam_reads_file}2

