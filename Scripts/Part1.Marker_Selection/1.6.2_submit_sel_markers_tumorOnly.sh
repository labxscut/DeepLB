
checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}
while getopts "i:g:t:v:c:n:" opt; do
  case $opt in
    i) project_dir="$OPTARG"
       ;;
	g) group="$OPTARG"
	   ;;
	t) tumor_type="$OPTARG"
	   ;;
	v) validate_time="$OPTARG"
	   ;;
	c) diff_min_list="$OPTARG"
	   ;;
	n) thre_type="$OPTARG"
	   ;;
    \?) echo "无效的选项: -$OPTARG" >&2
        exit 1
       ;;
  esac
done


#tumor_type="paad"  #stad lihc paad
#diff_min_list="0.1 0.15 0.2 0.25" #0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65"    #"0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.65 0.7"   #"0.1 0.15 0.2 0.25 0.35 0.45 0.55 0.6"        #"0.1 0.15 0.2 0.25 0.35 0.45 0.55" #"0.45 0.5 0.55 0.6" #"0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5" #"0.15 0.2 0.25"
#group="TH"    #"TH"     #"MH"
#thre_type="subsample"    #"top30"
#project_dir="DeepLB/"

all_marker_def_file=${project_dir}"/Result/1.2_define_bins/blocks_D200_S3.txt"
hg38_genome_file=${project_dir}"/Predata/reference/hg38.genome.refined.table"

in_dir=${project_dir}"/Result/1.4_gen_rand_runs/tumorOnly/"${tumor_type}"/"${group}"/"${thre_type}
out_dir=${project_dir}"/Result/1.6_select_marker/tumorOnly/mean-diff/"${tumor_type}"/"${group}"/"${thre_type}
mkdir -p $out_dir
log_dir=${project_dir}"/log/1.6/"${tumor_type}"/"${group}"/"${thre_type}
mkdir -p $log_dir

for diff_min in $diff_min_list; do
	method="mean-diff,$diff_min" #mean-diff
	echo $method
	for i in $(seq 1 $validate_time); do # train data index
		in_params_dir=${project_dir}"/Result/1.5_train_model_for_parameters/tumorOnly/"${tumor_type}"/"${group}"/"${thre_type}"/MoM_model_$i"
		alphas_file="$in_params_dir/alphas"
		betas_file="$in_params_dir/betas"

		out_file="$out_dir/${i}_${diff_min}"
		echo -e "\t${i}_${diff_min}"

		python ${project_dir}"/Scripts/Part1.Marker_Selection/1.6.2_sel_markers_tumorOnly.py" $method $tumor_type $in_dir/train_$i $in_dir/train_${i}.matched $all_marker_def_file $alphas_file $betas_file $out_file > $log_dir"/"${i}_${diff_min}".log"
		python ${project_dir}"/Scripts/Part1.Marker_Selection/1.6.3_add_complementary_bins.py"  --chromosome_sizes_file=$hg38_genome_file --biomarker_features_file=$out_file --output_file=${out_file}.full
	done
done

