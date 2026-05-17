# Summary: Submit parameter estimation jobs for each random run.
#updata 20240311
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

while getopts "i:g:t:s:v:n:" opt; do
  case $opt in
        i) root_dir="$OPTARG"
       ;;
        g) group="$OPTARG"
                ;;
        t) tumor_type="$OPTARG"
                ;;
        s) method_type="$OPTARG"
                ;;
        v) validate_time="$OPTARG"
                ;;
        n) thre_type="$OPTARG"
                ;;
        \?) echo "Invalid option: -$OPTARG" >&2
        exit 1
       ;;
  esac
done

#tumor_type="stad" #"npc""stad""lihc""paad"
#method_type="tumorOnly"    #"paired" 
#thre_type="top30"  
#group="TH"   #"TH" "subsample" "PH"
#project_dir="DeepLB/"

input_dir=$root_dir"/Result/1.4_gen_rand_runs/"${method_type}"/"${tumor_type}"/"${group}"/"${thre_type}
job_log_dir=${root_dir}"/log/1.5/"
checkMakeDirectory $job_log_dir

for repeat_i in $(seq 1 $validate_time);do
        echo $repeat_i
        JOB=$job_log_dir"a_para_est_${method_type}_${repeat_i}.pbs"

	out_dir=$root_dir"/Result/1.5_train_model_for_parameters/"${method_type}"/"${tumor_type}"/"${group}"/"${thre_type}"/MoM_model_"$repeat_i
	mkdir -p $out_dir

	train_file=$input_dir"/train_"$repeat_i
	echo "Rscript $root_dir/Scripts/Part1.Marker_Selection/1.5_para_est_mom_mc.R $train_file $out_dir 6" > $JOB

        chmod a+x $JOB
        nohup $JOB > ${JOB}.log 2>&1 &
        pids+=($!)
done

wait "${pids[@]}"