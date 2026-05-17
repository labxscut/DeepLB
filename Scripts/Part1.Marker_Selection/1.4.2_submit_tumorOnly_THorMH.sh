# Summary: Submit tumor-only random-run generation across repeats and purity settings.
# make dir "paired_v2" before running the follwing command
# An example of the partition file is below:
#run_id	test_list	train_list
#1	CTR101,CTR108,CTR110,CTR113,CTR114,CTR126,CTR129,CTR84	CTR103,CTR104,CTR106,CTR107,CTR111,CTR117,CTR118,CTR127,CTR128,CTR131,CTR132,CTR134,CTR147,CTR148,CTR149,CTR150,CTR151,CTR152,CTR153,CTR154,CTR85,CTR86,CTR97,CTR98
#2	CTR106,CTR111,CTR127,CTR147,CTR150,CTR151,CTR153,CTR98	CTR101,CTR103,CTR104,CTR107,CTR108,CTR110,CTR113,CTR114,CTR117,CTR118,CTR126,CTR128,CTR129,CTR131,CTR132,CTR134,CTR148,CTR149,CTR152,CTR154,CTR84,CTR85,CTR86,CTR97
#3	CTR103,CTR104,CTR107,CTR132,CTR134,CTR148,CTR152,CTR86	CTR101,CTR106,CTR108,CTR110,CTR111,CTR113,CTR114,CTR117,CTR118,CTR126,CTR127,CTR128,CTR129,CTR131,CTR147,CTR149,CTR150,CTR151,CTR153,CTR154,CTR84,CTR85,CTR97,CTR98
#"top30"
while getopts "i:g:t:s:v:n:" opt; do
  case $opt in
    i) root_dir="$OPTARG"
       ;;
	g) group="$OPTARG"
	   ;;
	t) tumor="$OPTARG"
	   ;;
	s) samplelist="$OPTARG"
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

output_dir=${root_dir}"/Result/1.4_gen_rand_runs/tumorOnly" 
echo $output_dir
#partition_file="/u/project/xjzhou/qingjiao/projects/CancerDetector/workspace/normal_plasma_rand_partitions/partition_10runs.for_paper.dennis_lo"
#partition_file="/u/project/xjzhou/qingjiao/projects/TissueDeconvolution_paper/May18_2018/results_from_CancerLocator/genome_biology_paper/partition_10runs.for_paper.dennis_lo"
#line=$(($SGE_TASK_ID+1))
#trainingNormalPlasmaList=`sed "${line}q;d" $partition_file | cut -f 3` # get the line and its 3rd item, which is a list of training normal samples. In this list, all items are deliminated by a comma

for runID in $(seq 1 $validate_time)
do
	for cancerType in $tumor #lihc stad paad
	 do
		mkdir -p ${output_dir}/${cancerType}/${group}/${thre_type}
		echo  ${output_dir}/${cancerType}/${group}/${thre_type}
		python 1.4.2_gen_rand_runs_tumorOnly.py ${runID} ${cancerType} ${output_dir}/${cancerType}/${group}/${thre_type} ${root_dir} ${samplelist} ${thre_type}
	done
done
