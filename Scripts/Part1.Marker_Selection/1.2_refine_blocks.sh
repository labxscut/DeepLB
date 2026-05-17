# Summary: Add complementary bins to block definitions using genome sizes.
code_dir=$PWD
current_path=$1 #DeepLB path /home/yinliang/PROJECT/DeepLB
public_dir="${current_path}/Predata/reference/"
biomarker_dir="${current_path}/Result/1.2_define_bins"
mkdir -p $biomarker_dir
echo $biomarker_dir
#parameters
biomarker_file=$biomarker_dir/"blocks_D200_S3.txt"
output_file=$biomarker_dir/"blocks_D200_S3.all_bins"

python $code_dir/1.2_add_complementary_bins.py --chromosome_sizes_file=$public_dir/hg38.genome.refined.table \
                        --biomarker_features_file=$biomarker_file --output_file=$output_file



