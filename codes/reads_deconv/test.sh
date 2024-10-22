#samtools view /home/cmb-panasas2/xjzhou/WGBS/rao_data/cfDNA/simu_data/combination1/bam_data/NIC1254A18_CTR111/theta0.2_cov1.0/NIC1254A18_CTR111.bam | ./reads_binning -b ./data/binning_data/biomarkers.all_bins -r stdin > test/NIC1254A18_CTR111.theta0.2_cov1.0.reads_bins

samples_list="TBR36prePlsm TBR36postPlsm TBR36post3M TBR36post6M TBR36post12M TBR34prePlsm TBR34postPlsm TBR34_Post_2M_plamsa"

#samtools view /home/cmb-05/xjz/shulikan/data/CUHK/bam_files/dedup_sorted/${sample}.sort.bam | ./reads_binning -b ./data/binning_data/biomarkers.all_bins -r stdin > test/${sample}.reads_bins

#./calc_reads_prob -r test/out.debug_reads.bins_updated -t ./data/binning_data/biomarkers.5tissues.only_bins > test/debug.reads.5tissues.prob
#./calc_reads_prob -r test/out.debug_reads.bins_updated -t ./data/binning_data/biomarkers.14tissues.only_bins > test/debug.reads.14tissues.prob

./em -x 100 -p ./test/tmp/debug.reads.5tissues.prob
./em -x 100 -p ./test/tmp/debug.reads.14tissues.prob
