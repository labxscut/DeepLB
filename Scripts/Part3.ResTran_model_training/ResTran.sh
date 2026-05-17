#!/bin/bash
# Summary: Run ResTran training, prediction, and risk scoring per iteration.

while getopts "i:s:l:t:g:m:c:v:u" flag
do
    case "${flag}" in
        i) root_dir=${OPTARG};;
        s) sample=${OPTARG};;
        t) tumor=${OPTARG};;
        g) group=${OPTARG};;
        m) method=${OPTARG};;
        c) cutoff=${OPTARG};;
        v) iterations=${OPTARG};;
        u) usesample=${OPTARG};;
        *) echo "Usage: $0 -i <root_dir> -s <sample> -t <tumor> -g <group> -m <method> -c <cutoff> -v <iterations>"; exit 1;;
    esac
done

if [ -z "$root_dir" ] || [ -z "$sample" ] || [ -z "$tumor" ] || [ -z "$group" ] || [ -z "$method" ] || [ -z "$cutoff" ] || [ -z "$iterations" ]; then
    echo "Usage: $0 -i <root_dir> -s <sample> -t <tumor> -g <group> -m <method> -c <cutoff> -v <iterations>"
    exit 1
fi



for (( i=$sample; i<=iterations; i++ ))
do
    log_dir="${root_dir}/log/3.Restran/${tumor}_${group}/${i}_${cutoff}_${iterations}/"
    mkdir -p $log_dir
    echo "Starting iteration $i"

    echo "Running training.py - iteration $i"
    python training.py -i "$root_dir" -s "$sample" -t "$tumor" -g "$group" -m "$method" -c "$cutoff" -u "$usesample"> "${log_dir}/training_step_${i}.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error in training.py during iteration $i"
        exit 1
    fi

    echo "Running predict_reads_source.py - iteration $i"
    python predict_reads_source.py -i "$root_dir" -s "$sample" -t "$tumor" -g "$group" -m "$method" -c "$cutoff" > "${log_dir}/predict_step_${i}.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error in predict_reads_source.py during iteration $i"
        exit 1
    fi

    echo "Running cal_risk.py - iteration $i"
    python cal_risk.py -i "$root_dir" -s "$sample" -t "$tumor" -g "$group" -m "$method" -c "$cutoff" > "${log_dir}/cal_risk_step_${i}.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error in cal_risk.py during iteration $i"
        exit 1
    fi

    echo "Iteration $i completed"
done

echo "All iterations completed successfully."
