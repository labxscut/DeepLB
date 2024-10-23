import numpy as np
import os
import argparse
import env_modelue



def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--root', help='root direction') #DeepLB path
    parser.add_argument('-s', '--samples', help='sample list') #1 - 10
    parser.add_argument('-t', '--tumor', help='tumor') 
    parser.add_argument('-c', '--cutoff',  help='cutoff') 
    parser.add_argument('-g', '--group',  help='fgroup name') #PH PLH TH TLH
    parser.add_argument('-m', '--marker',  help='marker type') 
    args = parser.parse_args()
    group = args.group
    tumor = args.tumor
    threshold = args.cutoff
    marker = args.marker
    samples = args.samples
    if group == "PH":
        method = "paired"
    else:
        method = "tumorOnly"
    dir_name =  samples + "_" + threshold 
    result_name = "/CS_" + samples + "_" + threshold + ".txt"

    store_dir = args.root + "/Result/3.ResTran_results/test_result/" + tumor + "_" + group + "/" + dir_name + "/d-scores/" # directory where the predicted d-scores are stored
    output_dir = args.root + "/Result/3.ResTran_results/test_result/" + tumor + "_" + group + "/" + dir_name + "/" + result_name
    output = open(output_dir , 'w') # file to store predicted risk
    files = file_name(store_dir)
    # cal risk by maximize posterior probability
    gaps = np.linspace(0, 1, 1001)
    score = np.vstack((gaps, 1 - gaps))
    score = score.T
    for file in files:
        print(file)
        likelihood_1 = np.loadtxt(store_dir + file)
        likelihood_2 = 1 - likelihood_1
        likelihood = np.vstack((likelihood_1, likelihood_2))
        val = np.log10(np.dot(score, likelihood))
        sum = np.sum(val, axis=1)
        result = gaps[np.argmax(sum)]
        output.write(file + '\t' + str(result) + '\n')
    output.close()

