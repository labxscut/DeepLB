import numpy as np
import random
import os
import os.path
import re
import argparse

from tensorflow.keras.models import Sequential
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras import optimizers
import tensorflow as tf

from models import DISMIR_deep,CNN_only,resnet1d,cnn_transformer,res_transformer,res_lstm,cnn_lstm_transformer,cnn_lstm_cnn_transformer,lstm,cnn_lstm,transformer_encoder

def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files


# transform sequence into number for storage (A/T/C/G to 0/1/2/3, methylated C to 4)
def lstm_seq(seq, methy):
    i = 0
    lstmseq = np.zeros((len(seq), 66), dtype='int')
    while i < len(seq):
        tmp = seq[i]
        j = 0
        while j < len(tmp):
            if tmp[j] == 'A':
                lstmseq[i, j] = 0
            elif tmp[j] == 'T':
                lstmseq[i, j] = 1
            elif tmp[j] == 'C':
                lstmseq[i, j] = 2
            else:
                lstmseq[i, j] = 3
            if int(methy[i][j]) == 1:
                lstmseq[i, j] = 4
            j = j + 1
        i = i + 1
    return lstmseq


# transform sequence into one-hot code (0/1/2/3 to one-hot) and add methylation state channel
def conv_onehot(seq):
    module = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 1]])
    onehot = np.zeros((len(seq), 66, 5), dtype='int')
    for i in range(len(seq)):
        tmp = seq[i]
        tmp_onehot = np.zeros((66, 5), dtype='int')
        for j in range(len(tmp)):
            if tmp[j] == 0:
                tmp_onehot[j] = module[0]
            elif tmp[j] == 1:
                tmp_onehot[j] = module[1]
            elif tmp[j] == 2:
                tmp_onehot[j] = module[2]
            elif tmp[j] == 3:
                tmp_onehot[j] = module[3]
            else:
                tmp_onehot[j] = module[4]
        onehot[i] = tmp_onehot
    return onehot


# deep learning model

#design  parameters
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
train_dir = args.root + "/Result/3.ResTran_results/train_result/" + tumor + "_" + group + "/" + dir_name + "/" # directory where the model is saved
file_dir =  args.root + "/Result/2.simulation_result/" + method + "/" + tumor + "-" + group + "/2-6_reads_to_test/" + marker + "/" + dir_name +"/"# directory where the sample to test is saved, format same as training data
print(file_dir)
output_dir =  args.root + "/Result/3.ResTran_results/test_result/" + tumor + "_" + group + "/" + dir_name + "/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
store_dir = output_dir + "d-scores/"# directory to save predicted d-scores of reads
if not os.path.exists(store_dir):
    os.makedirs(store_dir)

if __name__ == '__main__':
    gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    model = res_transformer()


    files = file_name(file_dir)
    for file in files:
        print(file)
        input = open(file_dir + file,'r')
        seq = []
        methy = []
        for item in input:
            item = item.split()
            cpg = 0
            for i in range(len(item[2]) - 1):
                if (item[2][i] == 'C') & (item[2][i + 1] == 'G'):
                    cpg = cpg + 1
            if cpg > 2:
                seq.append(item[2])
                methy.append(item[3])
        input.close()
        seq_lstm = lstm_seq(seq, methy)
        seq_3one_hot = conv_onehot(seq_lstm)
        if len(seq_3one_hot)==0:
            pass
        else:
            result = model.predict(seq_3one_hot, verbose=0)
            np.savetxt(store_dir + 'result_' + file + '.txt', result)








