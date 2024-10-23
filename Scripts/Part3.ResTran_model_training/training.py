"""
update:20241023
Yin Liang; ZhanYu Liang
"""

import numpy as np
import random
import gc
import os
import os.path
import re
import argparse
os.environ['KERAS_BACKEND']='tensorflow'

from keras.models import Sequential
from keras import layers
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers
import tensorflow as tf
from models import DISMIR_deep,CNN_only,resnet1d,cnn_transformer,res_transformer,res_lstm,cnn_lstm_transformer,cnn_lstm_cnn_transformer,lstm,cnn_lstm,transformer_encoder

#define fomulate
def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        return files
    
def data_prepare(file_dir1, file_dir2):
    #file_dir1:helthy
    #file_dir2:tumor
    normal_seq = []
    normal_methy = []
    normal_chrom = []
    normal_region = []
    normal_sample = []
    files1 = file_name(file_dir1)
    files2 = file_name(file_dir2)
    for file in files1:
        normal_sample.append(file)
        input_norm = open(file_dir1 + file, 'r')
        for item in input_norm:
            item = item.split()
            cpg = 0
            if len(item[2]) == 66:  # check the length
                for i in range(len(item[2]) - 1):
                    if (item[2][i] == 'C') & (item[2][i + 1] == 'G'):
                        cpg = cpg + 1
                if cpg > 2: # filter out reads with less than 3 CpG sites
                    normal_chrom.append(item[0])
                    normal_region.append(item[1])
                    normal_seq.append(item[2])
                    normal_methy.append(item[3])
        input_norm.close()

    # prepare reads from cancer tissue
    tumor_seq = []
    tumor_methy = []
    tumor_chrom = []
    tumor_region = []
    tumor_sample = []
    files2 = file_name(file_dir2)
    flag = 0
    for file in files2:
        tumor_sample.append(file)
        input_tumor = open(file_dir2 + file, 'r')
        for item in input_tumor:
            item = item.split()
            cpg = 0
            if len(item[2]) == 66:
                for i in range(len(item[2]) - 1):
                    if (item[2][i] == 'C') & (item[2][i + 1] == 'G'):
                        cpg = cpg + 1
                if cpg > 2: # filter out reads with less than 3 CpG sites
                    tumor_chrom.append(item[0])
                    tumor_region.append(item[1])
                    tumor_seq.append(item[2])
                    tumor_methy.append(item[3])
        input_tumor.close()
    flag = flag + 1
    return normal_seq, normal_methy, normal_chrom, normal_region, tumor_seq, tumor_methy, tumor_chrom, tumor_region, normal_sample,tumor_sample


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

def make_chrom_region(chrom_0, region_0):
    i = 0
    chrom = np.zeros(len(chrom_0), dtype='int')
    region = np.zeros(len(region_0), dtype='int')
    while i < len(chrom_0):
        chrom[i] = int(re.findall('(\d+)', chrom_0[i])[0])
        region[i] = region_0[i]
        i = i + 1
    return chrom, region


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




if __name__ == '__main__':
    dir_name =  samples + "_" + threshold 
    train_dir = args.root + "/Result/3.ResTran_results/train_result/" + tumor + "_" + group + "/" + dir_name + "/"# directory to store the results and data
    file_dir = args.root + "/Result/2.simulation_result/" + method + "/" + tumor + "-" + group + "/2-4_reads_to_train/" # directory of input data
    if not os.path.exists(train_dir):
        os.makedirs(train_dir)

    # preprocess of input data
    # the first column : the chromosome where the read is from (e.g. chr1)
    # the second column : the position of the read on chromosome (int format)
    # the third column of each row in the file should be the sequence of a read with a length of 66bp (A/T/C/G)
    # the fourth column of each row in the file should be the methylation state corresponding to the sequence
    # also with a length of 66bp (1:methylated/0:unmethylated) both columns in the form of 'str'
    file_dir_tumor = file_dir + "tumor/" + marker + "/" + dir_name + "/"
    file_dir_normal = file_dir + "healthy/" + marker + "/" + dir_name + "/"
    normal_seq, normal_methy, normal_chrom_0, normal_region_0, tumor_seq, tumor_methy, tumor_chrom_0, tumor_region_0,normal_sample,tumor_sample = data_prepare(file_dir_normal,file_dir_tumor)
    print("###########################loading sample########################")
    print(f"normal sample:{normal_sample}\n")
    print(f"tumor sample:{tumor_sample}\n")
    print(f"normal sample:{normal_seq[1]}\n")
    print(f"tumor sample:{tumor_seq[1]}\n")
    r = min(len(tumor_chrom_0),len(normal_chrom_0))

    normal_seq_lstm = lstm_seq(normal_seq, normal_methy)
    tumor_seq_lstm = lstm_seq(tumor_seq, tumor_methy)
    normal_chrom, normal_region = make_chrom_region(normal_chrom_0, normal_region_0)
    tumor_chrom, tumor_region = make_chrom_region(tumor_chrom_0, tumor_region_0)

    # randomize before class balance
    perm_0 = random.sample(range(len(normal_seq)), len(normal_seq))
    normal_seq_lstm = normal_seq_lstm[perm_0]
    normal_chrom = normal_chrom[perm_0]
    normal_region = normal_region[perm_0]

    normal_seq_3one_hot = conv_onehot(normal_seq_lstm[0:r])
    tumor_seq_3one_hot = conv_onehot(tumor_seq_lstm[0:r])
    perm = random.sample(range(len(normal_seq[0:r]) + len(tumor_seq[0:r])), len(normal_seq[0:r]) + len(tumor_seq[0:r]))
    data_lstm_all = np.vstack((normal_seq_lstm[0:r], tumor_seq_lstm[0:r]))
    data = np.vstack((normal_seq_3one_hot[0:r], tumor_seq_3one_hot[0:r]))
    label_all = np.array([0] * len(normal_seq[0:r]) + [1] * len(tumor_seq[0:r]))  # generating labels
    chrom_all = np.hstack((normal_chrom[0:r], tumor_chrom[0:r]))
    region_all = np.hstack((normal_region[0:r], tumor_region[0:r]))

    # randomly mixture two types of data
    chrom_all = chrom_all[perm]
    region_all = region_all[perm]
    data = data[perm]
    data_lstm_all = data_lstm_all[perm]
    label_all = label_all[perm]

    del normal_seq_lstm, tumor_seq_lstm, normal_seq_3one_hot, tumor_seq_3one_hot
    del normal_seq, normal_methy, normal_chrom_0, normal_region_0, tumor_seq, tumor_methy, tumor_chrom_0, tumor_region_0, perm_0
    gc.collect()

    train_num = int(len(data) * 0.7) # 80% as training set, 75% among them for training, 25% among them for validation
    test_num = int(len(data) * 0.3)
    train_data = data[0:train_num]
    train_label = label_all[0:train_num]
    test_data = data[train_num:(train_num + test_num)]
    test_label = label_all[train_num:(train_num + test_num)]
    test_chrom = chrom_all[train_num:(train_num + test_num)]
    test_region = region_all[train_num:(train_num + test_num)]

    data_lstm = data_lstm_all[0:(train_num + test_num)]
    label = label_all[0:(train_num + test_num)]
    chrom = chrom_all[0:(train_num + test_num)]
    region = region_all[0:(train_num + test_num)]
    #data_origin = data_lstm_all[(train_num + test_num)::]  # 20% reserved for other analysis such as motif analysis
    #label_origin = label_all[(train_num + test_num)::]
    #chrom_origin = chrom_all[(train_num + test_num)::]
    #region_origin = region_all[(train_num + test_num)::]

    #np.savetxt(train_dir + 'data_origin_lstm.txt', data_origin)
    #np.savetxt(train_dir + 'label_origin.txt', label_origin)
    #np.savetxt(train_dir + 'chrom_origin.txt', chrom_origin)
    #np.savetxt(train_dir + 'region_origin.txt', region_origin)

    # building and training the model
    gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    model = res_transformer()

    #tf2.2代码
    #sgd = optimizers.SGD(lr=0.05, decay=1e-6, momentum=0.9, nesterov=True)

    #tf2.16
    #sgd = optimizers.SGD(learning_rate=0.05, weight_decay=1e-6, momentum=0.9, nesterov=True)
    sgd = optimizers.AdamW(learning_rate=1e-4, weight_decay=1e-6)
    model.compile(optimizer=sgd, loss='binary_crossentropy', metrics=['accuracy'])
    early_stopping = EarlyStopping(monitor='val_loss', patience=10)
    history = model.fit(train_data, train_label, epochs=150, batch_size=128, validation_split=0.1,
                        callbacks=[EarlyStopping(patience=10), ModelCheckpoint(filepath=train_dir + 'weight.keras', save_best_only=True)],
                        shuffle=True, verbose=2)

    model.load_weights(train_dir + 'weight.keras')
    result = model.predict(test_data, verbose=0)
    np.savetxt(train_dir + 'test_label.txt', test_label)
    np.savetxt(train_dir + 'predict_result.txt', result)
    np.savetxt(train_dir + 'test_chrom.txt', test_chrom)
    np.savetxt(train_dir + 'test_region.txt', test_region)

    result_all = model.predict(data[0:(train_num + test_num)], verbose=0)
    #np.savetxt(train_dir + 'data_lstm.txt', data_lstm)
    np.savetxt(train_dir + 'label_all.txt', label)
    np.savetxt(train_dir + 'chrom_all.txt', chrom)
    np.savetxt(train_dir + 'region_all.txt', region)
    np.savetxt(train_dir + 'result_all.txt', result_all)



