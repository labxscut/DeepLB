"""Generate tumor-only random runs with purity filtering and data loading."""

import sys
import os

# import ConfigParser
import re
import numpy as np
import subprocess
import math


def readData(fileName):
    data = {}
    with open(fileName) as f:
        for line in f:
            line = re.sub("NA", "nan", line)
            info = line.rstrip().split("\t")
            index = info.pop(0)
            # data[index] = np.array(map(float,info))
            data[index] = info
    return data


def splitData(allMethy, trainRatio):
    # random splitting
    trainMethy = {}
    testMethy = {}
    sampleIndex = list(allMethy.keys())
    trainNum = int(len(sampleIndex) * trainRatio)
    trainIndex = set(np.random.choice(sampleIndex, size=trainNum, replace=False))
    for index in sampleIndex:
        if index in trainIndex:
            trainMethy[index] = allMethy[index]
        else:
            testMethy[index] = allMethy[index]
    return [trainMethy, testMethy]


def splitData_v2(allMethy, trainSampleNames, sampleIndex2Name):
    # retrive the splitting from previous setting
    trainMethy = {}
    testMethy = {}
    sampleIndex = allMethy.keys()
    for index in sampleIndex:
        if sampleIndex2Name[index] in trainSampleNames:
            trainMethy[index] = allMethy[index]
        else:
            testMethy[index] = allMethy[index]
    return [trainMethy, testMethy]


def readPairedData(fileName1, fileName2):
    pairedSampleList = []
    data1 = {}
    with open(fileName1) as f:
        for line in f:
            line = re.sub("NA", "nan", line)
            info = line.rstrip().split("\t")
            index = info.pop(0)
            data1[index] = info  # np.array(map(float,info))
            pairedSampleList.append([index, "nan"])
    data2 = {}
    i = 0
    with open(fileName2) as f:
        for line in f:
            line = re.sub("NA", "nan", line)
            info = line.rstrip().split("\t")
            index = info.pop(0)
            data2[index] = info  # np.array(map(float,info))
            pairedSampleList[i][1] = index
            i += 1
    return (data1, data2, pairedSampleList)


cancer2tissue = {
    "lihc": "liver",
    "luad": "lung",
    "lusc": "lung",
    "brca": "breast",
    "read": "rectum",
    "paad": "pancrea",
    "coad": "colon",
    "kich": "kidney",
    "kirc": "kidney",
    "kirp": "kidney",
    "prad": "prostate",
    "gbm": "brain",
    "blca": "bladder",
    "stad": "stomach",
    "ucec": "uterus",
    "ucs": "uterus",
    "npc": "naso",
}  # add naso carcinoma by QJ 6/14/2018


runID = sys.argv[1]
np.random.seed(20240220 + int(runID))
cancer = sys.argv[2]  # for example, 'lihc', 'brca', 'luad', 'lusc'
outDir = sys.argv[3]
print(outDir)
 #top30,top40,top50,bottom30,bottom40,bottom50
# trainNormalSampleNames = sys.argv[4].split(',')
rootDir = sys.argv[4]
SampleList = sys.argv[5]
thre_type = sys.argv[6] #top30
sampleAnnotFile = rootDir + "/Predata/metadata/" + SampleList
 

binNum = 36342  # 42374 #10531
normalTrainRatio = 1#0.75
normalType = "plasma_background"
tumorType = cancer + "_tumor"
tissueType = cancer2tissue[cancer]

# For seq data of Shuli's method
depthCut = 20  #

##############################################################
##### load Testing plasma data: Dennis 2013 and 2015 single-end seq data, zhou lab (methylation levels) #####
##############################################################
print(
    "load plasma seq data: Dennis Lo 2013 and 2015 (including longitudinal data) and zhou lab data..."
)
seqDir='/home/yinliang/PROJECT/mMTS/mMTS_result/1.3_cal_methylratio_bins/all_reads_1sites_with_pos_and_sam'
#seqDir = rootDir + "/Result/1.3_extract_reads/all_reads_1sites_with_pos_and_sam"
# only load plasma samples information from sample_index
# file format
# 1	TBR52	plasma	neuro	tumor
# 2	CTR101	plasma	background	normal
# ...
# 34	TBR51	plasma	sms	tumor
# 35	HOT151	plasma	lihc	tumor
seqFileType = {}
with open(sampleAnnotFile) as f:
    for line in f:
        info = line.rstrip().split("\t")
        if info[1] == "plasma":
            if info[2] not in seqFileType:
                seqFileType[info[2]] = []
            seqFileType[info[2]].append(info[0])
seqFileType_backgroud = {
    key: value for key, value in seqFileType.items() if key in ["background"]
}

realSamples = {}
for type in seqFileType_backgroud:
    if (type == "background") or (type == "normal"):
        typeStr = "plasma_background"
    else:
        typeStr = type + "_tumor"
    for sample in seqFileType[type]:
        print("\t", sample)
        fileName = sample + ".methy_bins"
        methyDepth = ["NA"] * binNum
        allDepth = ["NA"] * binNum
        # file format for single-end reads file:
        # marker_index	methylation_rate	total_count	methylation_count	unmethylation_count	#reads
        # 1	0	0	0	0	0
        # 2	0.722972972972973	148	107	41	72
        # 3	0.867256637168142	113	98	15	57
        with open(os.path.join(seqDir, fileName)) as f:
            # f.next() # skip the first header line
            next(
                f
            )  # File object in Python 3 doesn't support next() method. Python 3 has a built-in function next()
            for line in f:
                (
                    binIndex,
                    methyValue,
                    allCytosines,
                    methyCytosines,
                    unmethyCytosines,
                    allReads,
                ) = line.rstrip().split("\t")
                if int(allCytosines) < depthCut:
                    continue
                binIndex = int(binIndex) - 1  # from 1-based to 0-based
                methyDepth[binIndex] = methyCytosines
                allDepth[binIndex] = allCytosines
        readsFile = os.path.join(seqDir, sample + ".reads_bins.gz")
        realSamples[sample] = {
            "id": typeStr + "_real_" + sample,
            "methy": methyDepth,
            "all": allDepth,
            "file": readsFile,
        }


###############################################################
#####  array data and normal plasma (methylation levels)  #####
###############################################################
print("load tcga array data ...")
arrayDir_tumor = (
    rootDir + "/Result/1.3_cal_methylratio_bins/450K"
)  # this needs to be changed, when moving project directory to other places
arrayDir_normal_plasma = (
   rootDir + "/Result/1.3_cal_methylratio_bins/WGBS"
)

sampleName2Index = {}  # real blood sample name to index
sampleIndex2Name = {}  # real blood sample index to name
# file format
# 1	TBR52	plasma	neuro	tumor
# 2	CTR101	plasma	background	normal
# ...
# 34	TBR51	plasma	sms	tumor
# 35	HOT151	plasma	lihc	tumor
with open(sampleAnnotFile) as f:
    for ind,line in enumerate(f, start=1):
        info = line.rstrip().split("\t")
        # if info[2] == 'plasma':
        sampleName2Index[info[0]] = ind
        sampleIndex2Name[ind] = info[0]

methyFile1 = os.path.join(arrayDir_tumor, cancer+ "_" + thre_type + ".tsv")
print("\t", cancer.upper() + "_" + thre_type + ".tsv")
allMethyTumor = readData(methyFile1)


methyFile = os.path.join(arrayDir_normal_plasma, normalType+ "_all.tsv")
print("\t", normalType + "_all.tsv")
allMethyPlasma = readData(methyFile)

###############################################################
# generate train and test samples
#   train: all solid matched tumor/normal samples + 75% of normal plasma samples (normalTrainRatio=75%)
#   test:  25% of normal plasma samples and all other plasma real/simulated samples
###############################################################
print("split data to training and test partitions, and write to files ...")
trainMethy_tumor = allMethyTumor

[trainMethy_normal, testMethy_normal] = splitData(allMethyPlasma, normalTrainRatio)
print("Train nomal plasma:")
print(trainMethy_normal.keys())
print("Test nomal plasma:")
print(testMethy_normal.keys())

trainOut = open(os.path.join(outDir, "train_" + runID), "w")
# trainMatchedOut = open(os.path.join(outDir, 'train_'+runID+'.matched'),'w')
trainOut_SampleList = open(os.path.join(outDir, "train_" + runID + ".list"), "w")
# samplelist = list(trainMethy_tumor.keys())
# samplelist.sort()
# print(samplelist)

for sampleIndex in trainMethy_tumor:
    #    print(sampleIndex)
    trainOut.write(
        "\t".join(
            [tumorType]
            + [re.sub("nan", "NA", str(x)) for x in trainMethy_tumor[sampleIndex]]
        )
        + "\n"
    )
    trainOut_SampleList.write(sampleIndex + "\n")

# samplelist = trainMethy_normal.keys()
# samplelist.sort()
for sampleIndex in trainMethy_normal:
    trainOut.write(
        "\t".join(
            [normalType]
            + [re.sub("nan", "NA", str(x)) for x in trainMethy_normal[sampleIndex]]
        )
        + "\n"
    )
    trainOut_SampleList.write(sampleIndex + "\n")

#testMethyOut = open(os.path.join(outDir, "test_" + runID + ".methycnt"), "w")
#testAllOut = open(os.path.join(outDir, "test_" + runID + ".allcnt"), "w")
testOut_SampleList = open(os.path.join(outDir, "test_" + runID + ".list"), "w")
print(realSamples.keys())
for sample in testMethy_normal.keys():
    realID = realSamples[sample]["id"]
    useThisSample = False
    if sample not in sampleName2Index:
        useThisSample = True
    elif (
        sampleName2Index[sample] not in trainMethy_normal
    ):  # if this normal plasma is not used for training, put into testing list
        useThisSample = True
    if useThisSample:
        #testMethyOut.write("\t".join([realID] + realSamples[sample]["methy"]) + "\n")
        #testAllOut.write("\t".join([realID] + realSamples[sample]["all"]) + "\n")
        testOut_SampleList.write(realID + "\t" + realSamples[sample]["file"] + "\n")
# samplelist = simuSamples.keys()
# samplelist.sort()
# for sample in samplelist:
#    simuID = simuSamples[sample]['id']
#    testMethyOut.write('\t'.join([simuID]+simuSamples[sample]['methy'])+'\n')
#    testAllOut.write('\t'.join([simuID]+simuSamples[sample]['all'])+'\n')
#    testOut_SampleList.write(simuID+'\t'+simuSamples[sample]['file']+'\n')
print("done.")
