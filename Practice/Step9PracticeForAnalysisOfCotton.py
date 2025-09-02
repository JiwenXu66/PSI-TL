#
import subprocess #调用linux上的sed命令
import os #调用linux上的其他命令
import shutil
import argparse
import difflib   #比较相似度函数
import pandas as pd
import numpy as np
import re
import random #用来随机抽数
import shutil#用来copy文件夹
from sklearn.ensemble import RandomForestClassifier #用于分析分类问题
from sklearn.model_selection import train_test_split  
from sklearn.preprocessing import LabelEncoder  
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_score, recall_score, f1_score 
import matplotlib.pyplot as plt 
#from rdkit import Chem
#from rdkit.Chem import AllChem
import tensorflow as tf   
from tensorflow.keras.layers import Layer, Embedding, LayerNormalization, MultiHeadAttention, Dense, GlobalAveragePooling1D
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dropout  
from tensorflow.keras.initializers import TruncatedNormal  #初始化器，避免梯度消失或爆炸的问题
from tensorflow.keras.models import Model, Sequential, load_model   
import FunctionZero
import FunctionOne
import FunctionTwo
import FunctionThree
from joblib import dump, load  #用来储存和读取随机森林模型的
from tensorflow.keras.utils import register_keras_serializable
from keras.utils import custom_object_scope  

def read_files(filename,label):  
    data = []  
    labels = [] 
    AA_index = {'A': 0,  # Alanine  
                'R': 1,  # Arginine  
                'N': 2,  # Asparagine  
                'D': 3,  # Aspartic acid  
                'C': 4,  # Cysteine  
                'Q': 5,  # Glutamine  
                'E': 6,  # Glutamic acid  
                'G': 7,  # Glycine  
                'H': 8,  # Histidine  
                'I': 9,  # Isoleucine  
                'L': 10, # Leucine  
                'K': 11, # Lysine  
                'M': 12, # Methionine  
                'F': 13, # Phenylalanine  
                'P': 14, # Proline  
                'S': 15, # Serine  
                'T': 16, # Threonine  
                'W': 17, # Tryptophan  
                'Y': 18, # Tyrosine  
                'V': 19  # Valine   
                }
    with open(filename,'r') as f:
        lines = f.readlines()
        #lines = set(lines)
    for line in lines:
        line = line.strip().split(',')
        if len(line) >0 and 'p-value' not in line:
            column = line[1]
            column = column[-8:]
            if all(aa in AA_index for aa in column):
                AA_indices = [ AA_index[aa] for aa in column]
                tensor = tf.constant(AA_indices,dtype=tf.int32)
                data.append(tensor)    
                labels.append(int(label))    
    return data,labels             

modified_lines_1 = []; modified_lines_2 = []; modified_lines_3 = []
FunctionOne.AA_composition('520AP2_seq.txt',modified_lines_1)
FunctionOne.dipeptide_composition('520AP2_seq.txt',modified_lines_2)
FunctionOne.Character_select('520AP2_seq.txt',modified_lines_3,"1")
# 随机森林部分
with open('RF_520AP2_seq.txt','w',encoding='utf-8') as f:
    for column1, column2, column3 in zip(modified_lines_1, modified_lines_2, modified_lines_3):
        f.write(f"{column1},{column2},{column3}\n")
data = np.loadtxt('RF_520AP2_seq.txt',delimiter = ',',dtype = float)
data1 = data[:,0:20]
data2 = data[:,420:564]
data3 = data[:,708:]
print(data1.shape,data2.shape,data3.shape)
data = np.concatenate((data1, data2, data3), axis=1)
X,y_lable = np.split(data,(164,),axis=1)


########随机森林部分################
model_RF = load('ModelOfRandomForest.joblib')
y_pred_RF = model_RF.predict_proba(X)[:, 1] #接近阳性的概率
y_pred_RF = np.expand_dims(y_pred_RF, axis=1) 
print(y_pred_RF)

#神经网络部分##########
#import Test_CNN  #对待测数据提取特征矩阵
#读取数据
FunctionTwo.deal_folder('520AP2_Result')
positive_data, positive_labels = FunctionTwo.read_folder('520AP2_Result', 1)  
negative_data, negative_labels = FunctionTwo.read_folder('520AP2_Result', 0) 
print("完成了数据的读取")
#合并数据  
X = np.concatenate((positive_data, negative_data), axis=0)  
y = np.concatenate((positive_labels, negative_labels), axis=0)

X = X[0:1379]
X = X.astype(np.float32) / 10

########卷积神经网络部分#############
model_CNN = load_model('ModelOfCNN.keras')
y_pred_CNN = model_CNN.predict(X)   #预测概率
print(X.shape)
print(y_pred_CNN)





##########Transformer-encoder部分###############
from Step5Transformer import model_Transformer
#用Transformer模型预测
positive_data, positive_labels = read_files('520AP2_seq.txt','1')  
negative_data, negative_labels = read_files('520AP2_seq.txt','0') 
test_inputs = tf.concat((positive_data, negative_data), axis=0)
test_inputs = test_inputs[:1379]
test_labels = tf.concat((positive_labels, negative_labels), axis=0)
print(test_inputs.shape)
y_pred_Trans = model_Transformer.predict(test_inputs,batch_size=64)

print(y_pred_RF.shape)  
print(y_pred_CNN.shape)  
print(y_pred_Trans.shape)

pre_result = np.concatenate((y_pred_RF, y_pred_CNN, y_pred_Trans), axis=1)
result = np.concatenate((pre_result,y_lable),axis = 1)
print(result)
print(result.shape)

np.savetxt('Stacking_AP2.txt', result, delimiter=',') 