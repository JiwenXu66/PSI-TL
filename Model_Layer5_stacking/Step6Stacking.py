#!/usr/bin/python3
import subprocess #调用linux上的sed命令
import os #调用linux上的其他命令
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
from sklearn.metrics import accuracy_score, roc_curve, f1_score, recall_score, precision_score, auc   # 用于分类问题评估  
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

##进行MSA
#import Step1MSAforData
##对MSA之后的数据进行筛选
#import Step2CleanData


##RF模型对待测数据进行打分
#对原始的阴性数据和阳性数据进行分析：数据特征是‘基因名’+‘\n’+‘序列’

#FunctionThree.remove_lines('Positive_ProteinSeq.txt','Positive_Test.txt','>')
#FunctionThree.remove_lines('Negative_ProteinSeq.txt','Negative_Test.txt','>')


#导入数据
data = np.loadtxt('RF_ProteinCharacterFile.txt',delimiter = ',',dtype = float)  
data1 = data[:,0:20]
data2 = data[:,420:564]
data3 = data[:,708:]
print(data1.shape,data2.shape,data3.shape)
data = np.concatenate((data1, data2,data3), axis=1)
X,y = np.split(data,(164,),axis=1)
print(X.shape,y.shape)
y = y.ravel()
#######随机森林部分########
model_RF = load('ModelOfRandomForest.joblib')
y_pred_RF = model_RF.predict_proba(X)[:, 1] #接近阳性的概率
y_pred_RF = np.expand_dims(y_pred_RF, axis=1) 
print(y_pred_RF)

########神经网络部分#############
#import Test_CNN  #对待测数据提取特征矩阵
#读取数据
FunctionTwo.deal_folder('PositiveResult')
FunctionTwo.deal_folder('NegativeResult')
positive_data, positive_labels = FunctionTwo.read_folder('PositiveResult', 1)  
negative_data, negative_labels = FunctionTwo.read_folder('NegativeResult', 0)  
print("完成了数据的读取")
#合并数据  
X = np.concatenate((positive_data[:], negative_data[:]), axis=0)  
y = np.concatenate((positive_labels[:], negative_labels[:]), axis=0)  
print(positive_data.shape)
print(negative_data.shape)
X = X.astype(np.float32) / 10

########卷积神经网络###################
model_CNN = load_model('ModelOfCNN.keras')
y_pred_CNN = model_CNN.predict(X)   #预测概率
#print(y_pred_CNN)

##########Transformer-encoder部分###############
from Step5Transformer import model_Transformer
##用Transformer模型预测
positive_data, positive_labels = read_files('PositiveResult2_shaixuan.txt','1')  
negative_data, negative_labels = read_files('NegativeResult2_shaixuan.txt','0') 
test_inputs = tf.concat((positive_data, negative_data), axis=0)
test_labels = tf.concat((positive_labels, negative_labels), axis=0)
#print(test_inputs.shape)
y_pred_Trans = model_Transformer.predict(test_inputs,batch_size=64)

print(y_pred_RF.shape)  
print(y_pred_CNN.shape)  
print(y_pred_Trans.shape)  
y = y.reshape(-1, 1) #给y增加一个维度

pre_result = np.concatenate((y_pred_RF, y_pred_CNN, y_pred_Trans), axis=1)
result = np.concatenate((pre_result,y),axis = 1)
#print(result)
#print(result.shape)

np.savetxt('Stacking.txt', result, delimiter=',') 


data = np.loadtxt('Stacking.txt',delimiter = ',',dtype = float)
print(data.shape)
X,y = np.split(data,(3,),axis=1)
y = y.ravel()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=42)
model_Stacking = RandomForestClassifier(n_estimators=100, random_state=42)  
#训练模型
model_Stacking.fit(X_train, y_train) 
#保存模型
dump(model_Stacking, 'ModelOfStacking.joblib')
#model_Stacking = load('ModelOfStacking.joblib')
y_pred = model_Stacking.predict(X_test)  
#y_pred_1 = model_Stacking.predict_proba(X)[:,1]
#print(y_pred)
#print(y_pred_1)
#np.savetxt('result.txt',y_pred_1, delimiter='\n')
#检测
accuracy = accuracy_score(y_test, y_pred)  
print(f"Accuracy: {accuracy}")  
#计算精确率（需要指定正类的标签，这里假设正类标签为1）  
precision = precision_score(y_test, y_pred, pos_label=1)  
print(f'Precision: {precision}')  
#计算召回率（需要指定正类的标签，这里假设正类标签为1）  
recall = recall_score(y_test, y_pred, pos_label=1)  
print(f'Recall: {recall}')   
#计算F1分数（需要指定正类的标签，这里假设正类标签为1）  
f1 = f1_score(y_test, y_pred, pos_label=1)  
print(f'F1 Score: {f1}')  