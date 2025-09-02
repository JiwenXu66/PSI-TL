#!/usr/bin/python
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
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import LabelEncoder  
from tensorflow.keras.layers import SimpleRNN, Dense, Input, GRU
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_score, recall_score, f1_score 
import matplotlib.pyplot as plt 
#from rdkit import Chem
#from rdkit.Chem import AllChem
import tensorflow as tf   
from tensorflow.keras.layers import Layer, Embedding, LayerNormalization, MultiHeadAttention, Dense, GlobalAveragePooling1D
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dropout  
from tensorflow.keras.initializers import TruncatedNormal  #初始化器，避免梯度消失或爆炸的问题
from tensorflow.keras.models import Model, Sequential, load_model   
import FunctionTwo
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_lines(lines,i,folder): 
    file_index = i
    for line in lines:  
        modified_lines_1 = []; modified_lines_2 = []; modified_lines_3 = []
        modified_lines_4 = []; modified_lines_5 = []; modified_lines_6 = []     
        columns = line.strip().split(',')
        if len(columns) > 1:
            seq = str(columns[1])
            #print ("############################提取出来的序列是",seq)
            FunctionTwo.AA_one_hot(seq,modified_lines_1)
            FunctionTwo.Physicochemical_properties(seq,modified_lines_2)
            FunctionTwo.Physicol_properties(seq,modified_lines_3)
            FunctionTwo.scoring_matrix_2(seq,modified_lines_4,file_index,columns)
            FunctionTwo.second_structure_annotation(modified_lines_5,file_index,columns)
            FunctionTwo.Putative_RSA(modified_lines_6,columns,file_index)
        target_folder = folder
        output_filename = f'Seq_{file_index:02d}.txt'
        targetfile = os.path.join(target_folder,output_filename) 
        #print ("#######################################################################################")
        with open('record.txt','a') as re:
            re.write(str(file_index) +','+ "性质：" +','+ str(len(modified_lines_1))+','+str(len(modified_lines_2))+','+str(len(modified_lines_3))+','+str(len(modified_lines_4))+','+str(len(modified_lines_5))+','+str(len(modified_lines_6))+'\n')
        with open(targetfile,'w',encoding='utf-8') as f:
            for column1, column2, column3, column4, column5, column6 in zip(modified_lines_1, modified_lines_2, modified_lines_3, modified_lines_4, modified_lines_5, modified_lines_6):
                f.write(f"{column1},{column2},{column3},{column4},{column5},{column6}")
        file_index += 1

########卷积神经网络部分#############
#os.system('mkdir ScoringMatrix')
#os.system('mkdir SecondStructureAnnotation')
#os.system('mkdir PutativeRSA')
##阳性数据分析
#inputfile2 = 'PositiveResult2.txt'
#FunctionTwo.scoring_matrix_1(inputfile2) #创建pssm搜索数据库
#with open(inputfile2,'r',encoding='utf-8') as file:
#    lines = file.readlines()[:]
#    #lines = lines[86:]   
##########################################################333    
#chunk_size = len(lines) // 3000
#folder = 'PositiveResult'
#with ProcessPoolExecutor(max_workers=28) as executor:  # 创建一个包含5个线程的线程池   
#    futures = [executor.submit(process_lines, lines[i:i+chunk_size],i,folder) 
#               for i in range(0,len(lines),chunk_size)]  
#    for future in as_completed(futures):  
#        try:  
#            # 获取任务的结果  
#            result = future.result()  
#            # 处理结果（如果需要的话）  
#            #print(f"任务完成，结果: {result}")  
#        except Exception as exc:  
#            # 处理异常  
#            print(f"任务执行出错: {exc}")
###################################################################3
#os.system('mv ScoringMatrix PositiveResult/')
#os.system('mv SecondStructureAnnotation PositiveResult/')
#os.system('mv PutativeRSA PositiveResult/')

##阴性数据分析
#os.system('mkdir ScoringMatrix')
#os.system('mkdir SecondStructureAnnotation')
#os.system('mkdir PutativeRSA')
#inputfile2 = 'NegativeResult2.txt'
#FunctionTwo.scoring_matrix_1(inputfile2) #创建pssm搜索数据库
#with open(inputfile2,'r',encoding='utf-8') as file:
#    lines = file.readlines()[:]   
###############################################################
#chunk_size = len(lines) // 3000
#folder = 'NegativeResult'
#with ProcessPoolExecutor(max_workers=28) as executor:  # 创建一个包含5个线程的线程池   
#    futures = [executor.submit(process_lines, lines[i:i+chunk_size],i,folder) 
#               for i in range(0,len(lines),chunk_size)]  
#    for future in as_completed(futures):  
#        try:  
#            # 获取任务的结果  
#            result = future.result()  
#            # 处理结果（如果需要的话）  
#            #print(f"任务完成，结果: {result}")  
#        except Exception as exc:  
#            # 处理异常  
#            print(f"任务执行出错: {exc}")
##################################################################
#file_index = 0
#for i in range(0,len(lines)):
#    folderpath =  f'asaq.output_{file_index:02d}.txt'
#    folderpath1 = f'ASA_{file_index:02d}.txt'
#    folderpath2 = f'result1_{file_index:02d}.txt'
#    folderpath3 = f'result2_{file_index:02d}.txt'
#    shutil.rmtree(folderpath)
#    os.remove(folderpath1)
#    os.remove(folderpath2)
#    os.remove(folderpath3)
#    file_index += 1      
#os.system('mv ScoringMatrix NegativeResult/')
#os.system('mv SecondStructureAnnotation NegativeResult/')
#os.system('mv PutativeRSA NegativeResult/')
#CNN模型训练
FunctionTwo.deal_folder('PositiveResult')
FunctionTwo.deal_folder('NegativeResult')
#读取数据  
positive_data, positive_labels = FunctionTwo.read_folder('PositiveResult', 1)  
negative_data, negative_labels = FunctionTwo.read_folder('NegativeResult', 0)  
print("完成了数据的读取")
#合并数据  
X = np.concatenate((positive_data[:], negative_data[:]), axis=0)  
y = np.concatenate((positive_labels[:], negative_labels[:]), axis=0)  


# 数据类型转换和缩放
X = X.astype(np.float32) / 5
y = y.astype(np.int32)
 
# 定义KFold交叉验证
kf = KFold(n_splits=5, shuffle=True, random_state=42)
 
# 定义RNN模型构建函数
def create_model():
    model = Sequential([
        GRU(units=64, input_shape=(8, 58), return_sequences=False),
        Dense(32, activation='relu'),
        Dropout(0.2),
        Dense(units=1, activation='sigmoid') 
        ])  
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    return model
 
# 存储每次交叉验证的结果
accuracies, precisions, recalls, f1s, aucs = [], [], [], [], []
 
# 交叉验证过程
for train_index, val_index in kf.split(X):
    X_train, X_val = X[train_index], X[val_index]
    y_train, y_val = y[train_index], y[val_index]
 
    model = create_model()
    model.fit(X_train, y_train, epochs=10, batch_size=32, verbose=0)  # verbose=0以减少输出
 
    # 预测
    y_val_pred = (model.predict(X_val) > 0.5).astype("int32").flatten()
    y_val_pred_proba = model.predict(X_val).flatten()
 
    # 计算性能指标
    accuracy = accuracy_score(y_val, y_val_pred)
    precision = precision_score(y_val, y_val_pred)
    recall = recall_score(y_val, y_val_pred)
    f1 = f1_score(y_val, y_val_pred)
    fpr, tpr, thresholds = roc_curve(y_val, y_val_pred_proba)
    roc_auc = auc(fpr, tpr)
 
    # 存储结果
    accuracies.append(accuracy)
    precisions.append(precision)
    recalls.append(recall)
    f1s.append(f1)
    aucs.append(roc_auc)
 
# 输出每次交叉验证的结果
for i in range(5):
    print(f"Fold {i+1}:")
    print(f"  Accuracy: {accuracies[i]}")
    print(f"  Precision: {precisions[i]}")
    print(f"  Recall: {recalls[i]}")
    print(f"  F1 Score: {f1s[i]}")
    print(f"  AUC: {aucs[i]}")
    print()
 
# 如果需要，也可以输出平均性能指标
print(f"Average Accuracy: {np.mean(accuracies)}")
print(f"Average Precision: {np.mean(precisions)}")
print(f"Average Recall: {np.mean(recalls)}")
print(f"Average F1 Score: {np.mean(f1s)}")
print(f"Average AUC: {np.mean(aucs)}")