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
from sklearn.model_selection import train_test_split  
from sklearn.preprocessing import LabelEncoder  
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_score, recall_score, f1_score 
import matplotlib.pyplot as plt 
#from rdkit import Che
#from rdkit.Chem import AllChem
import tensorflow as tf   
from tensorflow.keras.layers import Layer, Embedding, LayerNormalization, MultiHeadAttention, Dense, GlobalAveragePooling1D
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dropout  
from tensorflow.keras.initializers import TruncatedNormal  #初始化器，避免梯度消失或爆炸的问题
from tensorflow.keras.models import Model, Sequential, load_model   
import FunctionTwo
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.model_selection import StratifiedKFold
from keras.optimizers import Adam
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

# 数据归一化
X = X.astype(np.float32) / 5
y = y.astype(np.int32)
 
# 5折交叉验证
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
fold_accuracies = []
fold_precisions = []
fold_recalls = []
fold_f1s = []
fold_aucs = []
 
for train_index, val_index in skf.split(X, y):
    X_train, X_val = X[train_index], X[val_index]
    y_train, y_val = y[train_index], y[val_index]
 
    # 构建CNN模型
    model_CNN = Sequential([
        Conv1D(4, 2, activation='relu', input_shape=(8, 58), padding='SAME'),
        MaxPooling1D(2),
        Conv1D(8, 2, activation='relu', padding='SAME'),
        MaxPooling1D(2),
        Flatten(),
        Dense(64, activation='relu'),
        Dropout(0.5),
        Dense(1, activation='sigmoid')
    ])
 
    model_CNN.compile(optimizer=Adam(), loss='binary_crossentropy', metrics=['accuracy'])
 
    # 训练模型
    model_CNN.fit(X_train, y_train, epochs=10, batch_size=32, verbose=0)
 
    # 预测验证集
    y_val_pred = (model_CNN.predict(X_val) > 0.5).astype("int32").flatten()
    y_val_pred_proba = model_CNN.predict(X_val).flatten()
 
    # 计算指标
    accuracy = accuracy_score(y_val, y_val_pred)
    precision = precision_score(y_val, y_val_pred)
    recall = recall_score(y_val, y_val_pred)
    f1 = f1_score(y_val, y_val_pred)
    fpr, tpr, thresholds = roc_curve(y_val, y_val_pred_proba)
    roc_auc = auc(fpr, tpr)
 
    # 保存结果
    fold_accuracies.append(accuracy)
    fold_precisions.append(precision)
    fold_recalls.append(recall)
    fold_f1s.append(f1)
    fold_aucs.append(roc_auc)
 
    print(f"Fold Accuracy: {accuracy}")
    print(f"Fold Precision: {precision}")
    print(f"Fold Recall: {recall}")
    print(f"Fold F1 Score: {f1}")
    print(f"Fold AUC: {roc_auc}")
 
