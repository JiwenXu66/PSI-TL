#!/usr/bin/python
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
from sklearn.linear_model import LogisticRegression 
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc   # 用于分类问题评估  
from joblib import dump, load
import FunctionOne
import FunctionZero
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict, cross_val_score

print("开始构建随机森林模型")
#########随机森林部分########
#随机森林-提取阳性数据的数据特征
#modified_lines_1 = []; modified_lines_2 = []; modified_lines_3 = []
#inputfile = 'PositiveResult2_shaixuan.txt'
#FunctionOne.AA_composition(inputfile,modified_lines_1)
#FunctionOne.dipeptide_composition(inputfile,modified_lines_2)
#FunctionOne.Character_select(inputfile,modified_lines_3,"1")
#FunctionZero.seq_select(inputfile)
#with open('RF_ProteinCharacterOfPositive.txt','w',encoding='utf-8') as f:
#    for column1, column2, column3 in zip(modified_lines_1, modified_lines_2, modified_lines_3):
#        f.write(f"{column1},{column2},{column3}\n")
#os.system('mkdir PositiveResult')
#os.system('mv WetTest PositiveResult/')
#随机森林-提取阴性数据的数据特征
#modified_lines_1 = []; modified_lines_2 = []; modified_lines_3 = []
#inputfile = 'NegativeResult2_shaixuan.txt'
#FunctionOne.AA_composition(inputfile,modified_lines_1)
#FunctionOne.dipeptide_composition(inputfile,modified_lines_2)
#FunctionOne.Character_select(inputfile,modified_lines_3,"0")
#FunctionZero.seq_select(inputfile)
#with open('RF_ProteinCharacterOfNegative.txt','w',encoding='utf-8') as f:
#    for column1, column2, column3 in zip(modified_lines_1, modified_lines_2, modified_lines_3):
#        f.write(f"{column1},{column2},{column3}\n")
#os.system('mkdir NegativeResult')
#os.system('mv WetTest NegativeResult/')
#随机森林-阴性和阳性数据整合  
#with open('RF_ProteinCharacterOfNegative.txt','r') as f1, open('RF_ProteinCharacterOfPositive.txt','r') as f2:
#    lines1 = f1.readlines()[:]
#    lines2 = f2.readlines()[:]
#    with open('RF_ProteinCharacterFile.txt','w') as f3:
#        f3.writelines(lines1)
#        f3.writelines(lines2)

#os.system('cat RF_ProteinCharacterOfPositive.txt RF_ProteinCharacterOfNegative.txt > RF_ProteinCharacterFile.txt') #linux 上
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


# 划分测试集
#X_train_full, X_test, y_train_full, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
 
# 设置逻辑回归模型
model_LR = LogisticRegression(solver='saga', max_iter=500, random_state=42)
 
# 使用StratifiedKFold进行5倍交叉验证
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
 
# 存储每次验证的结果
accuracies = []
precisions = []
recalls = []
f1_scores = []
aucs = []
 
for train_index, val_index in skf.split(X, y):
    X_train, X_val = X [train_index], X [val_index]
    y_train, y_val = y [train_index], y [val_index]
    
    # 训练模型
    model_LR.fit(X_train, y_train)
    
    # 预测验证集
    y_val_pred = model_LR.predict(X_val)
    y_val_pred_proba = model_LR.predict_proba(X_val)[:, 1]
    
    # 计算准确率
    accuracy = accuracy_score(y_val, y_val_pred)
    accuracies.append(accuracy)
    
    # 计算精确率
    precision = precision_score(y_val, y_val_pred, pos_label=1)
    precisions.append(precision)
    
    # 计算召回率
    recall = recall_score(y_val, y_val_pred, pos_label=1)
    recalls.append(recall)
    
    # 计算F1分数
    f1 = f1_score(y_val, y_val_pred, pos_label=1)
    f1_scores.append(f1)
    
    # 计算ROC曲线和AUC
    fpr, tpr, thresholds = roc_curve(y_val, y_val_pred_proba)
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
 
# 输出每一次验证的结果
for i in range(5):
    print(f"Fold {i+1} - Accuracy: {accuracies[i]}, Precision: {precisions[i]}, Recall: {recalls[i]}, F1 Score: {f1_scores[i]}, AUC: {aucs[i]}")
 
# 计算所有验证结果的平均值
mean_accuracy = np.mean(accuracies)
mean_precision = np.mean(precisions)
mean_recall = np.mean(recalls)
mean_f1 = np.mean(f1_scores)
mean_auc = np.mean(aucs)
 
print(f"Mean Accuracy: {mean_accuracy}")
print(f"Mean Precision: {mean_precision}")
print(f"Mean Recall: {mean_recall}")
print(f"Mean F1 Score: {mean_f1}")
print(f"Mean AUC: {mean_auc}")
 










