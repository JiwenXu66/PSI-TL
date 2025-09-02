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
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc   # 用于分类问题评估  
from joblib import dump, load
import FunctionOne
import FunctionZero
from sklearn.model_selection import KFold
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

# 搭建随机森林模型
model_RF = RandomForestClassifier(n_estimators=100, random_state=42)
 
# 设置5倍交叉验证
kf = KFold(n_splits=5, shuffle=True, random_state=42)
 
# 存储每次交叉验证的结果
accuracy_scores = []
precision_scores = []
recall_scores = []
f1_scores = []
auc_scores = []
 
# 进行交叉验证
for train_index, test_index in kf.split(X):
    X_train_fold, X_test_fold = X[train_index], X[test_index]
    y_train_fold, y_test_fold = y[train_index], y[test_index]
    
    # 训练模型
    model_RF.fit(X_train_fold, y_train_fold)
    
    # 预测
    y_pred = model_RF.predict(X_test_fold)
    y_pred_proba = model_RF.predict_proba(X_test_fold)[:, 1]
    
    # 计算准确率
    accuracy = accuracy_score(y_test_fold, y_pred)
    accuracy_scores.append(accuracy)
    
    # 计算精确率（需要指定正类的标签，这里假设正类标签为1）
    precision = precision_score(y_test_fold, y_pred, pos_label=1)
    precision_scores.append(precision)
    
    # 计算召回率（需要指定正类的标签，这里假设正类标签为1）
    recall = recall_score(y_test_fold, y_pred, pos_label=1)
    recall_scores.append(recall)
    
    # 计算F1分数（需要指定正类的标签，这里假设正类标签为1）
    f1 = f1_score(y_test_fold, y_pred, pos_label=1)
    f1_scores.append(f1)
    
    # 计算ROC曲线和AUC
    fpr, tpr, thresholds = roc_curve(y_test_fold, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    auc_scores.append(roc_auc)
 
# 输出每一次验证的结果
for i in range(5):
    print(f"Fold {i+1}:")
    print(f"  Accuracy: {accuracy_scores[i]}")
    print(f"  Precision: {precision_scores[i]}")
    print(f"  Recall: {recall_scores[i]}")
    print(f"  F1 Score: {f1_scores[i]}")
    print(f"  AUC: {auc_scores[i]}")
    print()
 