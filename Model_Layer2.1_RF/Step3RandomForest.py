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
##os.system('mkdir PositiveResult')
##os.system('mv WetTest PositiveResult/')
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
##os.system('mkdir NegativeResult')
##os.system('mv WetTest NegativeResult/')
#随机森林-阴性和阳性数据整合  
with open('RF_ProteinCharacterOfPositive.txt','r') as f1, open('RF_ProteinCharacterOfNegative.txt','r') as f2:
    lines1 = f1.readlines()[:]
    lines2 = f2.readlines()[:]
    with open('RF_ProteinCharacterFile.txt','w') as f3:
        f3.writelines(lines1)
        f3.writelines(lines2)

#os.system('cat RF_ProteinCharacterOfPositive.txt RF_ProteinCharacterOfNegative.txt > RF_ProteinCharacterFile.txt') #linux 上
#导入数据
data = np.loadtxt('RF_ProteinCharacterFile.txt',delimiter = ',',dtype = float)  
data1 = data[:,0:20]
data2 = data[:,420:564]
data3 = data[:,708:]

#data1 = data[:,0:420]
#data2 = data[:,420:708]
#data3 = data[:,708:]

print(data1.shape,data2.shape,data3.shape)
data = np.concatenate((data1, data2,data3), axis=1)
X,y = np.split(data,(164,),axis=1)
print(X.shape,y.shape)
y = y.ravel()

#划分测试，验证，训练集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42) 
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42) 
#搭建随机森林模型
model_RF = RandomForestClassifier(n_estimators=100, random_state=42)  
#训练模型
model_RF.fit(X_train, y_train) 
#保存模型
dump(model_RF, 'ModelOfRandomForest.joblib')
#使用验证集
score_val = model_RF.score(X_val, y_val)  
score_test = model_RF.score(X_test, y_test)
print(f"验证集得分：{score_val}" + "   " +f"测试集得分：{score_test}")
#验证测试集的accurancy
y_pred = model_RF.predict(X_test)  
y_pred_proba = model_RF.predict_proba(X_test)[:, 1]
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

#计算ROC曲线和AUC  
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)  
roc_auc = auc(fpr, tpr)  
print(f'AUC: {roc_auc}')  

#import matplotlib.pyplot as plt    
#plt.figure()  
#plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)  
#plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')  
#plt.xlim([0.0, 1.0])  
#plt.ylim([0.0, 1.05])  
#plt.xlabel('False Positive Rate')  
#plt.ylabel('True Positive Rate')  
#plt.title('Receiver Operating Characteristic')  
#plt.legend(loc="lower right")  
#plt.show()

