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
import FunctionTwo
import FunctionThree
from concurrent.futures import ProcessPoolExecutor, as_completed
from joblib import dump, load

data = np.loadtxt('Stacking_AP2.txt',delimiter = ',',dtype = float)
print(data.shape)
X,y = np.split(data,(3,),axis=1)
y = y.ravel()
# 加载模型
model_Stacking = load('ModelOfStacking.joblib')
y_pred = model_Stacking.predict(X)  
y_pred_1 = model_Stacking.predict_proba(X)[:,1]
print(y_pred)
print(y_pred_1)
np.savetxt('Stacking_cotton_result.txt',y_pred_1, delimiter='\n')
