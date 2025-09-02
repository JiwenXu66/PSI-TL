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
from sklearn.metrics import accuracy_score  # 用于分类问题评估  
from rdkit import Chem
from rdkit.Chem import AllChem
import tensorflow as tf   
from tensorflow.keras.layers import Layer, Embedding, LayerNormalization, MultiHeadAttention, Dense, GlobalAveragePooling1D
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dropout  
from tensorflow.keras.initializers import TruncatedNormal  #初始化器，避免梯度消失或爆炸的问题
from tensorflow.keras.models import Model, Sequential, load_model   
import FunctionOne
import FunctionTwo
import FunctionZero
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
os.system('mkdir ScoringMatrix')
os.system('mkdir SecondStructureAnnotation')
os.system('mkdir PutativeRSA')
os.system('mkdir 520AP2_Result')
#创建pssm搜索数据库
with open('PSSM.data.txt','r') as f:
    lines = f.readlines()
with open('output.fasta','w') as f1:
    for line in lines:
        column0 = line.strip().split(',')[0]
        column1 = line.strip().split(',')[1]
        f1.write('>' + column0 + '\n' + column1 + '\n')
os.system('makeblastdb -dbtype prot -in output.fasta -input_type fasta -out test.blastdb') 
inputfile2 = '520AP2_seq.txt'

with open(inputfile2,'r',encoding='utf-8') as file:
    lines = file.readlines()[579:]   
file_index = 0
for line in lines:  
    modified_lines_1 = []; modified_lines_2 = []; modified_lines_3 = []
    modified_lines_4 = []; modified_lines_5 = []; modified_lines_6 = []     
    columns = line.strip().split(',')
    if len(columns) > 0:
        seq = str(columns[1])
        print ("############################提取出来的序列是",seq)
        FunctionTwo.AA_one_hot(seq,modified_lines_1)
        FunctionTwo.Physicochemical_properties(seq,modified_lines_2)
        FunctionTwo.Physicol_properties(seq,modified_lines_3)
        FunctionTwo.scoring_matrix_2(seq,modified_lines_4,file_index,columns)
        FunctionTwo.second_structure_annotation(modified_lines_5,file_index,columns)
        FunctionTwo.Putative_RSA(modified_lines_6,columns,file_index)
    target_folder = '520AP2_Result'
    output_filename = f'Seq_{file_index:02d}.txt'
    targetfile = os.path.join(target_folder,output_filename) 
    print ("#######################################################################################3")
    print (len(modified_lines_1), len(modified_lines_2), len(modified_lines_3), len(modified_lines_4), len(modified_lines_5), len(modified_lines_6))
    with open(targetfile,'w',encoding='utf-8') as f:
        for column1, column2, column3, column4, column5, column6 in zip(modified_lines_1, modified_lines_2, modified_lines_3, modified_lines_4, modified_lines_5, modified_lines_6):
            f.write(f"{column1},{column2},{column3},{column4},{column5},{column6}")
    file_index += 1
os.system('mv ScoringMatrix 520AP2_Result/')
os.system('mv SecondStructureAnnotation 520AP2_Result/')
os.system('mv PutativeRSA 520AP2_Result/')
file_index = 0
for i in range(0,len(lines)):
    folderpath =  f'asaq.output_{file_index:02d}.txt'
    folderpath1 = f'ASA_{file_index:02d}.txt'
    folderpath2 = f'result1_{file_index:02d}.txt'
    folderpath3 = f'result2_{file_index:02d}.txt'
    shutil.rmtree(folderpath)
    os.remove(folderpath1)
    os.remove(folderpath2)
    os.remove(folderpath3)
    file_index += 1     
