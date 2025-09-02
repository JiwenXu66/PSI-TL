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
from rdkit import Chem
from rdkit.Chem import AllChem
import FunctionZero

print("函数导入完成")
experi_seq = 'Positive_ProteinSeq.txt'
namefile = 'Positive_ProteinName.txt'
result_searching = 'PositiveSearchingResult.txt'
FunctionZero.SeqSearching(experi_seq,namefile)
with open('ProteinName.txt','r',encoding='utf-8') as name:
    for line in name:
        with open('name.txt','w',encoding='utf-8') as linshifile:
            linshifile.write(line)
        os.system('seqkit grep -f name.txt Positive_ProteinSeq.txt > seq.txt')
        os.system('echo ">1" >> seq1.txt')
        os.system('tail -1 seq.txt >> seq1.txt')
        os.system('cat seq1.txt >> seq.txt')
        os.system('meme seq.txt -protein -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 20 -maxw 50 ')
        os.system('fimo --oc fimo_result --verbosity 2 --bgfile --nrdb-- --thresh 1.0E-6 meme.txt DatabaseOfPanGreenPlant.fa')
        os.system('cat fimo_result/fimo.txt > PositiveSearchingResult.txt')         
        FunctionZero.seq_analyse(experi_seq,result_searching) #产生一个名为 “Analys_result.txt” 分析结果
        os.system('rm -f seq.txt')
        os.system('rm -f seq1.txt')
        os.system('cat PositiveSearchingResult.txt >> PositiveResult1.txt')
        os.system('cat Analys_result.txt >> PositiveResult2.txt')
        os.system('rm -f Analys_result.txt')
FunctionZero.seq_treatment('PositiveResult1.txt','PositiveResult2.txt')  #去除掉重复的结果
print ("阳性序列分析结束")
experi_seq = 'Negative_ProteinSeq.txt'
namefile = 'Negative_ProteinName.txt'
result_searching = 'NegativeSearchingResult.txt'
FunctionZero.SeqSearching(experi_seq,namefile)
with open('ProteinName.txt','r',encoding='utf-8') as name:
    for line in name:
        with open('name.txt','w',encoding='utf-8') as linshifile:
            linshifile.write(line)
        os.system('seqkit grep -f name.txt Negative_ProteinSeq.txt > seq.txt')
        os.system('echo ">1" >> seq1.txt')
        os.system('tail -1 seq.txt >> seq1.txt')
        os.system('cat seq1.txt >> seq.txt')
        os.system('meme seq.txt -protein -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 20 -maxw 50 ')
        os.system('fimo --oc fimo_result --verbosity 2 --bgfile --nrdb-- --thresh 1.0E-6 meme.txt DatabaseOfPanGreenPlant.fa')
        os.system('cat fimo_result/fimo.txt > NegativeSearchingResult.txt')
        FunctionZero.seq_analyse(experi_seq,result_searching)
        os.system('rm -f seq.txt')
        os.system('rm -f seq1.txt')
        os.system('cat NegativeSearchingResult.txt >> NegativeResult1.txt')
        os.system('cat Analys_result.txt >> NegativeResult2.txt')
        os.system('rm -f Analys_result.txt')   
FunctionZero.seq_treatment('NegativeResult1.txt','NegativeResult2.txt')
print ("阴性序列分析结束")