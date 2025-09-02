#!/usr/bin/python
with open('NegativeResult2.txt','r') as f_n, open('PositiveResult2.txt','r') as f_p:
    lines1 = f_n.readlines()
    lines2 = f_p.readlines()
seq_p = []
seq_n = []
for line1 in lines1:
    line1 = line1.strip().split(',')
    line1 = line1[1]
    seq_n.append(line1)
for line2 in lines2:
    line2 = line2.strip().split(',')
    line2 = line2[1]
    seq_p.append(line2)
seq_n=set(seq_n)
print(len(seq_n))
seq_p=set(seq_p)
print(len(seq_p))
intersection = seq_n.intersection(seq_p)
intersection = list(intersection)      #获得阳性和阴性数据交集列表
print(len(intersection))


#去掉交集数据
with open('NegativeResult2.txt','r') as f_n:
    lines1 = f_n.readlines()
with open('NegativeResult2_shaixuan.txt','w') as f1:
    for line1 in lines1:
        column1 = line1.strip().split(',')[1]
        column2 = line1.strip().split(',')[2]
        column3 = line1.strip().split(',')[3]
        column4 = line1.strip().split(',')[4]
        column5 = line1.strip().split(',')[5]
        if column1 not in intersection and float(column5) < 0.9 and float(column4) > 0.9 :
            f1.writelines(line1)
with open('NegativeResult2_shaixuan.txt','r') as f1:
    lines1 = f1.readlines()
print(len(lines1))

with open('PositiveResult2.txt','r') as f_n:
    lines1 = f_n.readlines()
with open('PositiveResult2_shaixuan.txt','w') as f1:
    for line1 in lines1:
        column1 = line1.strip().split(',')[1]
        column2 = line1.strip().split(',')[2]
        column3 = line1.strip().split(',')[3]
        column4 = line1.strip().split(',')[4]
        column5 = line1.strip().split(',')[5]
        if column1 not in intersection and float(column5) < 0.9 and float(column4) > 0.9 :
            f1.writelines(line1)
with open('PositiveResult2_shaixuan.txt','r') as f1:
    lines1 = f1.readlines()
print(len(lines1))

