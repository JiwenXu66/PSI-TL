#关于该函数文件的内容介绍
#1.similar_upstream（上游相似性）,similar_core（核心相似性）,calculate_AA_similarity（氨基酸的匹配率）,char_similarity（全长匹配率）
#1.四个函数用来对MSA获得的序列进行分析，返回一个“大于0，小于1”的数值
#2.seq_analyse函数调用上面四个函数，返回5个计算值和一个源数据中的p-value值
#3.SeqSearching,seq_select,seq_treatment主要用来对以上数据进行筛选，过滤掉MSA之后的一部分数据，从而获得初始版训练数据集
def similar_upstream(s1, s2):  #这个函数用于计算upstream区域的相似度
    import difflib
    return difflib.SequenceMatcher(None, s1, s2).quick_ratio()
def similar_core(s1, s2):  #这个函数用于计算core区域的相似度
    import difflib
    return difflib.SequenceMatcher(None, s1, s2).quick_ratio()
def calculate_AA_similarity(s1,s2):
    if len(s1) != len(s2):
        return 0
    matches = sum(1 for a, b in zip(s1,s2) if a==b)
    similarity = matches/ len(s1)
    return similarity
def char_similarity(str1,str2):
    if len(str1) != len(str2):
        return 0
    same_count = 0
    for i in range(len(str1)):
        if str1[i] == str2[i]:
            same_count += 1
    similarity = same_count / len(str1)
    return similarity
def seq_analyse(experi_seq,result_searching):
    import difflib
    with open('seq.txt','r',encoding='utf-8') as file:
        next(file)
        second_line = next(file).strip()
        if second_line:
            num_chars_in_seq = len(second_line) #氨基酸的长度
            positions_L = [i for i, char in enumerate(second_line) if char == 'L'] #核心基序L出现的位置
            AA_core = [second_line[idx] for idx in positions_L]
            upstream = second_line[:positions_L[0]]
            core = second_line[positions_L[0]:]
            #print(num_chars_in_seq,positions_L,AA_core,upstream,core)          #
            a = num_chars_in_seq - 1
            if a > positions_L[-1]:
                positions_bwteen_L = [x + 1 for x in positions_L]
                AA_softcore = [second_line[idx] for idx in positions_bwteen_L]  
                #print(positions_bwteen_L,AA_softcore) 
            if a == positions_L[-1]:
                positions_bwteen_L = [x + 1 for x in positions_L]
                positions_bwteen_L = positions_bwteen_L[:-1]
                AA_softcore = [second_line[idx] for idx in positions_bwteen_L]
    with open(result_searching,'r',encoding='utf-8') as result, open('Analys_result.txt','w') as outfile:
        next(result) #丢弃第一行
        for line in result:
            columns = line.strip().split() #去除行尾换行符，并按照空格分成列
            if len(columns) >= 9:
                nine_column = columns[8]
                AA_core_new = [nine_column[idx] for idx in positions_L]
                #print(nine_column,positions_L)
                a = num_chars_in_seq - 1
                if a > positions_L[-1]:
                    positions_bwteen_L = [x + 1 for x in positions_L]
                    #print(nine_column,positions_bwteen_L)
                    AA_softcore_new = [nine_column[idx] for idx in positions_bwteen_L]        
                if a == positions_L[-1]:
                    positions_bwteen_L = [x + 1 for x in positions_L]
                    positions_bwteen_L = positions_bwteen_L[:-1]
                    AA_softcore_new = [nine_column[idx] for idx in positions_bwteen_L]
                data_upstream = similar_upstream(upstream,nine_column[:positions_L[0]]) #计算上游相似度
                data_core = similar_core(core,nine_column[positions_L[0]:]) #计算核心部位相似度
                data_AA_core = calculate_AA_similarity(AA_core,AA_core_new)
                data_AA_softcore =  calculate_AA_similarity(AA_softcore,AA_softcore_new)
                data_AA_seq = char_similarity(second_line,nine_column)
                outfile.write(str(columns[1]) +','+ str(columns[8]) +','+ str(data_upstream) +','+ str(data_core) +','+ str(data_AA_core) +','+ str(data_AA_softcore) +','+ str(data_AA_seq) +','+ str(columns[5]) + '\n')


def SeqSearching(experi_seq,name_file):
    with open(experi_seq,'r',encoding='utf-8') as seq, open(name_file,'w',encoding='utf-8') as name:
        line_num = 0
        for line in seq:
           if line_num % 2 == 0:
             name.write(line)       
           line_num += 1
    with open(name_file,'r',encoding='utf-8') as file_in, open('ProteinName.txt','w',encoding='utf-8') as file_out:
        for line in file_in:
            line = line.lstrip('>')
            file_out.write(line)
def seq_select(inputfile1): #此步结果用于湿实验验证,超过80%互作
    import os
    import pandas as pd
    target_folder = 'WetTest'
    target_file = os.path.join(target_folder,'WetTest.txt')
    target_file2 = os.path.join(target_folder,'WetTestDeletDuplication.txt')
    os.system('mkdir WetTest')
    df = pd.read_csv(inputfile1,header=None,sep=',')
    for col in range(2,8):
        df.iloc[:,col] = df.iloc[:,col].astype(float)  #将所有数据修改为浮点数
    filtered_df = df[(df.iloc[:,2]>0.5) & (df.iloc[:,3]<0.6) & (df.iloc[:,4]>0.4) & (df.iloc[:,5]<0.6) & (df.iloc[:,6]<0.6) & (df.iloc[:,7]>1e-10)] #这一步是筛选的关键！！！！！
    filtered_df.to_csv(target_file, index=False, header=False)   
    with open(target_file,'r',encoding='utf-8') as file:
        lines = file.readlines()
    modified_lines = []
    for line in lines:
        columns = line.strip().split(',')
        modified_line = columns[0] + '\n'
        modified_lines.append(modified_line)
    with open(target_file2,'w',encoding='utf-8') as file1:
        file1.writelines(modified_lines)
def seq_treatment(file1,file2): #对file1 和file2 进行去重处理，均根据基序所在序列进行去重
    import pandas as pd  
    with open(file1,'r') as f1:
        lines = f1.readlines()
    with open('linshi.txt','w') as f2:
        for line in lines:
            line = line.strip().split()
            if "p-value" not in line:
                f2.write(line[0] +'\t'+ line[1] +'\t'+ line[2] +'\t'+ line[3] +'\t'+ line[4] +'\t'+ line[5] +'\t'+ line[6] +'\t'+ line[7] +'\t'+ line[8] + '\n')
            
    df = pd.read_csv('linshi.txt', header=None, sep = '\t', names=[f'Col{i+1}' for i in range(9)])
    rows,cols = df.shape
    print(rows,cols)
    df['Col6'] = pd.to_numeric(df['Col6'], errors='coerce') 
    idx_to_keep = df.groupby(df.iloc[:, 8])['Col6'].idxmax()  #对同一个序列，根据第七列p值进行筛选，保留最大的那一个
    filtered_df = df.loc[idx_to_keep] 
    filtered_df.to_csv(file1, index=False, header=False,sep = ',')

    df = pd.read_csv(file2,header=None, sep = ',')
    df = pd.DataFrame(df)
    df_unique = df.drop_duplicates(subset=1, keep='first')  #删除第一列中出现的重复项，保留第一个出现的数据所在行
    df_unique.to_csv(file2, index=False, header=False,sep = ',')