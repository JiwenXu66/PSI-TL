#该函数文件用来卷积神经网络分析

def AA_one_hot(seq,modified_lines_1): #将氨基酸用one-hot表示，并用向量标注每个位置的氨基酸
    library =  {'A':(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'C':(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'D':(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'E':(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'F':(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'G':(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'H':(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                'I':(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                'K':(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                'L':(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),  
                'M':(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
                'N':(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
                'P':(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                'Q':(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
                'R':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                'S':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
                'T':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0), 
                'V':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
                'W':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
                'Y':(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)}  
    for AA in seq:
        if AA not in library:  
            raise ValueError(f"Invalid amino acid letter: {AA}") 
        one_hot = library[AA]  
        modified_line = one_hot
        modified_line_str = ','.join(str(item) for item in modified_line) 
        modified_lines_1.append(modified_line_str)
    return modified_lines_1
def Physicochemical_properties(seq,modified_lines_2): #原子数，静电荷数，潜在氢键
    library = {'A':(5,0,2), 'C':(6,0,2), 'D':(8,-1,4), 'E':(9,-1,4), 'F':(11,0,2), 
               'G':(4,0,2), 'H':(10,0,4), 'I':(8,0,2), 'K':(9,1,2), 'L':(8,0,2), 
               'M':(8,0,2), 'N':(8,0,4), 'P':(7,0,2), 'Q':(9,0,4), 'R':(11,1,4),
               'S':(6,0,4), 'T':(7,0,4), 'V':(7,0,2), 'W':(14,0,3), 'Y':(12,0,3)}
    for AA in seq:
        if AA not in library:  
            raise ValueError(f"Invalid amino acid letter: {AA}") 
        properties = library[AA]
        modified_line = properties
        modified_line_str = ','.join(str(item) for item in modified_line) 
        modified_lines_2.append(modified_line_str)
    return modified_lines_2
def Physicol_properties(seq,modified_lines_3): #六维顺序：graph shape index;Polarizability;volume;Hydrophobicity;Isoelectric point;Helix probability;Sheet probability 
    library = { 'A':(1.28,0.05,1.00,0.31,6.11,0.42,0.23),'G':(0.00,0.00,0.00,0.00,6.07,0.13,0.15), 'V':(3.67,0.14,3.00,1.22,6.02,0.27,0.49), 
                'L':(2.59,0.19,4.00,1.70,6.04,0.39,0.31),'I':(4.19,0.19,4.00,1.80,6.04,0.30,0.45),'F':(2.94,0.29,5.89,1.79,5.67,0.30,0.38),
                'Y':(2.94,0.30,6.47,0.96,5.66,0.25,0.41),'W':(3.21,0.41,8.08,2.25,5.94,0.32,0.42),'T':(3.03,0.11,2.60,0.26,5.60,0.21,0.36),
                'S':(1.31,0.06,1.60,-0.04,5.70,0.20,0.28),'R':(2.34,0.29,6.13,-1.01,10.74,0.36,0.25),'K':(1.89,0.22,4.77,-0.99,9.99,0.32,0.27),
                'H':(2.99,0.23,4.66,0.13,7.69,0.27,0.30),'D':(1.60,0.11,2.78,-0.77,2.95,0.25,0.20),'E':(1.56,0.15,3.78,-0.64,3.09,0.42,0.21),
                'N':(1.60,0.13,2.95,-0.60,6.52,0.21,0.22),'Q':(1.56,0.18,3.95,-0.22,5.65,0.36,0.25),'M':(2.35,0.22,4.43,1.23,5.71,0.38,0.32),
                'P':(2.67,0.00,2.72,0.72,6.80,0.13,0.34),'C':(1.77,0.13,2.43,1.54,6.35,0.17,0.41)}
    for AA in seq:
        if AA not in library:  
            raise ValueError(f"Invalid amino acid letter: {AA}") 
        properties = library[AA]
        modified_line = properties
        modified_line_str = ','.join(str(item) for item in modified_line) 
        modified_lines_3.append(modified_line_str)
    return modified_lines_3
#构建PSSM数据库和矩阵
def scoring_matrix_1(inputfile2): #先构建了矩阵pssm搜索数据库
    import os
    target_folder = 'ScoringMatrix'
    with open(inputfile2,'r',encoding='utf-8') as file:
        lines = file.readlines()
    #将输入文件中基因名和基因序列提取出来构建搜索库，并用于后续进行pssm搜索
    modified_lines = []
    for line in lines:
        columns = line.strip().split(',')
        if len(columns) > 0 and 'OG' not in columns:
            modified_line = '>' + columns[0] + '\n' + columns[1] + '\n'
            modified_lines.append(modified_line)
        else:
            modified_lines.append(line)
    with open('output.fasta', 'w', encoding='utf-8') as file:
        file.writelines(modified_lines)
    os.system('makeblastdb -dbtype prot -in output.fasta -input_type fasta -out test.blastdb') #使用全部的搜索到的阳性/阴性池中的数据构建一个搜索库
    with open('output.fasta', 'r', encoding='utf-8') as file:
        lines = file.readlines()
    #对每一个基因和序列搜索的结果进行处理和分析，有多少基因，就会形成多少个子文件
def scoring_matrix_2(seq,modified_lines_4,file_index,columns):
    import os
    import pandas as pd
    target_folder = 'ScoringMatrix'
    output_file = []
    output_filename = f'output_{file_index:02d}.txt'
    targetfile2 = os.path.join(target_folder,output_filename)
    with open(targetfile2, 'w', encoding='utf-8') as outfile:
        modified_line = '>' + columns[0] + '\n' + columns[1] + '\n'
        outfile.writelines(modified_line)
    var = str(output_filename)
    os.environ['var'] = var
    vas = str(file_index)
    os.environ['vas'] = vas
    os.system('nohup psiblast -db test.blastdb -query ScoringMatrix/$var -num_iterations 10 -out_ascii_pssm ScoringMatrix/$vas.pssm >> nohup_res2.txt')
    pssm_file = str(file_index) + '.pssm'
    targetfile3 = os.path.join(target_folder,pssm_file)
    pssm1 = '0A.pssm' #这是一个全是0的文件，跑程序之前自己设置！！！！！！！！！！！！！！！！！
    
    aa_num_1 = len(seq) + 2
    if os.path.exists(targetfile3):
        #将每一个pssm搜索结果进行处理，形成列表格式
        df = pd.read_csv(targetfile3, header = None )
    else: 
        df = pd.read_csv(pssm1, header = None )
    file1 = f'result1_{file_index:02}.txt'
    select = df.iloc[2:aa_num_1,0]   #pssm搜索结果中根据含有氨基酸的数量，提取行数
    select.to_csv(file1,sep=',',header = None)
    df = pd.read_csv(file1, header = None, sep='\s+')
    select1 = df.iloc[:,3:23]  #pssm搜索结果中提取的列数固定
    file2 = f'result2_{file_index:02}.txt'
    select1.to_csv(file2,index=False,header = None)
    with open(file2,'r') as re2:
        lines_5 = re2.readlines()
        for line_5 in lines_5:
            modified_line = line_5.rstrip('\n') #由于使用pandas，已经将数据直接str化，可以直接省略下一步。
            #modified_line_str = ','.join(item for item in modified_line) 
            modified_lines_4.append(modified_line)
    return modified_lines_4
#计算二级结构
def second_structure_annotation(modified_lines_5,file_index,columns):                 #使用python中的peptidebuilder函数
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import os
    L_seq = columns[1]
    L_mol = Chem.rdmolfiles.MolFromFASTA(L_seq,flavor=0)
    L_mol = Chem.AddHs(L_mol)
    #AllChem.EmbedMolToPDBBlock(L_mol)    生成三级结构
    L_pdbblock = Chem.MolToPDBBlock(L_mol)
    output_filename = f'output_{file_index:02d}.pdb'  
    target_folder = 'SecondStructureAnnotation'
    target_file = os.path.join(target_folder,output_filename)
    with open(target_file,'w',encoding='utf-8') as f:
        f.write(L_pdbblock)
    vas = str(output_filename)
    os.environ['vas'] = vas
    os.system('mkdssp -i SecondStructureAnnotation/$vas -o SecondStructureAnnotation/$vas.dssp')
    output_filename2 = f'output_{file_index:02d}.pdb.dssp'
    target_file2 = os.path.join(target_folder,output_filename2)
    with open (target_file2,'r') as file2:
        lines_1 = file2.readlines()[28:]
    for line_1 in lines_1:
        modified_line = str(line_1)
        modified_line_str = ','.join(str(item) for item in modified_line)
        select = modified_line_str[32]  #提取二级结构缩写的字母，出现在列表中的第32个str
        #print (select)
        library = { 'B':(1,0,0,0,0,0,0), 
                    'E':(0,1,0,0,0,0,0), 
                    'G':(0,0,1,0,0,0,0), 
                    'H':(0,0,0,1,0,0,0),
                    'I':(0,0,0,0,1,0,0),
                    'T':(0,0,0,0,0,1,0),
                    'S':(0,0,0,0,0,0,1) }
        if select not in library:
            structure = (0,0,0,0,0,0,0)
        else:
            structure = library[select]
        structure_str = ','.join(str(item) for item in structure)
        modified_lines_5.append(structure_str)
    return modified_lines_5
#推测RSA
def Putative_RSA(modified_lines_6,columns,file_index):  
    import os
    import pandas as pd
    target_folder = 'PutativeRSA'
    output_filename = f'output_{file_index:02d}.txt'
    targetfile1 = os.path.join(target_folder,output_filename)
    modified_line = '>' + columns[0] + '\n' + columns[1] + '\n'
    with open(targetfile1, 'w', encoding='utf-8') as file:
        file.writelines(modified_line)
    var = str(output_filename)
    os.environ['var'] = var
    #os.system('cd PutativeRSA')
    os.system('nohup ASAquick PutativeRSA/$var >> nohup_res1.txt')
    #os.system('cd -')
    folder = f'asaq.output_{file_index:02d}.txt'
    filename = 'rasaq.pred'
    targetfile2 = os.path.join(folder,filename)
    df = pd.read_csv(targetfile2,header = None, sep = '\s+')
    df = df.iloc[:,2]
    file = f'ASA_{file_index:02d}.txt'
    df.to_csv(file,header = None, index = None)
    with open(file,'r') as file1:
        lines_1 = file1.readlines()[:]
    for line_1 in lines_1:
        modified_lines_6.append(line_1)
    return modified_lines_6


#read_folder和deal_folder两个文件用来对神经网络模型训练原始文件夹的读取和处理
def read_folder(folder_path, label):  
    import os
    import numpy as np
    data = []     
    for filename in os.listdir(folder_path):  
        if filename.endswith('.txt'):  
            file_path = os.path.join(folder_path, filename)  
            with open(file_path, 'r') as file:  
                # 假设每行是一个逗号分隔的数值列表  
                lines = file.readlines()
                rows = []
                for line in lines:
                    line = line.strip().split(',')
                    colmuns1 = line[:19]
                    colmuns2 = line[19:]
                    colmuns = np.concatenate((colmuns1, colmuns2), axis=0)
                    print(len(colmuns))
                    rows.append(colmuns)
                #rows = [list(map(float, row.strip().split(','))) for row in file.readlines()]  
                # 将二维列表转换为numpy数组，并重塑为(1, 10, 40, 1)  
                # 注意：这里我们假设每个文件都是一个样本，你需要根据你的实际情况调整 
                #row = rows[:8] 
                #print (file_path)
                print(len(rows))
                if len(rows) > 7 :
                    #sample1 = np.array(rows[:19])
                    sample = np.array(rows).reshape(1, 8, 58, 1) 
                    #sample = np.concatenate((sample1, sample2), axis=1).reshape(1, 8, 38, 1) 
                    data.append(sample)  
    # 将所有样本合并为一个大的numpy数组  
    data = np.concatenate(data, axis=0)  
    return data, np.full(data.shape[0], label)   
# 遍历文件夹中的所有文件  
def deal_folder(folder_path):
    import os
    for filename in os.listdir(folder_path):  
        if filename.endswith('.txt'):  # 确保只处理txt文件  
            # 构建原文件和新文件的完整路径  
            file_path = os.path.join(folder_path, filename)  
            # 读取原文件的前8行  
            with open(file_path, 'r', encoding='utf-8') as file:  
                lines = [line.rstrip() for line in file.readlines()[-8:]]  # 读取前8行并去除行尾的换行符    
                # 如果不足8行，lines将包含实际读取的行数       
                # 将前8行（或实际存在的行数）写入新文件  
            with open(file_path, 'w', encoding='utf-8') as file:  
                file.writelines(line + '\n' for line in lines)  # 添加最后一个换行符（如果需要）  

