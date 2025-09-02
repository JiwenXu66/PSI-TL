#该函数文件获得的参数用于随机森林模型分析
###第1行到第218行的24个函数是用来分析氨基酸序列的特性，返回的结果是用数字描述氨基酸序列
def Hydrophobicity_1(seq):
    class1 = ('R','K','E','D','Q','N'); class2 = ('G','A','S','T','P','H','Y'); class3 = ('C','L','V','I','M','F','W')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Normalized_van_der_Waals_volume_2(seq):
    class1 = ('G','A','S','T','P','D'); class2 = ('N','V','E','Q','I','L'); class3 = ('M','H','K','F','R','Y','W')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq   
def Polarity_3(seq):
    class1 = ('L','I','F','W','C','M','V','Y'); class2 = ('P','A','T','G','S'); class3 = ('H','Q','R','K','N','E','D')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Polarizability_4(seq):
    class1 = ('G','A','S','D','T'); class2 = ('C','P','N','V','E','Q','I','L'); class3 = ('K','M','H','F','R','Y','W')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Charge_5(seq):
    class1 = ('K','R'); class2 = ('A','N','C','Q','G','H','I','L','M','F','P','S','T','W','Y','V'); class3 = ('D','E')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Secondary_structure_6(seq):
    class1 = ('E','A','M','L','Q','K','R','H'); class2 = ('V','I','Y','C','W','F','T'); class3 = ('G','N','P','S','D')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Solvent_accessibility_7(seq):
    class1 = ('A','L','F','C','G','I','V','W'); class2 = ('P','K','Q','E','N','D'); class3 = ('M','P','S','T','H','Y')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Surface_tension_8(seq):
    class1 = ('G','Q','D','N','A','H','R'); class2 = ('K','T','S','E','C'); class3 = ('I','L','M','F','P','W','Y','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_protein_interface_hotspot_propensity_Bogan_9(seq):
    class1 = ('D','H','I','K','N','P','R','W','Y'); class2 = ('E','Q','S','T','G','A','M','F'); class3 = ('C','L','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_protein_interface_propensity_Ma_10(seq):
    class1 = ('C','D','F','M','P','Q','R','W','Y'); class2 = ('A','G','H','V','L','N','S','T'); class3 = ('E','I','K')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_DNA_interface_propensity_Schneide_11(seq):
    class1 = ('G','K','N','Q','R','S','T','Y'); class2 = ('A','D','E','F','H','I','L','V','W'); class3 = ('C','M','P')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_DNA_interface_propensity_Ahmad_12(seq):
    class1 = ('G','H','K','N','Q','R','S','T','Y'); class2 = ('A','D','E','F','I','P','V','W'); class3 = ('C','M','L')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_RNA_interface_propensity_Kim_13(seq):
    class1 = ('H','K','M','R','Y'); class2 = ('F','G','I','L','N','P','Q','S','V','W'); class3 = ('C','D','E','A','T')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_RNA_interface_propensity_Ellis_14(seq):
    class1 = ('H','G','K','M','R','S','Y','W'); class2 = ('A','F','I','N','P','Q','T'); class3 = ('C','D','E','L','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_RNA_interface_propensity_Phipps_15(seq):
    class1 = ('H','K','M','R','S','Q'); class2 = ('A','D','E','F','G','L','N','P','V','Y'); class3 = ('C','I','T','W')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_ligand_binding_site_propensity_Khazanov_16(seq):
    class1 = ('C','F','H','W','Y'); class2 = ('G','I','L','N','M','S','T','R'); class3 = ('A','E','D','K','P','Q','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Protein_ligand_valid_binding_site_propensity_Khazanov_17(seq):
    class1 = ('C','F','H','W','Y','M'); class2 = ('D','G','I','L','N','S','T','V'); class3 = ('A','E','K','P','Q','R')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Propensity_for_protein_ligand_polar_aromatic_non_bonded_interactions_Imai_18(seq):
    class1 = ('D','E','H','R','Y'); class2 = ('C','F','K','M','N','Q','S','T','W'); class3 = ('A','G','I','L','P','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Molecular_Weight_19(seq):
    class1 = ('A','G','S'); class2 = ('C','D','E','H','I','K','L','M','N','Q','P','T','V'); class3 = ('F','R','W','Y')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def cLogP_20(seq):
    class1 = ('R','K','D','N','E','Q','H'); class2 = ('P','Y','S','T','G','A','C','V'); class3 = ('W','M','F','L','I')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def No_of_hydrogen_bond_donor_in_side_chain_21(seq):
    class1 = ('H','K','N','Q','R'); class2 = ('D','E','S','T','W','Y'); class3 = ('A','C','G','F','I','L','M','P','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def No_of_hydrogen_bond_acceptor_in_side_chain_22(seq):
    class1 = ('D','E','H','N','Q','R'); class2 = ('K','S','T','W','Y'); class3 = ('A','C','G','F','I','L','M','P','V')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Solubility_in_water_23(seq):
    class1 = ('A','C','G','K','R','T'); class2 = ('E','F','H','I','L','M','N','P','Q','S','V','W'); class3 = ('D','Y')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
def Amino_acid_flexibility_index_24(seq):
    class1 = ('E','G','K','N','Q','S'); class2 = ('A','D','H','I','P','R','T','V'); class3 = ('C','F','L','M','W','Y')
    ca_map = { char: 1 for char in class1 } #设置映射字典，并完成class1的映射
    for char in class2:  
        ca_map[char] = 2  
    for char in class3:  
        ca_map[char] = 3  
    transformed_seq = ''.join(str(ca_map.get(char, '?')) for char in seq)
    return transformed_seq
###composition（返回3个数值）,transition（返回3个数值）,distribution（返回6个数值）三个函数是用来计算组成成分的，和上面的24个函数结合起来使用。
def composition(transformed_seq):
    num = len(transformed_seq)
    count_C1 = transformed_seq.count('1'); count_C2 = transformed_seq.count('2'); count_C3 = transformed_seq.count('3')
    ratio_C1 = count_C1 / num; ratio_C2 = count_C2 / num; ratio_C3 = count_C3 / num
    ratio = f"{ratio_C1:.2f},{ratio_C2:.2f},{ratio_C3:.2f}" 
    return ratio
def transition(transformed_seq):
    num = len(transformed_seq)
    count_T12 = transformed_seq.count('12'); count_T13 = transformed_seq.count('13'); count_T23 = transformed_seq.count('23')
    count_T21 = transformed_seq.count('21'); count_T31 = transformed_seq.count('31'); count_T32 = transformed_seq.count('32')
    ratio_T12 = (count_T12 + count_T21) / num; ratio_T13 = (count_T13 + count_T31) / num; ratio_T23 = (count_T23 + count_T32) / num
    ratio = f"{ratio_T12:.2f},{ratio_T13:.2f},{ratio_T23:.2f}"
    return ratio
def distribution(transformed_seq): 
    #由于序列氨基酸数量较少，参考文献中分别识别0，25%，50%，75%，100%的位置需要1，2，3三种类型氨基酸每种至少5个。所以我们修改为只识别第一个和最后一个
    num = len(transformed_seq)
    position_1 = [idx for idx, char in enumerate(transformed_seq) if char == '1']
    position_2 = [idx for idx, char in enumerate(transformed_seq) if char == '2']
    position_3 = [idx for idx, char in enumerate(transformed_seq) if char == '3']
    #print (transformed_seq)
    #print (transformed_seq,len(position_1),len(position_2),len(position_3))
    if len(position_1) >=2:
        diatribution_1_start = position_1[0] / num
        diatribution_1_end = position_1[-1] / num 
    else: 
        diatribution_1_start = 0
        diatribution_1_end = 0
    if len(position_2) >=2:
        diatribution_2_start = position_2[0] / num
        diatribution_2_end = position_2[-1] / num
    else: 
        diatribution_2_start = 0
        diatribution_2_end = 0
    if len(position_3) >=2:
        diatribution_3_start = position_3[0] / num
        diatribution_3_end = position_3[-1] / num 
    else: 
        diatribution_3_start = 0
        diatribution_3_end = 0
    ratio = f"{diatribution_1_start:.2f},{diatribution_2_start:.2f},{diatribution_3_start:.2f},{diatribution_1_end:.2f},{diatribution_2_end:.2f},{diatribution_3_end:.2f}"
    return ratio
###下面两个个函数主要用来计算氨基酸比例和二肽比例，一共抽提出420个参数
def AA_ratio_caculate(seq,amino): #用于计算氨基酸的比例
    seq_num = len(seq)
    AA_count = seq.count(amino)
    AA_ratio = AA_count / seq_num
    return AA_ratio
def dipeptide_ratio_caculate(seq,i):   #用于计算二肽的比例
    dipeptide_num = len(seq) - 1
    dipeptide_count = seq.count(i)
    dipeptide_ratio = dipeptide_count / dipeptide_num
    return dipeptide_ratio
###下面三个函数用抽提对氨基酸序列描述的参数，AA_composition（20个参数），dipeptide_composition（400个参数），Character_select（288个参数+label）
def AA_composition(inputfile1,modified_lines_1): #计算序列中氨基酸的比例，共计20个
    with open(inputfile1,'r',encoding='utf-8') as file:
        lines = file.readlines()
    for line in lines:
        columns = line.strip().split(',')
        if len(columns) >=1:
            ratio_all = []
            seq = str(columns[1])
            aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
            for amino in aa:
                AA_ratio = AA_ratio_caculate(seq,amino)
                ratio_all.append(AA_ratio)
                ratio_all_str = ','.join(str(item) for item in ratio_all)  #将含有float的list形式转化为str格式，耗费较长时间和算力
            modified_line = ratio_all_str
            modified_lines_1.append(modified_line)
    return modified_lines_1
def dipeptide_composition(inputfile1,modified_lines_2):  #计算序列中两个氨基酸组成二肽的比例，共计400个组合
    from itertools import product
    aa = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
    dipeptides = [''.join(pair) for pair in product(aa, repeat=2)] #按顺序输出400个二肽        
    with open(inputfile1,'r',encoding='utf-8') as file:
        lines = file.readlines()
    for line in lines:
        columns = line.strip().split(',')
        if len(columns) >=1:
            ratio_all = []
            seq = str(columns[1])
            for i in dipeptides:
                dipeptide_ratio = dipeptide_ratio_caculate(seq,i)
                ratio_all.append(dipeptide_ratio)
                ratio_all_str = ','.join(str(item) for item in ratio_all)  #将含有float的list形式转化为str格式，耗费较长时间和算力
            modified_line = ratio_all_str
            modified_lines_2.append(modified_line)
    return modified_lines_2
def Character_select(inputfile1,modified_lines_3,lable):
    
    #提取氨基酸序列信息——1：阳性池数据信息提取
    with open(inputfile1,'r',encoding='utf-8') as file:
        lines = file.readlines()
        for line in lines:
            columns = line.strip().split(',')
            if len(columns) >=1:
                ratio_all = []
                seq = str(columns[1])
                #print (seq)
                transformed_seq = Hydrophobicity_1(seq); composition1 = composition(transformed_seq); transition1 = transition(transformed_seq); distribution1 = distribution(transformed_seq)
                transformed_seq = Normalized_van_der_Waals_volume_2(seq); composition2 = composition(transformed_seq); transition2 = transition(transformed_seq); distribution2 = distribution(transformed_seq)
                transformed_seq = Polarity_3(seq); composition3 = composition(transformed_seq); transition3 = transition(transformed_seq); distribution3 = distribution(transformed_seq)
                transformed_seq = Polarizability_4(seq); composition4 = composition(transformed_seq); transition4 = transition(transformed_seq); distribution4 = distribution(transformed_seq)
                transformed_seq = Charge_5(seq); composition5 = composition(transformed_seq); transition5 = transition(transformed_seq); distribution5 = distribution(transformed_seq)
                transformed_seq = Secondary_structure_6(seq); composition6 = composition(transformed_seq); transition6 = transition(transformed_seq); distribution6 = distribution(transformed_seq)
                transformed_seq = Solvent_accessibility_7(seq); composition7 = composition(transformed_seq); transition7 = transition(transformed_seq); distribution7 = distribution(transformed_seq)
                transformed_seq = Surface_tension_8(seq); composition8 = composition(transformed_seq); transition8 = transition(transformed_seq); distribution8 = distribution(transformed_seq)
                transformed_seq = Protein_protein_interface_hotspot_propensity_Bogan_9(seq); composition9 = composition(transformed_seq); transition9 = transition(transformed_seq); distribution9 = distribution(transformed_seq)
                transformed_seq = Protein_protein_interface_propensity_Ma_10(seq); composition10 = composition(transformed_seq); transition10 = transition(transformed_seq); distribution10 = distribution(transformed_seq)
                transformed_seq = Protein_DNA_interface_propensity_Schneide_11(seq); composition11 = composition(transformed_seq); transition11 = transition(transformed_seq); distribution11 = distribution(transformed_seq)
                transformed_seq = Protein_DNA_interface_propensity_Ahmad_12(seq); composition12 = composition(transformed_seq); transition12 = transition(transformed_seq); distribution12 = distribution(transformed_seq)
                transformed_seq = Protein_RNA_interface_propensity_Kim_13(seq); composition13 = composition(transformed_seq); transition13 = transition(transformed_seq); distribution13 = distribution(transformed_seq)
                transformed_seq = Protein_RNA_interface_propensity_Ellis_14(seq); composition14 = composition(transformed_seq); transition14 = transition(transformed_seq); distribution14 = distribution(transformed_seq)
                transformed_seq = Protein_RNA_interface_propensity_Phipps_15(seq); composition15 = composition(transformed_seq); transition15 = transition(transformed_seq); distribution15 = distribution(transformed_seq)
                transformed_seq = Protein_ligand_binding_site_propensity_Khazanov_16(seq); composition16 = composition(transformed_seq); transition16 = transition(transformed_seq); distribution16 = distribution(transformed_seq)
                transformed_seq = Protein_ligand_valid_binding_site_propensity_Khazanov_17(seq); composition17 = composition(transformed_seq); transition17 = transition(transformed_seq); distribution17 = distribution(transformed_seq)
                transformed_seq = Propensity_for_protein_ligand_polar_aromatic_non_bonded_interactions_Imai_18(seq); composition18 = composition(transformed_seq); transition18 = transition(transformed_seq); distribution18 = distribution(transformed_seq)
                transformed_seq = Molecular_Weight_19(seq); composition19 = composition(transformed_seq); transition19 = transition(transformed_seq); distribution19 = distribution(transformed_seq)
                transformed_seq = cLogP_20(seq); composition20 = composition(transformed_seq); transition20 = transition(transformed_seq); distribution20 = distribution(transformed_seq)
                transformed_seq = No_of_hydrogen_bond_donor_in_side_chain_21(seq); composition21 = composition(transformed_seq); transition21 = transition(transformed_seq); distribution21 = distribution(transformed_seq)
                transformed_seq = No_of_hydrogen_bond_acceptor_in_side_chain_22(seq); composition22 = composition(transformed_seq); transition22 = transition(transformed_seq); distribution22 = distribution(transformed_seq)
                transformed_seq = Solubility_in_water_23(seq); composition23 = composition(transformed_seq); transition23 = transition(transformed_seq); distribution23 = distribution(transformed_seq)
                transformed_seq = Amino_acid_flexibility_index_24(seq); composition24 = composition(transformed_seq); transition24 = transition(transformed_seq); distribution24 = distribution(transformed_seq)
                modified_line = [composition1 +','+ composition2 +','+ composition3 +','+ composition4 +','+ composition5 +','+ composition6 +','+ composition7 +','+ composition8 +','+ composition9 +','+ composition10 +','+ 
                                composition11 +','+ composition12 +','+ composition13 +','+ composition14 +','+ composition15 +','+ composition16 +','+ composition17 +','+ composition18 +','+ composition19 +','+ 
                                composition20 +','+ composition21 +','+ composition22 +','+ composition23 +','+ composition24 +','+ transition1 +','+ transition2 +','+ transition3 +','+ transition4 +','+ 
                                transition5 +','+ transition6 +','+ transition7 +','+ transition8 +','+ transition9 +','+ transition10 +','+ transition11 +','+ transition12 +','+ transition13 +','+ transition14
                                +','+ transition15 +','+ transition16 +','+ transition17 +','+ transition18 +','+ transition19 +','+ transition20 +','+ transition21 +','+ transition22 +','+ transition23 
                                +','+ transition24 +','+ distribution1 +','+ distribution2 +','+ distribution3 +','+ distribution4 +','+ distribution5 +','+ distribution6 +','+ distribution7 +','+ distribution8
                                +','+ distribution9 +','+ distribution10 +','+ distribution11 +','+ distribution12 +','+ distribution13 +','+ distribution14 +','+ distribution15 +','+ distribution16 +','+ distribution17
                                +','+ distribution18 +','+ distribution19 +','+ distribution20 +','+ distribution21 +','+ distribution22 +','+ distribution23 +','+ distribution24 +','+ str(lable)]
                modified_line_str = ','.join(str(item) for item in modified_line)
                modified_lines_3.append(modified_line_str)
    return modified_lines_3

