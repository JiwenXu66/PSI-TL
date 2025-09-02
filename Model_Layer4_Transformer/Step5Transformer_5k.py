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
from sklearn.ensemble import RandomForestClassifier #用于分析分类问题
from sklearn.model_selection import train_test_split  
from sklearn.preprocessing import LabelEncoder  
from sklearn.metrics import accuracy_score, roc_curve, f1_score, recall_score, precision_score, auc   # 用于分类问题评估  
#from rdkit import Chem
#from rdkit.Chem import AllChem
import tensorflow as tf   
from tensorflow.keras.layers import Layer, Embedding, LayerNormalization, MultiHeadAttention, Dense, GlobalAveragePooling1D
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dropout  
from tensorflow.keras.initializers import TruncatedNormal  #初始化器，避免梯度消失或爆炸的问题
from tensorflow.keras.models import Model, Sequential, load_model   
import FunctionOne
import FunctionTwo
import FunctionZero
from tensorflow.keras.utils import register_keras_serializable 
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold

########Transformer框架的应用
class PositionalEncoding(Layer):   #为输入序列添加位置编码，定义了一个层layer，可以和keras中其他层一样使用的类，这个类是Layer的子类
    def __init__(self, position_size, d_model, **kwargs):  #初始化方法，三个参数：位置编码长度，模型的特征维度，任意数量的关键字参数
        super(PositionalEncoding, self).__init__(**kwargs)  #调用了父类Layer中__init__函数，并捕获关键字参数**kwargs
        self.position_size = position_size  #将传入的参数，保存到类的实例变量里，后续可以直接调用实例变量
        self.d_model = d_model  
        self.pos_encoding = self.add_weight(  
            name = 'pos_encoding',  
            shape=(self.position_size, self.d_model),  
            initializer=self._positional_encoding_init,  
            trainable=False,  
        )  #add_weight属于Layer的基类，矩阵性状shape由两个变脸确定，权重被设定为不可以训练
    def _positional_encoding_init(self, shape, dtype=None):  #初始孵化器，接收性状shape和stype作为参数
        encoding = tf.range(start=0, limit=shape[0], delta=1.0)  
        encoding = tf.expand_dims(encoding, 1)  
        encoding = tf.math.sin(encoding / (10000 ** (tf.range(0, shape[1], 2.0) / tf.cast(shape[1], tf.float32))))  
        encoding = tf.concat([encoding, tf.math.cos(encoding)], axis=1)  #！这里有一个错误：应该首先计算余弦编码，然后将其与正弦编码拼接，因为encoding变量在拼接前已被正弦编码覆盖。
        return tf.cast(encoding, dtype=tf.float32)  #确保位置编码矩阵的数据类型是tf.float32
    def call(self, inputs):  #这是PositionalEncoding层的call方法，它定义了层的前向传播逻辑，
        return inputs + self.pos_encoding#[:tf.shape(inputs)[0], :]  #每个输入位置都会根据其位置获得一个独特的位置编码
  
class TransformerEncoderBlock(Layer):  #定义了transformer框架
    def __init__(self, d_model, num_heads, dff, rate=0.1, **kwargs):  #dff 前反馈网络中的维度
        super(TransformerEncoderBlock, self).__init__(**kwargs)  # 可以理解为这一步在调用tensorflow中Layer函数，并在这个基础上进行增加
        self.mha = MultiHeadAttention(num_heads=num_heads, key_dim=d_model)  
        self.ffn = tf.keras.Sequential(  
            [  
                Dense(dff, activation='relu'),  
                Dense(d_model)  
            ]  
        )  #transformer框架中前反馈网络由两个全连接层组成，第一个连接层进行线性变换，使用relu函数，第二个全连接层是将relu函数处理后的特征映射回原始维度
        self.layernorm1 = LayerNormalization(epsilon=1e-6)  #归一化层
        self.layernorm2 = LayerNormalization(epsilon=1e-6)  
        self.dropout1 = tf.keras.layers.Dropout(rate)  #丢弃神经元
        self.dropout2 = tf.keras.layers.Dropout(rate)  
    def call(self, inputs, training, mask):  
        attn_output = self.mha(inputs, inputs, attention_mask=mask)  
        attn_output = self.dropout1(attn_output, training=training)  
        out1 = self.layernorm1(inputs + attn_output)  
        ffn_output = self.ffn(out1)  
        ffn_output = self.dropout2(ffn_output, training=training)  
        return self.layernorm2(out1 + ffn_output)  
  
class TransformerEncoder(Layer):  
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size, maximum_position_encoding, rate=0.1, **kwargs):  
        super(TransformerEncoder, self).__init__(**kwargs)  
        self.num_heads = num_heads  
        self.dff = dff  
        self.maximum_position_encoding = maximum_position_encoding
        self.d_model = d_model  
        self.num_layers = num_layers  

        self.embedding = Embedding(input_vocab_size, d_model)  
        self.pos_encoding = PositionalEncoding(maximum_position_encoding, d_model)  
        self.enc_layers = [TransformerEncoderBlock(d_model, num_heads, dff, rate) for _ in range(num_layers)]  
        self.dropout = tf.keras.layers.Dropout(rate)  
        
        #self.pool = GlobalAveragePooling1D()  # 添加全局平均池化层  
        #self.output_layer = Dense(1, activation='sigmoid')  # 添加输出层 
    def call(self, inputs, training, mask):  
        seq_len = tf.shape(inputs)[1]  
        # 添加位置编码  
        inputs = self.pos_encoding(self.embedding(inputs))  
        # 添加dropout  
        inputs = self.dropout(inputs, training=training)  
        # 编码器层堆叠  
        for i in range(self.num_layers):  
            inputs = self.enc_layers[i](inputs, training=training, mask=mask)  
        #pooled_output = self.pool(inputs)   # 全局平均池化
        #inputs = self.output_layer(inputs)
        return inputs  
    def get_config(self):  
        config = super(TransformerEncoder, self).get_config()  
        config.update({  
            'num_layers': self.num_layers,  
            'd_model': self.d_model,  
            'num_heads': self.num_heads,
            'dff': self.dff, 
            'input_vocab_size': self.embedding.input_dim,  # Embedding 层的 input_dim 属性  
            'maximum_position_encoding': self.maximum_position_encoding,
            'rate': self.dropout.rate,  # Dropout 层的 rate 属性  
        })  
        return config  

class TransformerEncoderForBinaryClassification(Model):  
    def __init__(self, transformer_encoder, num_classes=2, **kwargs):  
        super(TransformerEncoderForBinaryClassification, self).__init__(**kwargs)  
        self.transformer_encoder = transformer_encoder #(inputs, training=training, mask=mask)
        self.classifier = Dense(1, activation='sigmoid')  # 对于二分类，通常是sigmoid，但softmax也可以（取第一个输出）  
  
    def call(self, inputs, training=False, mask=None):  
        encoder_outputs = self.transformer_encoder(inputs, training=training, mask=mask)  
        print(encoder_outputs.shape)
        cls_token = encoder_outputs[:, 0, :]  # 假设CLS token在第一个位置  
        classifier_output = self.classifier(cls_token)
        print(classifier_output.shape)
        return classifier_output   
    def get_build_config(self):  
        config = super(TransformerEncoderForBinaryClassification, self).get_config()  
        config.update({  
            'transformer_encoder_config': self.transformer_encoder.get_config() if hasattr(self.transformer_encoder, 'get_config') else None,  
            'num_classes': self.classifier.units  # 或者直接使用 1，因为对于二分类它是固定的  
        })  
        return config  
 
    @classmethod  
    def build_from_config(cls, config):  
        transformer_encoder = TransformerEncoder.from_config(config['transformer_encoder_config'])    
        return cls(transformer_encoder=transformer_encoder, num_classes=config['num_classes'])

def read_files(filename,label):  
    data = []  
    labels = [] 
    AA_index = {'A': 0,  # Alanine  
                'R': 1,  # Arginine  
                'N': 2,  # Asparagine  
                'D': 3,  # Aspartic acid  
                'C': 4,  # Cysteine  
                'Q': 5,  # Glutamine  
                'E': 6,  # Glutamic acid  
                'G': 7,  # Glycine  
                'H': 8,  # Histidine  
                'I': 9,  # Isoleucine  
                'L': 10, # Leucine  
                'K': 11, # Lysine  
                'M': 12, # Methionine  
                'F': 13, # Phenylalanine  
                'P': 14, # Proline  
                'S': 15, # Serine  
                'T': 16, # Threonine  
                'W': 17, # Tryptophan  
                'Y': 18, # Tyrosine  
                'V': 19  # Valine   
                }
    with open(filename,'r') as f:
        lines = f.readlines()
        #lines = set(lines)
    for line in lines:
        line = line.strip().split(',')
        if len(line) >0 and 'p-value' not in line:
            column = line[1]
            column = column[-8:]
            if all(aa in AA_index for aa in column):
                AA_indices = [ AA_index[aa] for aa in column]
                tensor = tf.constant(AA_indices,dtype=tf.int32)
                data.append(tensor)    
                labels.append(int(label))    
    return data,labels



# 读取数据  
positive_data, positive_labels = read_files('PositiveResult2_shaixuan.txt', '1')  
negative_data, negative_labels = read_files('NegativeResult2_shaixuan.txt', '0')  
# 合并数据  
inputs = tf.concat((positive_data[:], negative_data[:]), axis=0)  #这时候有13000多个数据
labels = tf.concat((positive_labels[:], negative_labels[:]), axis=0)
#print(inputs.shape)

# 5倍交叉验证
kf = KFold(n_splits=5, shuffle=True, random_state=42)
results = []
 
for fold, (train_index, val_index) in enumerate(kf.split(inputs)):
    print(f"Training for fold {fold + 1}")
    # 将 NumPy 索引数组转换为 TensorFlow 张量
    train_index_tensor = tf.convert_to_tensor(train_index, dtype=tf.int32)
    val_index_tensor = tf.convert_to_tensor(val_index, dtype=tf.int32)
    # 使用 tf.gather 根据索引收集数据
    train_inputs = tf.gather(inputs, train_index_tensor)
    val_inputs = tf.gather(inputs, val_index_tensor)
    train_labels = tf.gather(labels, train_index_tensor)
    val_labels = tf.gather(labels, val_index_tensor)
    #参数设计
    num_layers = 2  #编码器堆叠层数
    d_model = 128    #数据扩增后的维度,可以自己设定
    num_heads = 6    #多头自注意力的头数
    dff = 1024       #前反馈神经网络中的维数
    input_vocab_size = 20  #输入词汇的大小
    maximum_position_encoding = 8  #输入词汇的最大位置编码啊
    # 创建transformer_encoder框架  
    transformer_encoder = TransformerEncoder(num_layers, d_model, num_heads, dff, input_vocab_size, maximum_position_encoding) 
    # 在transformer_encoder的基础上构建二分类模型  
    model_Transformer = TransformerEncoderForBinaryClassification(transformer_encoder)
    # 编译二分类模型  
    model_Transformer.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])   # 因为是二分类问题 
    # 训练模型
    mask = None
    model_Transformer.fit(train_inputs, train_labels, epochs=10, batch_size=64, validation_split=0.2)  #batch_size adam分析中，每次迭代的样本量，32的倍数
    
    # 评估模型
    val_loss, val_accuracy = model_Transformer.evaluate(val_inputs, val_labels, verbose=0)
    y_probs = model_Transformer.predict(val_inputs).ravel()
    y_pred = (y_probs > 0.5).astype("int32")
    
    f1 = f1_score(val_labels, y_pred)
    recall = recall_score(val_labels, y_pred)
    precision = precision_score(val_labels, y_pred)
    accuracy = accuracy_score(val_labels, y_pred)
    fpr, tpr, thresholds = roc_curve(val_labels, y_probs)
    roc_auc = auc(fpr, tpr)
    
    results.append({
        'fold': fold + 1,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'auc': roc_auc
    })
    
    print(f"Fold {fold + 1} - Accuracy: {accuracy}, Precision: {precision}, Recall: {recall}, F1: {f1}, AUC: {roc_auc}")
 
# 打印所有结果
for result in results:
    print(f"Fold {result['fold']} - Accuracy: {result['accuracy']}, Precision: {result['precision']}, Recall: {result['recall']}, F1: {result['f1']}, AUC: {result['auc']}")







