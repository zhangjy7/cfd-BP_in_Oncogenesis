import pandas as pd
import numpy as np
import sys
import os
import re  
import json
import random
import time
import collections
import math
from scipy.stats import poisson
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from sklearn import preprocessing

print('starting')
ii = sys.argv[1]
ii = int(ii)-1
print(ii)


Human_Clust_path = 'Human_ClustGO/'
GO_path = 'GO_annotation/'
path_TPM = 'TPM/'
T_path = 'TPM/TCGA_TPM2/'
G_path = 'TPM/GTEx_TPM2/'

TCGA_project = ['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']
GTEx_tissues = ['Adipose-Subcutaneous', 'Adipose-Visceral', 'Adrenal_Gland', 'Artery-Aorta','Artery-Coronary', 'Artery-Tibial', 'Bladder', 'Brain-Amygdala','Brain-Anterior_cingulate_cortex', 'Brain-Caudate', 'Brain-Cerebellar_Hemisphere', 'Brain-Cerebellum','Brain-Cortex','Brain-Frontal_Cortex', 'Brain-Hippocampus', 'Brain-Hypothalamus', 'Brain-Nucleus_accumbens', 'Brain-Putamen','Brain-Spinal_cord', 'Brain-Substantia_nigra', 'Breast-Mammary_Tissue', 'Cells-Cultured_fibroblasts', 'Cells-EBV-transformed_lymphocytes', 'Cervix-Ectocervix', 'Cervix-Endocervix', 'Colon-Sigmoid', 'Colon-Transverse','Esophagus-Gastroesophageal_Junction', 'Esophagus-Mucosa','Esophagus-Muscularis', 'Fallopian Tube', 'Heart-Atrial_Appendage', 'Heart-Left_Ventricle', 'Kidney-Cortex', 'Kidney-Medulla', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle-Skeletal', 'Nerve-Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin-Not_Sun_Exposed_Suprapubic', 'Skin-Sun_Exposed_Lower_leg','Small_Intestine-Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']

Tissues = ['Adipose', 'Adrenal_Gland', 'Bladder', 'Blood_Vessel', 'Blood', 'Brain', 'Breast', 'Cervix_Uteri', 'Colon', 'Esophagus', 'Fallopian_Tube', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pacreas', 'Pituitary', 'Prostate', 'Salivary_Gland', 'Skin', 'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina']

GTEx_idx = [[1, 2], [3], [7], [4, 5, 6], [23, 54], [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21], [24, 25], [26, 27], [28, 29, 30], [31], [32, 33], [34, 35], [36], [37], [39], [40], [41], [42], [43], [44], [38], [22, 45, 46], [47], [48], [49], [50], [51], [52], [53]]
TCGA_idx = [[0], [1], [2], [0], [7], [9, 15], [3], [4], [6, 24], [8], [0], [0], [11, 12, 13], [5, 16], [17, 18], [0], [0], [20], [21], [0], [23], [10], [26], [0], [0], [27], [28], [29], [31, 32], [0]]

name_idx_TCGA = np.loadtxt('name_idx_TCGA.txt', delimiter='\t')
name_idx_TCGA = [int(x) for x in name_idx_TCGA]
name_idx_GTEx = np.loadtxt('name_idx_GTEx.txt', delimiter='\t')
name_idx_GTEx = [int(x) for x in name_idx_GTEx]
gName = np.genfromtxt('gName_TGinsect.txt', dtype='str')  

Human_Clust = pd.read_csv(Human_Clust_path + 'HC_All_mir3.csv', skiprows = None, sep = ',', header = 0) # (395668, 34)
print(Human_Clust.shape)
print(Human_Clust.columns)

comb_BPs = pd.read_csv('comb_BPs.txt',sep = ',',header = None)
comb_BPs = list(comb_BPs.iloc[0])
print(len(comb_BPs))

with open(GO_path + 'BPtoGene2.json') as json_file:
    BPtoGene = json.load(json_file)
    print(BPtoGene.keys())# dict_keys(['BP_to_Gene'])
    BP_to_Gene = BPtoGene['BP_to_Gene']
    
with open(Human_Clust_path + 'Human_Clust_mir.json') as json_file:
    Human_Clust_JS = json.load(json_file)
    print(Human_Clust_JS.keys())# dict_keys(['Clustered_GO', 'Clustered_Gene', 'Clustered_TF'])
    
def FoundInBP(json_object, name, output_type, search_type = 'GO_index'):
    # json_object is the dataset that is going to be searched, in this case is BP_to_Gene
    # name is the terms that is looking for, e.g. BP = 'GO:0000002'
    # search_type is the searching itemes, e.g. search_type = 'GO_index' (default)
    # output_type is the output itemes, e.g. output_type = 'Genes'
    return [obj for obj in json_object if obj[search_type]==name][0][output_type]

def GO_to_gene(BP, data_set, output_type, Gene_Norm): # Pick genes in targeted GO term
    # BP is GO term index, e.g. BP = 'GO:0000002'
    # data_set is the json file that has GO and genes information, e.g., here is 'BP_to_Gene'
    # output_type is the item that is looking for, e.g. output_type = 'Genes', or 'TFs'
    # Gene_Norm is the gene list (in ORIGINAL order) of the RNA-seq data
    Gene_set = FoundInBP(data_set, BP, output_type)
    
    Gene_set = list(set(Gene_set))
    GO_Gene_idx = [i for i, item in enumerate(Gene_Norm) if item in set(Gene_set)]
    return(GO_Gene_idx)

def isInside(circle_x, circle_y, rad1, rad2, x, y): 
      
    # Compare radius of circle 
    # with distance of its center 
    # from given point 
    if ((x - circle_x) * (x - circle_x) + 
        (y - circle_y) * (y - circle_y) <= rad1 * rad2): 
        return True; 
    else: 
        return False; 
    
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

# def confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', **kwargs):
def confidence_ellipse(x, y, n_std=2.0):

    if len(x) != len(y):
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y) # covariance of x and y 
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1]) # Correlation coefficient of x,y
    
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))# Slope of ellipse
    if theta < (0):
        theta = (theta + 180)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x) 
    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)
    
    
    diff_mean = abs(mean_x - mean_y)
    diff_scale = abs(scale_x - scale_y)
    Distance = sum(np.square(x-mean_x) + np.square(y-mean_y))/len(x)
    
    return pearson, diff_mean, diff_scale, theta, Distance

'''
Calculate CorCoef2
'''

path_norm = 'Normalized_data/Norm_sepTiss2/Norm_mean_All/'
Chson_Tiss = [1,2,6,17,18,20,21,25,26,27] # 13 liver and 14 lung calculate seperately
Tiss = Chson_Tiss[ii]
lgn = []
print(Tiss, Tissues[Tiss])
G_idex = GTEx_idx[Tiss]
T_idex = TCGA_idx[Tiss]

file_TCGA = path_TPM + 'TCGA_TPM2/'
fn = T_idex[0]-1
TPM = pd.read_csv(file_TCGA + TCGA_project[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None, nrows=2)
print(TPM.shape)
Label_tmp = pd.read_csv('Labels/Label_' + TCGA_project[fn] + '2.csv')
caca = [i for i in np.where(Label_tmp['Sample Type']!='Primary solid Tumor')[0]]
caca = np.array(caca)
lgn.append(TPM.shape[1])

file_GTEx = path_TPM + 'GTEx_TPM2/'
fn = G_idex[0]-1
temp = pd.read_csv(file_GTEx + GTEx_tissues[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None, nrows=2)
print(temp.shape)
lgn.append(temp.shape[1])

Norm = pd.read_csv(path_norm + 'Norm_' + Tissues[Tiss]+'.txt',sep = '\t',index_col = [0],header=None)
Gene_Norm = list(Norm.index)
print(Norm.shape)

gty = ['Coef.','diff_mean','diff_scale','Slope','Distance']
xx = []      
xx = xx + ['Cancer.' + i for i in gty]
xx = xx + ['Health.' + i for i in gty]
xx =['BP1'] + ['BP2'] + ['BP_name1'] + ['BP_name2'] + xx
print(xx)
CorCoef = pd.DataFrame(columns=xx) # Creat file

for Pair in range(len(comb_BPs)):
# for Pair in range(5):
    if Pair % 1000 == 0:
        print(Pair)
    tmp = Human_Clust.iloc[comb_BPs[Pair]]
    GO1 = tmp[0]
    GO2 = tmp[2]
    tmp_cc = [GO1, GO2] 
    tmp_cc.append(tmp[1])
    tmp_cc.append(tmp[3])

    ## Check if the original data indexes are correct
    if (GO1!= Human_Clust_JS['Clustered_GO'][comb_BPs[Pair]]['clust_GOA'])|(GO2!= Human_Clust_JS['Clustered_GO'][comb_BPs[Pair]]['clust_GOB']):
        print('ERROR: GO terms ' + str(comb_BPs[Pair]) + ' different')
        break

    ## Corr. of all genes in BP pairs
    GeneSet1 = GO_to_gene(GO1, BP_to_Gene, 'Genes', Gene_Norm)
    GeneSet2 = GO_to_gene(GO2, BP_to_Gene, 'Genes', Gene_Norm)    

    meanGO1 = np.mean(Norm.iloc[GeneSet1,caca], axis = 0)
    meanGO2 = np.mean(Norm.iloc[GeneSet2,caca], axis = 0)
    [pearson, diff_mean, diff_scale, theta, Distance] = confidence_ellipse(meanGO1, meanGO2)
    cor_tmp = [pearson, diff_mean, diff_scale, theta, Distance]
    tmp_cc = tmp_cc + cor_tmp

    meanGO1 = np.mean(Norm.iloc[GeneSet1,lgn[0]:(sum(lgn))], axis = 0)
    meanGO2 = np.mean(Norm.iloc[GeneSet2,lgn[0]:(sum(lgn))], axis = 0)
    [pearson, diff_mean, diff_scale, theta, Distance] = confidence_ellipse(meanGO1, meanGO2)
    cor_tmp = [pearson, diff_mean, diff_scale, theta, Distance]
    tmp_cc = tmp_cc + cor_tmp

    CorCoef.loc[len(CorCoef)] = tmp_cc
print('Saving')
CorCoef.to_csv(Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_' + Tissues[Tiss] + '.csv', index=False)
print('Saved'+Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_' + Tissues[Tiss] + '.csv')
