import pandas as pd
import numpy as np
import sys
import os
import re  
import json
import random
import time
# import collections
# import math
# from scipy.stats import poisson
# from scipy import stats
# import matplotlib.pyplot as plt
# from matplotlib.patches import Ellipse
# import matplotlib.transforms as transforms
# from sklearn import preprocessing

print('starting')
ii = sys.argv[1]
ii = int(ii)-1
print(ii)


Human_Clust_path = 'Human_ClustGO/'
GO_path = 'GO_annotation/'
path_TPM = 'TPM/'
T_path = 'TPM/TCGA_TPM2/'
path_norm = 'Normalized_data/Norm_sepTiss2/Norm_mean_All/'

Human_Clust = pd.read_csv(Human_Clust_path + 'HC_All_mir3.csv', skiprows = None, sep = ',', header = 0) # (395668, 34)
print(Human_Clust.shape)
print(Human_Clust.columns)

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

caca = range(54)
hyhy = range(54,87)
'''
Calculate mean of all samples in one tissue
'''
Norm = pd.read_csv(path_norm + 'Norm_Bladder.txt',sep = '\t',index_col = [0],header=None)
Gene_Norm = list(Norm.index)
print(Norm.shape)
Norm = pd.read_csv(path_norm + 'GeneMean.csv', skiprows = None, sep = '\t', header = None)
xx =['BP1'] + ['BP2'] + ['BP_name1'] + ['BP_name2'] + ['Coef.healthy'] + ['Coef.cancer']
print(xx)
CorCoef = pd.DataFrame(columns=xx) # Creat file

for Pair in range(Human_Clust.shape[0]):
#for Pair in range(ii*10000,(ii+1)*10000):
    if Pair % 1000 == 0:
        print(Pair)
    tmp = Human_Clust.iloc[Pair]
    GO1 = tmp[0]
    GO2 = tmp[1]
    tmp_cc = [GO1, GO2] 
    tmp_cc.append(tmp[1])
    tmp_cc.append(tmp[3])

    ## Check if the original data indexes are correct
    if (GO1!= Human_Clust_JS['Clustered_GO'][Pair]['clust_GOA'])|(GO2!= Human_Clust_JS['Clustered_GO'][Pair]['clust_GOB']):
        print('ERROR: GO terms ' + str(Pair) + ' different')
        break

    ## Corr. of all genes in BP pairs
    GeneSet1 = GO_to_gene(GO1, BP_to_Gene, 'Genes', Gene_Norm)
    GeneSet2 = GO_to_gene(GO2, BP_to_Gene, 'Genes', Gene_Norm)  
    # print(len(GeneSet1),len(GeneSet2))
    if Human_Clust['Pair_Gene'][Pair]!=0:
        rmv = []
        for jj in Human_Clust_JS['Clustered_Gene'][Pair]['clust_GeneA']:
            kk = np.where(np.array(Gene_Norm)==jj)[0]
            if len(kk)!=0:
                rmv.append(kk[0])
        GeneSet1 = [i for i in GeneSet1 if i not in rmv]
        rmv = []
        for jj in Human_Clust_JS['Clustered_Gene'][Pair]['clust_GeneB']:
            kk = np.where(np.array(Gene_Norm)==jj)[0]
            if len(kk)!=0:
                rmv.append(kk[0])
        GeneSet2 = [i for i in GeneSet2 if i not in rmv]
    # print(len(GeneSet1),len(GeneSet2))

    olGene = list(set(GeneSet1) & set(GeneSet2))
    rmvA = []
    rmvB = []
    if len(olGene) >0:
        for k in range(len(olGene)):
            rmvA.append(GeneSet1.index(olGene[k]))
            rmvB.append(GeneSet2.index(olGene[k]))
    GeneSet11 = list(set(range(len(GeneSet1))) - set(rmvA))
    GeneSet11 = list(np.array(GeneSet1)[GeneSet11])
    GeneSet22= list(set(range(len(GeneSet2))) - set(rmvB))
    GeneSet22 = list(np.array(GeneSet2)[GeneSet22])
    
    meanGO1 = np.mean(Norm.iloc[GeneSet11,caca], axis = 0)
    meanGO2 = np.mean(Norm.iloc[GeneSet22,caca], axis = 0)
    cov = np.cov(meanGO1, meanGO2) # covariance of x and y 
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    tmp_cc = tmp_cc + [pearson]

    meanGO1 = np.mean(Norm.iloc[GeneSet11,hyhy], axis = 0)
    meanGO2 = np.mean(Norm.iloc[GeneSet22,hyhy], axis = 0)
    cov = np.cov(meanGO1, meanGO2) # covariance of x and y 
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    tmp_cc = tmp_cc + [pearson]

    CorCoef.loc[len(CorCoef)] = tmp_cc


print('Saving')
CorCoef.to_csv(Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_All2_'+str(ii)+'.csv', index=False)
print('Saved'+Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_All2_'+str(ii)+'.csv')
