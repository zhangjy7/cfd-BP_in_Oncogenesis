# This is a python3 file
print('starting')
import sys
import numpy as np
import os
import re
import pandas as pd
import json
from random import sample 
import time
import collections
from sklearn.neighbors import KNeighborsClassifier
from sklearn import decomposition
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC

print('starting')
ii = sys.argv[1]
ii = int(ii)-1
print(ii)


TCGA_project = ['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']
GTEx_tissues = ['Adipose-Subcutaneous', 'Adipose-Visceral', 'Adrenal_Gland', 'Artery-Aorta','Artery-Coronary', 'Artery-Tibial', 'Bladder', 'Brain-Amygdala','Brain-Anterior_cingulate_cortex', 'Brain-Caudate', 'Brain-Cerebellar_Hemisphere', 'Brain-Cerebellum','Brain-Cortex','Brain-Frontal_Cortex', 'Brain-Hippocampus', 'Brain-Hypothalamus', 'Brain-Nucleus_accumbens', 'Brain-Putamen','Brain-Spinal_cord', 'Brain-Substantia_nigra', 'Breast-Mammary_Tissue', 'Cells-Cultured_fibroblasts', 'Cells-EBV-transformed_lymphocytes', 'Cervix-Ectocervix', 'Cervix-Endocervix', 'Colon-Sigmoid', 'Colon-Transverse','Esophagus-Gastroesophageal_Junction', 'Esophagus-Mucosa','Esophagus-Muscularis', 'Fallopian Tube', 'Heart-Atrial_Appendage', 'Heart-Left_Ventricle', 'Kidney-Cortex', 'Kidney-Medulla', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle-Skeletal', 'Nerve-Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin-Not_Sun_Exposed_Suprapubic', 'Skin-Sun_Exposed_Lower_leg','Small_Intestine-Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']

### Combined sample and PCA/HC analysis based on tissues
Tissues = ['Adipose', 'Adrenal_Gland', 'Bladder', 'Blood_Vessel', 'Blood', 'Brain', 'Breast', 'Cervix_Uteri', 'Colon', 'Esophagus', 'Fallopian_Tube', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pacreas', 'Pituitary', 'Prostate', 'Salivary_Gland', 'Skin', 'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina']

GTEx_idx = [[1, 2], [3], [7], [4, 5, 6], [23, 54], [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21], [24, 25], [26, 27], [28, 29, 30], [31], [32, 33], [34, 35], [36], [37], [39], [40], [41], [42], [43], [44], [38], [22, 45, 46], [47], [48], [49], [50], [51], [52], [53]]
TCGA_idx = [[0], [1], [2], [0], [7], [9, 15], [3], [4], [6, 24], [8], [0], [0], [11, 12, 13], [5, 16], [17, 18], [0], [0], [20], [21], [0], [23], [10], [26], [0], [0], [27], [28], [29], [31, 32], [0]]

GO_path = 'GO_annotation/'
Human_Clust_path = 'Human_ClustGO/'
path_label = 'Labels/'
path_norm = 'Normalized_data/Norm_sepTiss2/Norm_mean_All/'
T_path = 'TPM/TCGA_TPM2/'
G_path = 'TPM/GTEx_TPM2/'


with open(GO_path + 'BPtoGene2.json') as json_file:
    BPtoGene = json.load(json_file)
    print(BPtoGene.keys())# dict_keys(['BP_to_Gene'])
    BP_to_Gene = BPtoGene['BP_to_Gene']

## Import Normalized sample TPM ####
    
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

def isInside(circle_x, circle_y, rad1, x, y): 
      
    # Compare radius of circle 
    # with distance of its center 
    # from given point 
    if ((x - circle_x) * (x - circle_x) + 
        (y - circle_y) * (y - circle_y) <= rad1): 
        return True; 
    else: 
        return False; 

classifiers = [KNeighborsClassifier(5)]
names = ["Nearest Neighbors"]

'''
Calcualte distance <dis_hd> to cancer and healthy regions from chosen BP (T_BP)
'''
thresh_PCA = pd.read_csv('thresh_tissuesPCA_mir_learn.csv', skiprows = None, header = 0, index_col=0 )

Chson_Tiss = [1,2,6,17,18,20,21,25,26,27] # Except Lung
Tiss = Chson_Tiss[ii]
TG= ['KICH','KIRC','KIRP']
# TG= ['LUAD','LUSC']
# Tiss = 12
nms = Tissues[Tiss]+'-'+TG[ii]
CorCoef = pd.read_csv(Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_' + nms+ '_learn.csv', skiprows = None, sep = ',', header = 0)
print(Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_' + Tissues[Tiss] + '_learn.csv')
Y = CorCoef.values[:,4:14]
X = preprocessing.scale(Y,axis=0)
pca = decomposition.PCA(n_components=5)
pca.fit(X) 
PCAX = pca.transform(X)
# plt.scatter(PCAX[:,0], PCAX[:,1], c='grey', marker='.',s = 1)

thresh = np.where(thresh_PCA.index==nms)
thresh_si = thresh_PCA.values[thresh][0]

T_BP = []
circle_x = thresh_si[0]
circle_y = thresh_si[1]
rad1 = thresh_si[2]
for i in range(0,PCAX.shape[0]):
    if isInside(circle_x, circle_y, rad1, PCAX[i,0], PCAX[i,1])==False:
        T_BP.append(i)
print(len(T_BP))
circle_x = thresh_si[3]
circle_y = thresh_si[4]
rad1 = thresh_si[5]
for i in range(0,PCAX.shape[0]):
    if isInside(circle_x, circle_y, rad1, PCAX[i,0], PCAX[i,1])==False:
        T_BP.append(i)
print(len(T_BP))
TP_tmp = np.where(abs(CorCoef['Cancer.Coef.']-CorCoef['Health.Coef.'])>0.5)[0]
T_BP = T_BP+ list(TP_tmp)
T_BP = list(set(T_BP))
print(len(T_BP))

if len(T_BP)!=thresh_si[-1]:
    print('Error')
    sys.exit()

print(Tiss, Tissues[Tiss])
Norm = pd.read_csv( path_norm+'Norm_' + Tissues[Tiss] + '.txt', skiprows = None, sep = '\t', header = None, index_col = [0])# 
print(Norm.shape)
Gene_Norm = list(Norm.index)
Norm_learn = pd.read_csv( 'Normalized_data/Norm_split9-1/Norm_' + nms + '_learn.csv', skiprows = None, sep = ',', header = 0, index_col = None)# 
print(Norm.shape)
Norm_veri = pd.read_csv( 'Normalized_data/Norm_split9-1/Norm_' + nms + '_veri.csv', skiprows = None, sep = ',', header = 0, index_col = None)# 
print(Norm.shape)


### Seperate different normal/cancer tissue #start###
T_idx = TCGA_idx[Tiss]
G_idx = GTEx_idx[Tiss] 

lgn = []
fn = T_idx[ii] - 1
print(T_path + TCGA_project[fn] +'.txt.TPM')   
temp = pd.read_csv(T_path + TCGA_project[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None, nrows = 2)
lgn.append(temp.shape[1])

fn = G_idx[0] - 1
temp = pd.read_csv(G_path + GTEx_tissues[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None, nrows = 2)
print(GTEx_tissues[fn], temp.shape)
lgn.append(temp.shape[1])

caca = np.where(np.array([int(i) for i in Norm_learn.columns])<lgn[0])[0]
hyhy = np.where(np.array([int(i) for i in Norm_learn.columns])>=lgn[0])[0]
y_learn = np.array([1]*len(caca) + [0]*len(hyhy))
dis_hd_learn = [1]*len(caca) + [0]*len(hyhy)
print(len(dis_hd_learn))

caca = np.where(np.array([int(i) for i in Norm_veri.columns])<lgn[0])[0]
hyhy = np.where(np.array([int(i) for i in Norm_veri.columns])>=lgn[0])[0]
y_veri = np.array([1]*len(caca) + [0]*len(hyhy))
dis_hd_veri = [1]*len(caca) + [0]*len(hyhy)
print(len(dis_hd_veri))

atmp = pd.DataFrame(columns=['GO1','GO2','GO1_name','GO2_name','score_learn','score_veri'])
a = 0
start_time = time.time()
# for Pair in T_BP:
for Pair in T_BP[:5]:
    a = a + 1
    if a % 1000==0:
        print(a)
    tmp = CorCoef.iloc[Pair]
    GO1 = tmp[0]
    GO2 = tmp[1]
    aatmp = [tmp[0], tmp[1], tmp[2], tmp[3]]

    GeneSet1 = GO_to_gene(GO1, BP_to_Gene, 'Genes', Gene_Norm)
    GeneSet2 = GO_to_gene(GO2, BP_to_Gene, 'Genes', Gene_Norm)
    # Learn
    Genetmp = Norm_learn.values[(GeneSet1 + GeneSet2)]
    
    X = np.transpose(Genetmp)
    X_train, X_test, y_train, y_test = train_test_split(X, y_learn, test_size=0.33, random_state=42)
    clf = SVC(kernel="linear", C=1)
    clf.fit(X_train, y_train)
    score = clf.score(X_test, y_test)
    aatmp.append(score)
    
    Genetmp_veri = Norm_veri.values[(GeneSet1 + GeneSet2)]
    X_veri = np.transpose(Genetmp_veri)
    X_train, X_test, y_train, y_test = train_test_split(X_veri, y_veri, test_size=0.33, random_state=42)
    score = clf.score(X_test, y_test)
    
    Z = clf.predict(X)
    health = np.where((Z == 0) & (y_learn ==0))[0]
    cancer = np.where((Z == 1) & (y_learn ==1))[0]
    
    if (len(cancer)==0) | (len(health)==0):
        aatmp.append('na')
        atmp.loc[len(atmp)] = aatmp
        dis_hd_learn = np.column_stack((dis_hd_learn,np.array(dis_hd_learn.shape[0]*[0])))
        dis_hd_learn = np.column_stack((dis_hd_learn,np.array(dis_hd_learn.shape[0]*[0])))
    else:
        aatmp.append(score)
        atmp.loc[len(atmp)] = aatmp
        
        healthy_tmp = X[health,:]
        dish = np.linalg.norm(np.transpose(Genetmp)-healthy_tmp[0,:],axis = 1)
        for tt in range(1,healthy_tmp.shape[0]):
            dis_tmp = np.linalg.norm(np.transpose(Genetmp)-healthy_tmp[tt,:],axis = 1)
            dish = dish + dis_tmp
        dis_hd_learn = np.column_stack((dis_hd_learn,dish/len(health)))
        cancerh_smp = X[cancer,:]
        disc = np.linalg.norm(np.transpose(Genetmp)-cancerh_smp[0,:],axis = 1)
        for tt in range(1,cancerh_smp.shape[0]):
            dis_tmp = np.linalg.norm(np.transpose(Genetmp)-cancerh_smp[tt,:],axis = 1)
            disc = disc + dis_tmp
        dis_hd_learn = np.column_stack((dis_hd_learn,-disc/len(cancer)))

        # Veri
        dish = np.linalg.norm(np.transpose(Genetmp_veri)-healthy_tmp[0,:],axis = 1)
        for tt in range(1,healthy_tmp.shape[0]):
            dis_tmp = np.linalg.norm(np.transpose(Genetmp_veri)-healthy_tmp[tt,:],axis = 1)
            dish = dish + dis_tmp
        dis_hd_veri = np.column_stack((dis_hd_veri,dish/len(health)))
        disc = np.linalg.norm(np.transpose(Genetmp_veri)-cancerh_smp[0,:],axis = 1)
        for tt in range(1,cancerh_smp.shape[0]):
            dis_tmp = np.linalg.norm(np.transpose(Genetmp_veri)-cancerh_smp[tt,:],axis = 1)
            disc = disc + dis_tmp
        dis_hd_veri = np.column_stack((dis_hd_veri,-disc/len(cancer)))
    
np.savetxt(Human_Clust_path + 'Distances/dis_hd_mir_'+nms+'_learn.csv', dis_hd_learn, delimiter=",")
np.savetxt(Human_Clust_path + 'Distances/dis_hd_mir_'+nms+'_veri.csv', dis_hd_veri, delimiter=",")
atmp.to_csv(Human_Clust_path + 'Distances/dis_hd_mir_'+nms+'_go.csv', index=False)
print("--- %s seconds ---" % (time.time() - start_time))
