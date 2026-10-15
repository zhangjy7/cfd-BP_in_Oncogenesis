#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings("ignore")
import json
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing, decomposition
from sklearn.svm import SVC
from sklearn.metrics import ConfusionMatrixDisplay

# ---- reference tables (notebook cell 0) ----
TCGA_project = ['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LAML', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']
GTEx_tissues = ['Adipose-Subcutaneous', 'Adipose-Visceral', 'Adrenal_Gland', 'Artery-Aorta','Artery-Coronary', 'Artery-Tibial', 'Bladder', 'Brain-Amygdala','Brain-Anterior_cingulate_cortex', 'Brain-Caudate', 'Brain-Cerebellar_Hemisphere', 'Brain-Cerebellum','Brain-Cortex','Brain-Frontal_Cortex', 'Brain-Hippocampus', 'Brain-Hypothalamus', 'Brain-Nucleus_accumbens', 'Brain-Putamen','Brain-Spinal_cord', 'Brain-Substantia_nigra', 'Breast-Mammary_Tissue', 'Cells-Cultured_fibroblasts', 'Cells-EBV-transformed_lymphocytes', 'Cervix-Ectocervix', 'Cervix-Endocervix', 'Colon-Sigmoid', 'Colon-Transverse','Esophagus-Gastroesophageal_Junction', 'Esophagus-Mucosa','Esophagus-Muscularis', 'Fallopian Tube', 'Heart-Atrial_Appendage', 'Heart-Left_Ventricle', 'Kidney-Cortex', 'Kidney-Medulla', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle-Skeletal', 'Nerve-Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin-Not_Sun_Exposed_Suprapubic', 'Skin-Sun_Exposed_Lower_leg','Small_Intestine-Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']
Tissues = ['Adipose', 'Adrenal_Gland', 'Bladder', 'Blood_Vessel', 'Blood', 'Brain', 'Breast', 'Cervix_Uteri', 'Colon', 'Esophagus', 'Fallopian_Tube', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pacreas', 'Pituitary', 'Prostate', 'Salivary_Gland', 'Skin', 'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina']
GTEx_idx = [[1, 2], [3], [7], [4, 5, 6], [23, 54], [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [21], [24, 25], [26, 27], [28, 29, 30], [31], [32, 33], [34, 35], [36], [37], [39], [40], [41], [42], [43], [44], [38], [22, 45, 46], [47], [48], [49], [50], [51], [52], [53]]
TCGA_idx = [[0], [1], [2], [0], [7], [9, 15], [3], [4], [6, 24], [8], [0], [0], [11, 12, 13], [5, 16], [17, 18], [0], [0], [20], [21], [0], [23], [10], [26], [0], [0], [27], [28], [29], [31, 32], [0]]

# ---- paths (notebook cell 0) ----
Human_Clust_path = 'Human_ClustGO/'
GO_path = 'GO_annotation/'
path_norm = 'Norm_mean_All/'
path_label = 'TCGA_clinic/'
T_path = 'TPM/TCGA_TPM2/'
G_path = 'TPM/GTEx_TPM2/'

# GO term -> genes table (needed by Figure 5B)
with open(GO_path + 'BPtoGene2.json') as json_file:
    BP_to_Gene = json.load(json_file)['BP_to_Gene']

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


# ====================================================================================================
# BLOCK 1 - Cancer identification: PCA -> SVM, confusion matrices
# Fig 4, Fig S5 (scatter rows: train / test / validation / validation + grey biopsy-normal), Fig S6 (confusion
# matrices), Fig 5A (Prostate).  Body = notebook cell 36 (cells 38 and 40 differ only in SVM settings and file names).
# ====================================================================================================
SEED = 16

# (index in Tissues, file-name stem, index in TCGA_idx[Tiss], title, SVM kernel, SVM C)
# SVM kernel / C per cancer are the values recorded in the notebook comments (cells 36, 38, 40).
# Testis has no recorded setting; the notebook default (linear, C = 1) is used.
RUNS = [
    ( 1, 'Adrenal_Gland', 0, 'Adrenal_Gland', 'poly',   1),
    ( 2, 'Bladder',       0, 'Bladder',       'poly',   1),
    ( 6, 'Breast',        0, 'Breast',        'poly',  10),
    (12, 'Kidney-KICH',   0, 'Kidney-KICH',   'linear', 1),
    (12, 'Kidney-KIRC',   1, 'Kidney-KIRC',   'rbf',    3),
    (12, 'Kidney-KIRP',   2, 'Kidney-KIRP',   'poly',  15),
    (13, 'Liver',         1, 'Liver-LIHC',    'rbf',  100),   # TCGA_idx[13] = [CHOL, LIHC]; LIHC is the second
    (14, 'Lung-LUAD',     0, 'Lung-LUAD',     'poly',  10),
    (14, 'Lung-LUSC',     1, 'Lung-LUSC',     'poly',  10),
    (17, 'Ovary',         0, 'Ovary',         'poly',   1),
    (18, 'Pacreas',       0, 'Pancreas',      'poly', 100),   # 'Pacreas' is the spelling used in the file names
    (20, 'Prostate',      0, 'Prostate',      'linear', 1),   # -> Figure 5A
    (25, 'Stomach',       0, 'Stomach',       'linear', 1),
    (26, 'Testis',        0, 'Testis',        'linear', 1),
    (27, 'Thyroid',       0, 'Thyroid',       'linear', 1),
]
ONLY = None      # e.g. ['Pancreas', 'Prostate'] to run a subset

def cancer_identification(Tiss, nms, ii, label, kernel, C):
    random.seed(SEED)   # the notebook split was unseeded; seeding makes the 2:1 train/test split repeatable
    dis_hd_learn = pd.read_csv(Human_Clust_path + 'Distances/dis_hd_mir_'+nms+'_learn.csv', skiprows = None, sep = ',', header = None)
    print(dis_hd_learn.shape)
    dis_hd_veri = pd.read_csv(Human_Clust_path + 'Distances/dis_hd_mir_'+nms+'_veri.csv', skiprows = None, sep = ',', header = None)
    print(dis_hd_veri.shape)
    
    Norm_learn = pd.read_csv( 'Norm_split9-1/Norm_' + nms + '_learn.csv', skiprows = None, sep = ',', header = 0, index_col = None)# 
    print(Norm_learn.shape)
    Norm_veri = pd.read_csv( 'Norm_split9-1/Norm_' + nms + '_veri.csv', skiprows = None, sep = ',', header = 0, index_col = None)# 
    print(Norm_veri.shape)
    
    T_idex = TCGA_idx[Tiss]
    fn = T_idex[ii]-1
    Label_tmp = pd.read_csv(path_label+'Labels/Label_' + TCGA_project[fn] + '2.csv')
    print(TCGA_project[fn],Label_tmp.shape)
    lgn = Label_tmp.shape[0]

    # 1:3 test and train sample splits
    sampling = random.sample(range(dis_hd_learn.shape[0]), k=round(dis_hd_learn.shape[0]*0.667))
    dis_hd_train = dis_hd_learn.values[sampling,:]
    dis_hd_test = np.delete(dis_hd_learn.values, sampling, 0)
    y_train = dis_hd_train[:,0]
    y_test = dis_hd_test[:,0]
    
    Y = dis_hd_train[:,1::]
    X = preprocessing.scale(Y,axis=0)
    # X = X.transpose()
    pca = decomposition.PCA(n_components=5)
    pca.fit(X)  
    PCAX1 = pca.transform(X)
    print(PCAX1.shape)
    
    cols = ['green','red','grey']
    col1 = np.where(y_train==0)[0]
    col2 = np.where(y_train==1)[0]
    col3 = np.where(y_train==-1)[0]
    print(len(col1),len(col2),len(col3))
    data = (col1,col2,col3)
    
    clf = SVC(kernel=kernel, degree=3, C=C).fit(PCAX1[:, :2], y_train)   # kernel / C: see RUNS below
    print(clf.score(PCAX1[:,:2],y_train))
    np.set_printoptions(precision=2)
# Plot non-normalized confusion matrix
    titles_options = [
        (label+"-train", None),
        (label+"-train", "true"),
    ]
    for title, normalize in titles_options:
        disp = ConfusionMatrixDisplay.from_estimator(
            clf,
            PCAX1[:, :2],
            y_train,
    #         display_labels=class_names,
            cmap=plt.cm.Blues,
            normalize=normalize,
            colorbar=False
        )
        disp.ax_.set_title(title)

        print(title)
        print(disp.confusion_matrix)
    plt.show()
        # create a mesh to plot in
    h = 4 # step size in the mesh
    x_min, x_max = PCAX1[:, 0].min() - 20, PCAX1[:, 0].max() + 20
    y_min, y_max = PCAX1[:, 1].min() - 20, PCAX1[:, 1].max() + 20
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))

    # title for the plots
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    cols = ['limegreen','lightcoral','grey']
    plt.figure(figsize=(10,10))
    plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.2)
    plt.title(label)
    plt.xlim([min(PCAX1[:,0])-20, max(PCAX1[:,0])+20])
    plt.ylim([min(PCAX1[:,1])-20, max(PCAX1[:,1])+20])
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX1[tmp,0], PCAX1[tmp,1],c = col,alpha =1)

    # Plot test sample on to the train area    
    plt.figure(figsize=(10,10))
    plt.title(label)
    plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.2)
    plt.xlim([min(PCAX1[:,0])-20, max(PCAX1[:,0])+20])
    plt.ylim([min(PCAX1[:,1])-20, max(PCAX1[:,1])+20])
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX1[tmp,0], PCAX1[tmp,1],c = col,alpha = 0.8)

    Y_test = dis_hd_test[:,1::]
    X_test = preprocessing.scale(Y_test,axis=0)
    PCAX = pca.transform(X_test)
    print(PCAX.shape)
    print('test')
    print(clf.score(PCAX[:,:2],y_test))
    
    cols = ['forestgreen','red','grey']
    col1 = np.where(y_test==0)[0]
    col3 = np.where(y_test==-1)[0]
    col2 = np.where(y_test==1)[0]
    print(len(col1),len(col2),len(col3))
    data = (col1,col2,col3)
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX[tmp,0], PCAX[tmp,1],c = col,alpha =1)
    plt.title(label+' '+str(round(clf.score(PCAX[:,:2],y_test),2)))
    
    np.set_printoptions(precision=2)
# Plot non-normalized confusion matrix
    titles_options = [
        (label+"-test", None),
        (label+"-test", "true"),
    ]
    for title, normalize in titles_options:
        disp = ConfusionMatrixDisplay.from_estimator(
            clf,
            PCAX[:,:2],
            y_test,
    #         display_labels=class_names,
            cmap=plt.cm.Blues,
            normalize=normalize,
            colorbar=False
        )
        disp.ax_.set_title(title)

        print(title)
        print(disp.confusion_matrix)

    plt.show()

    # plot veri samples on learn samples
    y_learn = dis_hd_learn.values[:,0]
    y_veri = dis_hd_veri.values[:,0]
    Y = dis_hd_learn.values[:,1::]
    X = preprocessing.scale(Y,axis=0)
    PCAX_lr = pca.transform(X)
    print(PCAX_lr.shape)

    cols = ['limegreen','lightcoral','grey']
    col1 = np.where(y_train==0)[0]
    col2 = np.where(y_train==1)[0]
    col3 = np.where(y_train==-1)[0]
    print(len(col1),len(col2),len(col3))
    data = (col1,col2,col3)

    plt.figure(figsize=(10,10))
    plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.2)
    plt.xlim([min(PCAX1[:,0])-20, max(PCAX1[:,0])+20])
    plt.ylim([min(PCAX1[:,1])-20, max(PCAX1[:,1])+20])
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX1[tmp,0], PCAX1[tmp,1],c = col,alpha = 0.5)

    Y_veri = dis_hd_veri.values[:,1::]
    X_veri = preprocessing.scale(Y_veri,axis=0)
    PCAX = pca.transform(X_veri)
    print(PCAX.shape)
    print('veri')
    
    col1 = np.where(y_veri==0)[0]
    col2 = np.where(y_veri==1)[0]
    stage00 = np.array([int(Norm_veri.columns[i]) for i in range(len(Norm_veri.columns))])
    stage00 = stage00[np.where(stage00<lgn)[0]]
    col3 = np.where(Label_tmp['Sample Type'][stage00]!='Primary solid Tumor')[0]
    col1 = np.array([i for i in col1 if i not in col3])
    col2 = np.array([i for i in col2 if i not in col3])
    
    print(len(col1),len(col2),len(col3))
    data = (col1,col2)
    cols = ['darkgreen','purple']
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX[tmp,0], PCAX[tmp,1],c = col,alpha =1)
    plt.title(label+' '+str(round(clf.score(np.delete(PCAX,col3,0)[:,:2],np.delete(y_veri,col3)),2)))
    titles_options = [
        (label+"-veri", None),
        (label+"-veri", "true"),
    ]
    for title, normalize in titles_options:
        disp = ConfusionMatrixDisplay.from_estimator(
            clf,
            np.delete(PCAX,col3,0)[:,:2],
            np.delete(y_veri,col3),
    #         display_labels=class_names,
            cmap=plt.cm.Blues,
            normalize=normalize,
            colorbar=False
        )
        disp.ax_.set_title(title)

        print(title)
        print(disp.confusion_matrix)

    plt.show()
    
    # add grey dots here
    cols = ['limegreen','lightcoral','grey']
    col1 = np.where(y_train==0)[0]
    col2 = np.where(y_train==1)[0]
    col3 = np.where(y_train==-1)[0]
    print(len(col1),len(col2),len(col3))
    data = (col1,col2,col3)

    plt.figure(figsize=(10,10))
    plt.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.2)
    plt.xlim([min(PCAX1[:,0])-20, max(PCAX1[:,0])+20])
    plt.ylim([min(PCAX1[:,1])-20, max(PCAX1[:,1])+20])
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX1[tmp,0], PCAX1[tmp,1],c = col,alpha = 0.5)
    
    col1 = np.where(y_veri==0)[0]
    col2 = np.where(y_veri==1)[0]
    stage00 = np.array([int(Norm_veri.columns[i]) for i in range(len(Norm_veri.columns))])
    stage00 = stage00[np.where(stage00<lgn)[0]]
    col3 = np.where(Label_tmp['Sample Type'][stage00]!='Primary solid Tumor')[0]
    
    print(clf.score(np.delete(PCAX,col3,0)[:,:2],np.delete(y_veri,col3)))
    cols = ['darkgreen','purple','grey']
    data = (col1,col2,col3)
    for col,tmp in zip(cols,data):
        plt.scatter(PCAX[tmp,0], PCAX[tmp,1],c = col,alpha =1)
    plt.title(label+' '+str(round(clf.score(np.delete(PCAX,col3,0)[:,:2],np.delete(y_veri,col3)),2)))
    plt.show()


for (Tiss, nms, ii, label, kernel, C) in RUNS:
    if ONLY is None or label in ONLY:
        print('=====', label)
        cancer_identification(Tiss, nms, ii, label, kernel, C)


# Representative tissue: Liver 
Tiss = 13
CorCoef = pd.read_csv(Human_Clust_path + 'CorCoef/CorCoef_cbBP_mir_' + Tissues[Tiss] + '_learn.csv', skiprows = None, sep = ',', header = 0)
Ca_large = np.where(abs(CorCoef['Cancer.Coef.']-CorCoef['Health.Coef.'])>0.5)[0]     # pairs whose R changes by > 0.5
print(len(Ca_large))

# PCA of the five features of each BP pair
pca = decomposition.PCA(n_components=5)
Y = CorCoef.values[:,4:14]
X = preprocessing.scale(Y,axis=0)
pca.fit(X)
PCAX = pca.transform(X)
plt.figure()
plt.scatter(PCAX[:,0], PCAX[:,1], c='grey', marker='.',s = 1)
plt.scatter(PCAX[Ca_large,0], PCAX[Ca_large,1], c='blue', marker='.',s = 1)

# distribution of the change in correlation coefficient
plt.figure()
plt.hist(CorCoef['Cancer.Coef.']-CorCoef['Health.Coef.'])

Tiss = 20
print(Tiss, Tissues[Tiss])
Norm = pd.read_csv(path_norm + 'Norm_' + Tissues[Tiss] + '.txt', skiprows = None, sep = '\t', header = None, index_col = [0])# import normalzied RNA-seq data
print(Norm.shape)
Gene_Norm = list(Norm.index)
T_idx = TCGA_idx[Tiss]
G_idx = GTEx_idx[Tiss] 

lgn = []
fn = T_idx[0] - 1
print(T_path + TCGA_project[fn] +'.txt.TPM')   
TPM = pd.read_csv(T_path + TCGA_project[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None)
Label_tmp = pd.read_csv(path_label+'Labels/Label_TCGA-PRAD2.csv')
lgn.append(TPM.shape[1])
Gleason = pd.read_csv(path_label+'TCGA-PRAD/TCGA-PRAD_clinical.csv')
Gleason_index = []
for i in Label_tmp['Sample_index']:
    tmp = np.where(Gleason['submitter_id']==i)[0]
    Gleason_index.append(int(tmp[0]))   # int(array) breaks on numpy >= 1.25
Gleason_score1 = [int(Gleason['primary_gleason_grade'][ii].split(" ")[1]) for ii in Gleason_index]
Gleason_scoret = list(Gleason_score1)
for ii in range(len(Gleason_scoret)):
    if Gleason_scoret[ii]<=3:
        Gleason_scoret[ii] = 3
    if Gleason_scoret[ii]>4:
        Gleason_scoret[ii] = 5
Tumor_sample = np.where(Label_tmp['Sample Type']=='Primary solid Tumor')[0]
Gleason_score = [Gleason_scoret[i] for i in Tumor_sample]
# Gleason_score = Gleason['primary_gleason_grade'][Gleason_index]

fn = G_idx[0] - 1
temp = pd.read_csv(G_path + GTEx_tissues[fn]+'.txt.TPM', skiprows = None, sep = '\t', header = None)
print(GTEx_tissues[fn], temp.shape)
lgn.append(temp.shape[1])
print(lgn)

Pri3 = np.where(np.array(Gleason_score)==3)[0]
Pri4 = np.where(np.array(Gleason_score)==4)[0]
Pri5 = np.where(np.array(Gleason_score)==5)[0]
caca = [i for i in np.where(Label_tmp['Sample Type']=='Primary solid Tumor')[0]]
caca = np.array(caca)
stage001 = [i for i in np.where(Label_tmp['Sample Type']!='Primary solid Tumor')[0]]
stage001 = np.array(stage001)
hyhy = [i for i in range(lgn[0],sum(lgn))]
hyhy = np.array(hyhy)
data = (Pri3, Pri4, Pri5,hyhy,caca,stage001)
cols = ['orange','orchid','purple','green','red','grey']

# primary Gleason grade per sample (<=3 -> 3, 4, >4 -> 5), non-primary-tumour samples dropped
Gleason_score1 = [int(Gleason['primary_gleason_grade'][ii].split(" ")[1]) for ii in Gleason_index]
Gleason_scoret = list(Gleason_score1)
for ii in range(len(Gleason_scoret)):
    if Gleason_scoret[ii]<=3:
        Gleason_scoret[ii] = 3
    if Gleason_scoret[ii]>4:
        Gleason_scoret[ii] = 5
Gleason_score1 = [i for j, i in enumerate(Gleason_scoret) if j not in stage001]

Norm1 = np.delete(Norm.values, stage001,1)

Pri3 = np.where(np.array(Gleason_score1)==3)[0]
Pri4 = np.where(np.array(Gleason_score1)==4)[0]
Pri5 = np.where(np.array(Gleason_score1)==5)[0]
caca1 = np.array(range(len(Gleason_score1)))
hyhy1 = hyhy-len(stage001)
cols = ['orange','red','purple','green']
data = [Pri3, Pri4, Pri5, hyhy1]

# ---- Figure 5B: angiogenesis (GO:0001525, 226 genes) vs regulation of mitotic cell cycle (GO:0007346, 41 genes) ----
# one scatter per group: primary Gleason grade <=3 (orange), 4 (red), 5 (purple) and healthy (green); title = Pearson R
GO1, GO2 = "GO:0001525", "GO:0007346"
print(GO1, GO2)
print(FoundInBP(BP_to_Gene, GO1,'GO_Name', search_type = 'GO_index'))
print(FoundInBP(BP_to_Gene, GO2,'GO_Name', search_type = 'GO_index'))
GeneSet1 = GO_to_gene(GO1, BP_to_Gene, 'Genes', Gene_Norm)
GeneSet2 = GO_to_gene(GO2, BP_to_Gene, 'Genes', Gene_Norm)
meanGO1 = np.mean(Norm1[GeneSet1,:], axis = 0)
meanGO2 = np.mean(Norm1[GeneSet2,:], axis = 0)
cols = ['orange','red','purple','green']
fig, axs = plt.subplots(1, 4,figsize=(15,3.5))
for col, dd, ii in zip(cols, data, range(4)):
    cc1 = round(np.corrcoef(meanGO1[dd],meanGO2[dd])[1,0],2)
    axs[ii].scatter(meanGO1[dd],meanGO2[dd],c=col, marker='.',s = 3)
    axs[ii].set_title(str(cc1))
    axs[ii].set_xlim([min(meanGO1)-0.2, max(meanGO1)+0.2])
    axs[ii].set_ylim([min(meanGO2)-0.2, max(meanGO2)+0.2])
fig.text(0.5, 0.02, GO1+' ('+str(len(GeneSet1))+')', ha='center', va='center')
fig.text(0.09, 0.5, GO2+' ('+str(len(GeneSet2))+')', ha='center', va='center', rotation='vertical')
plt.show()


# ====================================================================================================
# prostate tumour-grade prediction (multi-layer perceptron, regression + threshold)
# ====================================================================================================
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# ---- settings -------------------------------------------------------------------------------------------
# Paper Methods: 32 cfd-BP pairs (Table S4), Gleason grade (PGG) 4 removed, 60 % train / 40 % test.
# Values used in the notebook cells 157 / 159 (kept for reference): frac = 0.7 / 0.8, cut-off 0.98 / 0.99, PGG 4 kept.
PAIR_CUTOFF   = 0.98      # pairs with sc_pri3 > cut-off and sc_pri5 > cut-off  (0.98 gave 32 pairs in cell 146)
TRAIN_FRAC    = 0.6
EXCLUDE_PGG4  = True
# ---------------------------------------------------------------------------------------------------------

# distance features made by dis_hd_PG.py: column 0 = label (0 healthy, 3 / 4 / 5 = primary Gleason grade), then 3 columns per BP pair
dis_hd = pd.read_csv(Human_Clust_path + 'Distances/dis_hd_mir_Prostate-PG_learn.csv', skiprows = None, sep = ',', header = None)
atmp = pd.read_csv(Human_Clust_path + 'Distances/dis_hd_mir_Prostate-PG_go.csv', skiprows = None, sep = ',', header = 0)
print(dis_hd.shape, atmp.shape)

# select the BP pairs that separate PGG <=3 from PGG 5 (the pairs of Table S4)
dis_hd2 = dis_hd.values[:,0]
tmp = np.where((atmp['sc_pri3']>PAIR_CUTOFF) & (atmp['sc_pri5']>PAIR_CUTOFF))[0]
print(len(tmp), 'BP pairs selected')
atmp.iloc[tmp].to_csv('Prostate_selected_BPpairs.csv', index=False)      # source of Table S4
for ii in tmp:
    dis_hd2 = np.column_stack((dis_hd2,dis_hd.values[:,(ii*3+1):(ii*3+4)]))
if EXCLUDE_PGG4:
    dis_hd2 = dis_hd2[dis_hd2[:,0]!=4,:]

# train / test split
dataset = pd.DataFrame(dis_hd2)
train_dataset = dataset.sample(frac=TRAIN_FRAC, random_state=0)
test_dataset = dataset.drop(train_dataset.index)
train_features = train_dataset.copy()
test_features = test_dataset.copy()
train_labels = train_features.pop(0)
test_labels = test_features.pop(0)
normalizer = tf.keras.layers.Normalization(axis=-1)
normalizer.adapt(np.array(train_features))

col0 = np.where(test_labels==0)[0]      # healthy
col1 = np.where(test_labels==3)[0]      # PGG <= 3
col3 = np.where(test_labels==5)[0]      # PGG 5

# multi-layer perceptron: normalisation layer + 3 ReLU hidden layers, MAE loss, Adam
tf.random.set_seed(0)
def build_and_compile_model(norm):
    model = keras.Sequential([
        norm,
        layers.Dense(10, activation='relu'),
        layers.Dense(25, activation='relu'),
        layers.Dense(25, activation='relu'),
        layers.Dense(1)
    ])
    model.compile(loss='mean_absolute_error',
                  optimizer=tf.keras.optimizers.Adam(0.001))
    return model
dnn_model = build_and_compile_model(normalizer)
history = dnn_model.fit(
    train_features,
    train_labels,
    validation_split=0.2,
    verbose=0, epochs=100)
test_predictions = dnn_model.predict(test_features).flatten()

# Figure 5C: distribution of the predicted grade score; green = healthy, orange = PGG <=3, purple = PGG 5
data = (col0,col1,col3)
cols = ['green','orange','purple']
plt.figure(figsize=(10,5))
for col,tmp in zip(cols,data):
    sns.distplot(test_predictions[tmp], color=col)
    plt.xlim(-0.5,5.5)
plt.axvline(x=4)                 # decision threshold between low and high grade
plt.axvline(x=3.5,ls = '-.')
plt.axvline(x=4.5,ls = '-.')
plt.xlabel('Predicted scores')
plt.show()
