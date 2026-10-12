# This is a python3 file
print('starting')
import sys
import pandas as pd
import numpy as np
import os
import re
import json
from random import sample
import time 

print('starting')
ii = sys.argv[1]
ii = int(ii)
print(ii)

GO_path = 'GO_annotation/'
GO_BP = pd.read_csv(GO_path + 'GO_human_ENSG_P.txt', skiprows = None, sep = ';', header = 0)
Gene_coord = pd.read_csv(GO_path + 'Gene_coord.csv', skiprows = None, sep = ',', header = 0)
Gene_coord2 = np.where(Gene_coord['external_gene_name'].notnull())
Gene_coord2 = Gene_coord.iloc[Gene_coord2]
GO_Def = pd.read_csv(GO_path + 'GO_Def.txt', skiprows = None, sep = '\t', header = 0)

Gene_All = list(set(Gene_coord2['ensembl_gene_id']) & set(GO_BP['ensmbl']))
print(len(Gene_All))
with open(GO_path + 'BPtoGene2.json') as json_file:
    BP_to_Gene = json.load(json_file)
    BP_to_Gene = BP_to_Gene['BP_to_Gene']
    
with open(GO_path + 'BPtoTF2.json') as json_file:
    BP_to_TF = json.load(json_file)
    BP_to_TF = BP_to_TF['BP_to_TF']

NumGene = []
for i in BP_to_Gene:
    NumGene.append(len(i['Genes']))
#Sel_BP = np.where(np.array(NumGene)>25)[0] #len(Sel_BP) = 929
Sel_BP = np.where((np.array(NumGene)>25) & (np.array(NumGene)<=250))[0] #len(Sel_BP) = 816
print(len(Sel_BP))

NumTF = []
for i in BP_to_TF:
    NumTF.append(len(i['TF']))
Sel_TF = np.where(np.array(NumTF)>0)[0]#len(Sel_TF) = 3330

#len(set(Sel_BP)& set(Sel_TF)) = 315

print('All required files are uploaded')
def Gene_clust(coord1, coord2, dist_thresh):
    ## coord1 and coord2 are two sets of coordinator with every of them is a n*3 matrix: n is gene number, 3 is start, end and middle coordinator of a gene.
    
    clustGeneA = []
    clustGeneB = []
    clustArea = 0
    clust_A = []
    for i in range(coord1.shape[0]):
        dis_tmp = abs(coord2 -coord1[i,:])
        x = np.where(dis_tmp<dist_thresh)[0]
        x = list(set(x))
        if len(x)>0:
            clust_A.append(x)
            clustGeneA.append(i)
            clustGeneB = clustGeneB + x
        if len(x)>1:
            clustArea = clustArea + len(x)-1
    PairGene = len(clustGeneB)
    ClustA = len(clustGeneA)
    ClustB = list(set(clustGeneB))
    
    tmp_index = []
    aa = 0
    for tmp_A in clust_A:
        for xx in tmp_A:
            lgn2 = len(clust_A)
            for yy in range(aa+1,lgn2):
                zz = clust_A[yy]                    
                if xx in zz:
                    tmp_index.append(clust_A[yy:lgn2].index(zz) + yy)
        aa = aa+1
    
    clustArea = clustArea + len(set(tmp_index))            

    clustArea = PairGene - clustArea
    
    return(PairGene, clustArea, ClustA, len(ClustB), clustGeneA, ClustB)


def RandomSampling(Gene_sel, Coord_sel, AllGeneSet, number, iteration):
    # Gene_sel: The selected set of genes, list
    # Coord_sel: The coordinators of selected genes, np array
    # AllGeneSet: All gene pool for analysis, list
    # number of gene generated, int
    # iteration: iteration times, int
    AllGeneSet2 = list(AllGeneSet)
    for i in Gene_sel: AllGeneSet2.remove(i)
#     print(len(AllGeneSet))
    
    tmp_matrix = np.zeros((iteration,3))
    for ilt in range(iteration):
        Gene_tmp = sample(AllGeneSet2, number)
        #idx_tmp = [i for i,val in enumerate(Gene_coord2['ensembl_gene_id']) if val in Gene_tmp]
        idx_tmp = [list(Gene_coord2['ensembl_gene_id']).index(i) for i in Gene_tmp]
        coord_tmp = Gene_coord2.values[idx_tmp,4:7]
        [PairGene, clustArea, ClustA, ClustB, clustGeneA, clustGeneB] = Gene_clust(Coord_sel, coord_tmp, 1000000)
        tmp_matrix[ilt, :] = [PairGene, clustArea, ClustA]
        
    return list(tmp_matrix.mean(axis = 0))


iteration = 1000

for iBP in range((ii-1)*10,(ii*10)):
#for iBP in [1]:
    print(iBP)
    #if os.path.exists('Human_Clust/Human_Clust'+ str(iBP) +'_b2022.csv'):
    #    print(iBP, 'exsited')
    #    continue


    Human_Clust = pd.DataFrame(columns=['BP_indexA', 'BP_nameA','BP_indexB', 'BP_nameB', 'Sum_Genes', 'Pos_Gpairs', 'Sum_TFs', 'Pos_Tpairs', 'Pair_Gene', 'AreaGene', 'ClustGeneA', 'ClustGeneB', 'Pair_GeneR', 'AreaGeneR', 'ClustGeneAR', 'ClustGeneBR', 'Pair_TF', 'AreaTF', 'ClustTFA', 'ClustTFB', 'Pair_TFR', 'AreaTFR', 'ClustTFAR', 'ClustTFBR'])

    Clust_list = {}
    Clust_list['Clustered_GO'] = []
    Clust_list['Clustered_Gene'] = []
    Clust_list['Clustered_TF'] = []

    start_time = time.time()
    BP_indexA = Sel_BP[iBP] #the index number of BP term in Sel_BP
    
    ## Record the index and name of GO_termA to both csv and json file ##
    GO_idxA = BP_to_Gene[BP_indexA]['GO_index']
    GO_nameA = BP_to_Gene[BP_indexA]['GO_Name']
    
    geneSetA = BP_to_Gene[BP_indexA]['Genes'] ## Genes included into GO termA
    coord1 = np.transpose(np.array((BP_to_Gene[BP_indexA]['Gene_st'], 
                                     BP_to_Gene[BP_indexA]['Gene_ed'], 
                                     BP_to_Gene[BP_indexA]['Gene_mid'])))## Coordinator of genes included into GO termA
    print(iBP, len(geneSetA), coord1.shape)
    
    TFSetA = BP_to_TF[BP_indexA]['TF']## TFs included into GO termA
    coordTF1 = np.transpose(np.array((BP_to_TF[BP_indexA]['TF_st'], ## TFs included into GO termA
                                       BP_to_TF[BP_indexA]['TF_ed'], 
                                       BP_to_TF[BP_indexA]['TF_mid'])))## Coordinator of TFs included into GO termA
        
    for jBP in range(iBP+1, 5):
    #for jBP in range(iBP+1, len(Sel_BP)):
        print(iBP, len(Sel_BP))
        Human_Clust_tmp = [GO_idxA, GO_nameA]
        BP_indexB = Sel_BP[jBP] #the index number of BP term in Sel_BP
        
        ## Record the index and name of GO_termB to both csv and json file ##
        GO_idxB = BP_to_Gene[BP_indexB]['GO_index']
        GO_nameB = BP_to_Gene[BP_indexB]['GO_Name']
        Human_Clust_tmp.append(GO_idxB)
        Human_Clust_tmp.append(GO_nameB)
        
        
        geneSetB = BP_to_Gene[BP_indexB]['Genes']  ## Genes included into GO termB
        coord2 = np.transpose(np.array((BP_to_Gene[BP_indexB]['Gene_st'], 
                                         BP_to_Gene[BP_indexB]['Gene_ed'], 
                                         BP_to_Gene[BP_indexB]['Gene_mid'])))## Coordinator of genes included into GO termB
        
        TFSetB = BP_to_TF[BP_indexB]['TF']## TFs included into GO termB
        coordTF2 = np.transpose(np.array((BP_to_TF[BP_indexB]['TF_st'], 
                                           BP_to_TF[BP_indexB]['TF_ed'], 
                                           BP_to_TF[BP_indexB]['TF_mid'])))## Coordinator of TFs included into GO termA
    
        
        ## Check and remove the overlapped genes from both BP terms ##
        olGene = list(set(geneSetA) & set(geneSetB))
        numGA = len(list(set(geneSetA)))
        numGB = len(list(set(geneSetB)))
        if len(olGene)> numGA*(1/3) or len(olGene)> numGB*(1/3):
            continue # if the overlap between GO term-A and GO term-B over 1/3, then the results are not recorded, continue directly to the next GO term-B.
        
        Clust_list['Clustered_GO'].append({
            'clust_GOA': GO_idxA,
            'clust_GOB': GO_idxB,
            'clust_NameA': GO_nameA,
            'clust_NameB': GO_nameB
        })
        
        Human_Clust_tmp.append(numGA+numGB)
        Human_Clust_tmp.append(numGA*numGB)
        
        numTA = len(list(set(TFSetA)))
        numTB = len(list(set(TFSetB)))
        
        Human_Clust_tmp.append(numTA+numTB)
        Human_Clust_tmp.append(numTA*numTB)
        
        ### remove the overlaid genes 
        rmvA = []
        rmvB = []        
        if len(olGene) >0:
            for k in range(len(olGene)):
                rmvA.append(geneSetA.index(olGene[k]))
                rmvB.append(geneSetB.index(olGene[k]))
            
        coord11 = np.delete(coord1, rmvA, 0)
        coord22 = np.delete(coord2, rmvB, 0)        
        [PairGene, clustArea, ClustA, ClustB, clustGeneA, clustGeneB] = Gene_clust(coord11, coord22, 1000000)
        print(jBP, len(olGene), len(rmvA), coord1.shape, coord11.shape, coord2.shape, coord22.shape)
        geneSetAA = list(set(range(len(geneSetA))) - set(rmvA))
        geneSetAA = list(np.array(geneSetA)[geneSetAA])
        geneSetBB= list(set(range(len(geneSetB))) - set(rmvB))
        geneSetBB = list(np.array(geneSetB)[geneSetBB])
            
        Human_Clust_tmp = Human_Clust_tmp + [PairGene, clustArea, ClustA, ClustB]
        
        Clust_list['Clustered_Gene'].append({
            'clust_GeneA': [geneSetAA[xx] for xx in clustGeneA],
            'clust_GeneB': [geneSetBB[xx] for xx in clustGeneB]
        })
        
        ### Random gene check for BP pair A and B
        # Gene_All = list(set(Gene_coord2['ensembl_gene_id']) & set(GO_BP['ensmbl']))
        rndA = RandomSampling(geneSetA, coord1, Gene_All, len(geneSetB), iteration)
        print(len(Gene_All))
        # Gene_All = list(set(Gene_coord2['ensembl_gene_id']) & set(GO_BP['ensmbl']))
        rndB = RandomSampling(geneSetB, coord2, Gene_All, len(geneSetA), iteration)
        print(len(Gene_All))
        Human_Clust_tmp = Human_Clust_tmp + [max(rndA[0], rndB[0]), max(rndA[1], rndB[1]), rndA[2], rndB[2]]
        ## Clustering situation of TF from two BPs
        if numTA*numTB==0:
            Human_Clust_tmp = Human_Clust_tmp + [0, 0, 0, 0, 0, 0, 0, 0]
            Clust_list['Clustered_TF'].append({
                'clust_GeneA': [],
                'clust_GeneB': []
            })
            # print(Human_Clust_tmp)
            Human_Clust.loc[len(Human_Clust)] = Human_Clust_tmp
        else:
            olGene = list(set(TFSetA) & set(TFSetB))
            rmvA = []
            rmvB = []
            if len(olGene) >0:
                for k in range(len(olGene)):
                    rmvA.append(TFSetA.index(olGene[k]))
                    rmvB.append(TFSetB.index(olGene[k]))

            coordTF11 = np.delete(coordTF1, rmvA, 0)
            coordTF22 = np.delete(coordTF2, rmvB, 0)
            [PairTF, TFclustArea, ClustTA, ClustTB, clustTFA, clustTFB] = Gene_clust(coordTF11, coordTF22, 1000000)
            
            Human_Clust_tmp = Human_Clust_tmp + [PairTF, TFclustArea, ClustTA, ClustTB]
            
            geneSetAA = list(set(range(len(TFSetA))) - set(rmvA))
            geneSetAA = list(np.array(TFSetA)[geneSetAA])
            geneSetBB= list(set(range(len(TFSetB))) - set(rmvB))
            geneSetBB = list(np.array(TFSetB)[geneSetBB])
            

            Clust_list['Clustered_TF'].append({
            'clust_GeneA': [geneSetAA[xx] for xx in clustTFA],
            'clust_GeneB': [geneSetBB[xx] for xx in clustTFB]
            })
            
            ### Random TF check for BP pair A and B
            # Gene_All = list(set(Gene_coord2['ensembl_gene_id']) & set(GO_BP['ensmbl']))
            rndA = RandomSampling(TFSetA, coordTF1, Gene_All, len(TFSetB), iteration)
            # Gene_All = list(set(Gene_coord2['ensembl_gene_id']) & set(GO_BP['ensmbl']))
            rndB = RandomSampling(TFSetB, coordTF2, Gene_All, len(TFSetA), iteration)
            Human_Clust_tmp = Human_Clust_tmp + [max(rndA[0], rndB[0]), max(rndA[1], rndB[1]), rndA[2], rndB[2]]
            Human_Clust.loc[len(Human_Clust)] = Human_Clust_tmp

    print('Saving document: '+'Human_Clust/Human_Clust'+ str(iBP) +'_b.csv')
    Human_Clust.to_csv('Human_ClustGO/Human_Clust'+ str(iBP) +'_b.csv', index=False)

    print('Saving document: '+'Human_Clust/Human_Clust' + str(iBP) + '_b.json')
    with open('Human_ClustGO/Human_Clust' + str(iBP) + '_b.json', 'w') as outfile:
        json.dump(Clust_list, outfile)
    print("--- %s seconds ---" % (time.time() - start_time))

