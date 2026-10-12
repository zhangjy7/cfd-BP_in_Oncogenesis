library(TCGAbiolinks) ##https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/index.html
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(mgsa)

##### Parse gene coordination #####
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
View(listFilters(ensembl))
for (chr in 1:1) {
  chrN_genes <- getBM(attributes=c('ensembl_gene_id', 'uniprot_gn_symbol','uniprot_gn_id', 'chromosome_name','start_position','end_position','external_gene_name'), filters ='chromosome_name', values = paste(chr), mart = ensembl)
  gene_len = abs(chrN_genes[,5] - chrN_genes[,6])
  
  gene_name = unique(chrN_genes[,1])
  gene_tmp = matrix(0, nrow = length(gene_name), ncol = 8)
  for (rw in 1:length(gene_name)) {
    tmp = which(chrN_genes[,1]==gene_name[rw])
    tmp1 = max(gene_len[tmp])
    tmp1 = which(gene_len[tmp]==tmp1)[1]
    tmp = tmp[1]
    for (cl in c(1:4)) {
      gene_tmp[rw,cl] = chrN_genes[tmp+tmp1-1, cl]
    }
    gene_tmp[rw,8] = chrN_genes[tmp+tmp1-1, 7]
    gene_tmp[rw,5] = chrN_genes[tmp+tmp1-1, 5] + chr*1000000000
    gene_tmp[rw,6] = chrN_genes[tmp+tmp1-1, 6] + chr*1000000000
    gene_tmp[rw,7] = (chrN_genes[tmp+tmp1-1, 5] + chrN_genes[tmp+tmp1-1, 6])/2 + chr*1000000000
  }
  Gene_coord = gene_tmp
}

for (chr in 2:22) {
  print(chr)
  chrN_genes <- getBM(attributes=c('ensembl_gene_id', 'uniprot_gn_symbol','uniprot_gn_id', 'chromosome_name','start_position','end_position','external_gene_name'), filters ='chromosome_name', values = paste(chr), mart = ensembl)
  gene_len = abs(chrN_genes[,5] - chrN_genes[,6])
  
  gene_name = unique(chrN_genes[,1])
  gene_tmp = matrix(0, nrow = length(gene_name), ncol = 8)
  for (rw in 1:length(gene_name)) {
    tmp = which(chrN_genes[,1]==gene_name[rw])
    tmp1 = max(gene_len[tmp])
    tmp1 = which(gene_len[tmp]==tmp1)[1]
    tmp = tmp[1]
    for (cl in 1:4) {
      gene_tmp[rw,cl] = chrN_genes[tmp+tmp1-1, cl]
    }
    gene_tmp[rw,8] = chrN_genes[tmp+tmp1-1, 7]
    gene_tmp[rw,5] = chrN_genes[tmp+tmp1-1, 5] + chr*1000000000
    gene_tmp[rw,6] = chrN_genes[tmp+tmp1-1, 6] + chr*1000000000
    gene_tmp[rw,7] = (chrN_genes[tmp+tmp1-1, 5] + chrN_genes[tmp+tmp1-1, 6])/2 + chr*1000000000
  }
  Gene_coord = rbind(Gene_coord, gene_tmp)
}

for (chr in c("X")) {
  chrN_genes <- getBM(attributes=c('ensembl_gene_id', 'uniprot_gn_symbol','uniprot_gn_id', 'chromosome_name','start_position','end_position','external_gene_name'), filters ='chromosome_name', values = paste(chr), mart = ensembl)
  gene_len = abs(chrN_genes[,5] - chrN_genes[,6])
  
  gene_name = unique(chrN_genes[,1])
  gene_tmp = matrix(0, nrow = length(gene_name), ncol = 8)
  for (rw in 1:length(gene_name)) {
    tmp = which(chrN_genes[,1]==gene_name[rw])
    tmp1 = max(gene_len[tmp])
    tmp1 = which(gene_len[tmp]==tmp1)[1]
    tmp = tmp[1]
    for (cl in 1:4) {
      gene_tmp[rw,cl] = chrN_genes[tmp+tmp1-1, cl]
    }
    gene_tmp[rw,8] = chrN_genes[tmp+tmp1-1, 7]
    gene_tmp[rw,5] = chrN_genes[tmp+tmp1-1, 5] + 23*1000000000
    gene_tmp[rw,6] = chrN_genes[tmp+tmp1-1, 6] + 23*1000000000
    gene_tmp[rw,7] = (chrN_genes[tmp+tmp1-1, 5] + chrN_genes[tmp+tmp1-1, 6])/2 + 23*1000000000
  }
  Gene_coord = rbind(Gene_coord, gene_tmp)
}
colnames(Gene_coord) = c(colnames(chrN_genes[1:6]), 'mid','external_gene_name')

for (chr in c("Y")) {
  chrN_genes <- getBM(attributes=c('ensembl_gene_id', 'uniprot_gn_symbol','uniprot_gn_id', 'chromosome_name','start_position','end_position','external_gene_name'), filters ='chromosome_name', values = paste(chr), mart = ensembl)
  gene_len = abs(chrN_genes[,5] - chrN_genes[,6])
  
  gene_name = unique(chrN_genes[,1])
  gene_tmp = matrix(0, nrow = length(gene_name), ncol = 8)
  for (rw in 1:length(gene_name)) {
    tmp = which(chrN_genes[,1]==gene_name[rw])
    tmp1 = max(gene_len[tmp])
    tmp1 = which(gene_len[tmp]==tmp1)[1]
    tmp = tmp[1]
    for (cl in 1:4) {
      gene_tmp[rw,cl] = chrN_genes[tmp+tmp1-1, cl]
    }
    gene_tmp[rw,8] = chrN_genes[tmp+tmp1-1, 7]
    gene_tmp[rw,5] = chrN_genes[tmp+tmp1-1, 5] + 24*1000000000
    gene_tmp[rw,6] = chrN_genes[tmp+tmp1-1, 6] + 24*1000000000
    gene_tmp[rw,7] = (chrN_genes[tmp+tmp1-1, 5] + chrN_genes[tmp+tmp1-1, 6])/2 + 24*1000000000
  }
  Gene_coord = rbind(Gene_coord, gene_tmp)
}
colnames(Gene_coord) = c(colnames(chrN_genes[1:6]), 'mid','external_gene_name')
write.csv(Gene_coord, 'Gene_coord.csv', quote = F, col.names = T, row.names = F)

##### Parse gene name #####
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
View(listFilters(ensembl))
chrN_genes <- getBM(attributes=c( 'chromosome_name','ensembl_gene_id', 'uniprot_gn_symbol','entrezgene_id','external_gene_name','rnacentral'), mart = ensembl,filters="", values= "with_uniprotsptrembl")
colnames(chrN_genes)=c('chromosome_name','ensembl_gene_id', 'uniprot_gn_symbol','entrezgene_id','external_gene_name','rnacentral')
write.csv(chrN_genes, 'Gene_info.csv', quote = F, col.names = F, row.names = F)


Gene_info = read.csv('/Volumes/InMoLab/TCGA\ data/Analysis/Gene_info.csv')
Gene_coord1= read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/Gene_coord.csv')
external_gene_name = Gene_coord$ensembl_gene_id
tmp = match(external_gene_name,Gene_info$ensembl_gene_id)
Gene_coord = cbind(Gene_coord,Gene_info$external_gene_name[tmp])
write.csv(Gene_coord, 'Gene_coord.csv', quote = F, col.names = T, row.names = F)

#### Make microRNA list #####
mirRNA = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/mirRNA_mirbase.txt',sep = '\t',header = FALSE,skip = 13)
tmp = which(mirRNA[,2]=='miRNA_primary_transcript')
mirRNA = mirRNA[tmp,]
mirRNA_name = mirRNA[,6]
mirRNA_name2 = rep('NA',length(mirRNA_name))
for (ii in 1:length(mirRNA_name)) {
  mirRNA_name2[ii] = strsplit(mirRNA_name[ii],split = ";")[[1]][3]
}
for (ii in 2:length(mirRNA_name)) {
  tmp = strsplit(mirRNA_name2[ii],split = "mir")[[1]][2]
  tmp = paste(strsplit(tmp, "")[[1]][-1],collapse='')
  mirRNA_name2[ii] = paste0('MIR',tmp)
}
mirRNA_name2 = c('mirRNA_name',mirRNA_name2)
write.csv(mirRNA_name2, 'mirRNA_names.txt', quote = F, col.names = F, row.names = F)

#### Generate GO_human RNA + protein list for 2021 updated #####
Gene_info = read.csv('/Volumes/InMoLab/TCGA\ data/Analysis/Gene_info.csv')

GO_protein = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/goa_human.gaf',sep = '\t',header = T,skip = 13)
GO_protein_BP = GO_protein[which(GO_protein$Qualifier=='involved_in'),]
GO_rna = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/goa_human_rna.gaf',sep = '\t',header = T,skip = 9)
GO_rna_BP = GO_rna[which(GO_rna$Qualifier=='involved_in'),]

# GO_protein = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/goa_human2019.gaf',sep = '\t',header = T,skip = 31)
# GO_protein_BP = GO_protein[which((GO_protein$Aspect=='P')&(GO_protein$Qualifier=='')),]
# GO_rna = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/goa_human_rna2019.gaf',sep = '\t',header = T,skip = 15)
# GO_rna_BP = GO_rna[which((GO_rna$Aspect=='P')&(GO_rna$Qualifier=='')),]

DB_Object_ID = rep('na',nrow(GO_rna_BP))
for (ii in 1:nrow(GO_rna_BP)) {
  DB_Object_ID[ii] = strsplit(GO_rna_BP$DB_Object_ID[ii],split = "_")[[1]][1]
}
GO_rna_BP$DB_Object_ID = DB_Object_ID

  
GO_human_ENSG_P = matrix('na',nrow = nrow(GO_protein_BP)+nrow(GO_rna_BP),ncol=7)
colnames(GO_human_ENSG_P) = c('DB','GeneName','GO','P','GeneDescrib','DB_id','ensmbl')
GO_human_ENSG_P[,1] = c(rep('UniProtKB',nrow(GO_protein_BP)),rep('RNAcentral',nrow(GO_rna_BP)))
GO_human_ENSG_P[,3] = c(GO_protein_BP$GO_ID,GO_rna_BP$GO_ID)
GO_human_ENSG_P[,4] = c(GO_protein_BP$Aspect,GO_rna_BP$Aspect)
GO_human_ENSG_P[,5] = c(GO_protein_BP$DB_Object_Name,GO_rna_BP$DB_Object_Name)
GO_human_ENSG_P[,6] = c(GO_protein_BP$DB_Object_ID,DB_Object_ID)
tmp1 = match(GO_rna_BP$DB_Object_ID,Gene_info$rnacentral)
GO_human_ENSG_P[,2] = c(GO_protein_BP$DB_Object_Symbol,Gene_info$external_gene_name[tmp1])
tmp2 = match(GO_protein_BP$DB_Object_Symbol,Gene_info$uniprot_gn_symbol)
GO_human_ENSG_P[,7] = c(Gene_info$ensembl_gene_id[tmp2],Gene_info$ensembl_gene_id[tmp1])

rna_tmp = as.data.frame(GO_human_ENSG_P[(nrow(GO_protein_BP)+1):(nrow(GO_protein_BP)+nrow(GO_rna_BP)),])
rna_idx = rep('NA',nrow(rna_tmp))
for (ii in 1:length(rna_idx)) {
  tmp = strsplit(rna_tmp$GeneDescrib[ii],split = c(" "))[[1]]
  tmpp = 'rna'
  for (jj in tmp[4:length(tmp)]) {
    tmpp = c(tmpp,tolower(strsplit(jj,split = c("-"))[[1]]))
  }
  if('mir'%in%tmpp){
    tt = which(tmpp=='mir')
    rna_idx[ii] = paste(tmpp[tt:(tt+1)],collapse = "")
  }
}
rna_mir = unique(rna_idx)
for (ii in 2:length(rna_mir)) {
  tmp = which(rna_idx==rna_mir[ii])
  if (length(tmp)>1) {
    tmpp = rna_tmp[tmp,]
    tt = which(!is.na(tmpp$ensmbl))[1]
    tt1 = which(is.na(tmpp$ensmbl))
    rna_tmp$ensmbl[tmp[tt1]] = tmpp$ensmbl[tt]
    rna_tmp$GeneName[tmp[tt1]] = tmpp$GeneName[tt]
  }
}
GO_human_ENSG_P[(nrow(GO_protein_BP)+1):(nrow(GO_protein_BP)+nrow(GO_rna_BP)),] = as.matrix(rna_tmp)

write.table(GO_human_ENSG_P, '/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/GO_human_ENSG_P.txt', quote = F, col.names = T, row.names = F,sep = ';')

GO_human_ENSG_P = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/GO_human_ENSG_P2019.txt',header = T,sep = ';')
GO_human_ENSG_P2 = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/GO_human_UPort_ENSG_P2019.txt',header = F,sep = ';')

GO_all = unique(GO_human_ENSG_P$GO)
a = 0
for (GO in GO_all) {
  tmp1 = GO_protein_BP[which(GO_protein_BP$GO_ID==GO),]
  tmp2 = GO_rna_BP[which(GO_rna_BP$GO_ID==GO),]
  tmp3 = GO_human_ENSG_P2[which(GO_human_ENSG_P2$V3==GO),]
  if (abs(length(unique(tmp1$DB_Object_Symbol))-length(unique(tmp3$V2)))>10) {
    a = a+1
    # print(GO)
    # print(paste('pr:',length(unique(tmp1$DB_Object_Symbol))))
    # print(paste('rna:',length(unique(tmp2$DB_Object_Symbol))))
    # print(paste('old:',length(unique(tmp3$V2))))
  }
}

GO_def = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/GO_Def.txt',header = T,sep = '\t')
GO_name = rep('na',nrow(GO_human_ENSG_P))
for (ii in 1:nrow(GO_human_ENSG_P)) {
  tmp = which(GO_def$GO_ID==GO_human_ENSG_P$GO)
  GO_name = GO_def$GO_Name
}
GO_human_ENSG_P = cbind(GO_name, GO_human_ENSG_P)

## for obtain genes in specific GO #####
GO_human_ENSG_P = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/GO_human_ENSG_P2019.txt',header = T,sep = ';')
target_GO = c("GO:0008544","GO:0034599","GO:0030308","GO:0070374","GO:0007517","GO:2001244","GO:0006986","GO:0045454","GO:0009953","GO:0043434","GO:0043536","GO:1902895","GO:0006879","GO:0042472","GO:0048856","GO:0043588","GO:0070498")
target_GO = target_GO[order(target_GO)]
for (ii in target_GO) {
  tmp = GO_human_ENSG_P[which(GO_human_ENSG_P$GO==ii),]
  tmp = tmp[order(tmp$GeneName),]
  nm =   strsplit(ii,split = ":")[[1]]
  nm = paste(nm,collapse='')
  write.xlsx(tmp, "/Volumes/InMoLab/TCGA\ data/Human_ClustGO/GOpairGroups/Prostate.xlsx", sheetName = nm, col.names = TRUE, row.names = F, append = TRUE)
}

target = c("miR199b","miR30a","miR126","miR1906","miR29a","miR130a","miRlet7b","miR21","miR199a","miRlet7g","miR378a","miR181b","miR31","miR210","miR199a","miR20a","miR27b","miR143","miR132","miR99b")
for (ii in target) {
  tmp = which(tolower(GO_human_ENSG_P$GeneName)==tolower(ii))[1]
  print(GO_human_ENSG_P$ensmbl[tmp])
}

##### Generate clustered TF list #####
Gene_info = read.csv('/Volumes/InMoLab/TCGA\ data/Analysis/Gene_info.csv')
Gene_coord= read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/Gene_coord.csv')
TF_info = read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/TFs.txt',sep = '\t')
rmv = which(TF_info$DbTF=='yes')
TF_info = TF_info[rmv,]
tmp1 = which(TF_info$gene_symbol%in%Gene_info$external_gene_name)
tmp2 = which(TF_info$entrez_human%in%Gene_info$entrezgene_id)
tmp = unique(c(tmp1,tmp2))
Gene_coord2 = Gene_coord[match(TF_info$gene_symbol[tmp], Gene_coord$external_gene_name),]

ii = 1
tt = Gene_coord2$external_gene_name[ii]
chrt = Gene_coord2$chromosome_name[ii]
stst  = (Gene_coord2$start_position-Gene_coord2$start_position[ii])
sted = (Gene_coord2$start_position-Gene_coord2$end_position[ii])
eded  = (Gene_coord2$end_position-Gene_coord2$end_position[ii])
edst = (Gene_coord2$end_position-Gene_coord2$start_position[ii])
tmpp = which((stst<1000000)&(stst>0)|(sted<1000000)&(sted>0)|(eded<1000000)&(eded>0)|(edst<1000000)&(edst>0))
rmv = which(Gene_coord2$external_gene_name[tmpp]!=tt)
tmpp = tmpp[rmv]
TF_clust = cbind(Gene_coord2[rep(ii, length(tmpp)),c(8,4:7)],Gene_coord2[tmpp,c(8,4:7)])

for (ii in 2:nrow(Gene_coord2)) {
  tt = Gene_coord2$external_gene_name[ii]
  chrt = Gene_coord2$chromosome_name[ii]
  stst  = (Gene_coord2$start_position-Gene_coord2$start_position[ii])
  sted = (Gene_coord2$start_position-Gene_coord2$end_position[ii])
  eded  = (Gene_coord2$end_position-Gene_coord2$end_position[ii])
  edst = (Gene_coord2$end_position-Gene_coord2$start_position[ii])
  tmpp = which((stst<1000000)&(stst>0)|(sted<1000000)&(sted>0)|(eded<1000000)&(eded>0)|(edst<1000000)&(edst>0))
  rmv = which(Gene_coord2$external_gene_name[tmpp]!=tt)
  tmpp = tmpp[rmv]
  dbtmp = cbind(Gene_coord2[rep(ii, length(tmpp)),c(8,4:7)],Gene_coord2[tmpp,c(8,4:7)])
  TF_clust = rbind(TF_clust,dbtmp)
}
colnames(TF_clust) = c("uniprot_gn_symbol_A", "chromosome_name_A", "start_position_A", "end_position_A", "mid_A" , "uniprot_gn_symbo_B", "chromosome_name_B", "start_position_B" ,"end_position_B", "mid_B" )
write.csv(TF_clust, 'Clustered_TF.csv', quote = F, col.names = T, row.names = F)

# TF clustered regions
TF_clust=read.csv('/Volumes/InMoLab/TCGA\ data/Original\ Data/GO_annotation/Clustered_TF.csv')
TFs = unique(c(TF_clust$uniprot_gn_symbol_A, TF_clust$uniprot_gn_symbo_B))

TF_clust_single = c('symbol','chr','Start','End','Mid')
TF_clust_double = c('symbol','chr','Start','End','Mid')
TF_clust_triple = c('symbol','chr','Start','End','Mid')
TF_region = rbind(as.matrix(TF_clust[,1:5]),as.matrix(TF_clust[,6:10]))
TF_region = TF_region[match(TFs, TF_region[,1]),]
TF_region = TF_region[order(as.numeric(TF_region[,5]),decreasing = F),]
a = 1
c = 1
for (ii in 1:dim(TF_region)[1]) {
  print(a)
  if(a>=dim(TF_region)[1]){
    break
  }
  TF_tmp= c('symbol','chr','Start','End','Mid')
  TF_tmp = rbind(TF_tmp,TF_region[a,])
  for (jj in (a):(dim(TF_region)[1]-1)) {
    b = 1
    dis_tmp =c((as.numeric(TF_region[(jj+1),3])-as.numeric(TF_region[jj,3])),as.numeric(TF_region[(jj+1),4])-as.numeric(TF_region[jj,4]),as.numeric(TF_region[(jj+1),3])-as.numeric(TF_region[jj,4]),as.numeric(TF_region[(jj+1),4])-as.numeric(TF_region[jj,3]))
    dis_tt = length(which(dis_tmp<1000000))
    if(dis_tt==0) {
      b = b+1
      print(b)
      break
    }
    else {
      TF_tmp = rbind(TF_tmp,TF_region[jj+1,])
    }
  }
  TF_tmp[2:nrow(TF_tmp),3:5] = as.numeric(TF_tmp[2:nrow(TF_tmp),3:5])-round(as.numeric(TF_tmp[2,3])/1000000000)*1000000000
  print(c('jj',jj))
  a = jj+1
  if (dim(TF_tmp)[1]==2){
    print('Error')
    TF_clust_single = rbind(TF_clust_single,TF_tmp[2,])
  }
  else if (dim(TF_tmp)[1]==3){
    TF_clust_double = rbind(TF_clust_double,TF_tmp[2:3,],rep(';',5))
  }
  else if (dim(TF_tmp)[1]==4){
    TF_clust_triple = rbind(TF_clust_triple,TF_tmp[2:4,],rep(';',5))
  }
  else {
    TF_tmp = TF_tmp[-1,]
    # write.xlsx(TF_tmp, "/Volumes/InMoLab/TCGA\ data/TF_clusts.xlsx", sheetName = paste0('TFcl_',c), col.names = TRUE, row.names = F, append = TRUE)
    c = c+1
  }
}
# TF_clust_single = TF_clust_single[-1,]
# write.xlsx(TF_clust_single, "/Volumes/InMoLab/TCGA\ data/TF_clusts.xlsx", sheetName = 'single', col.names = TRUE, row.names = F, append = TRUE)
TF_clust_double = TF_clust_double[-1,]
write.xlsx(TF_clust_double, "/Volumes/InMoLab/TCGA\ data/TF_clusts.xlsx", sheetName = 'double', col.names = TRUE, row.names = F, append = TRUE)
TF_clust_triple = TF_clust_triple[-1,]
write.xlsx(TF_clust_triple, "/Volumes/InMoLab/TCGA\ data/TF_clusts.xlsx", sheetName = 'triple', col.names = TRUE, row.names = F, append = TRUE)


# Convert hg19 to hg38 ####
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Define your input coordinates in hg19
input_coordinates <- data.frame(chr = c("chr1", "chr2"), start = c(1000, 2000), end = c(1500, 2500))

# Convert coordinates to hg38
converted_coordinates <- getLiftOver(ensembl, filters = "chromosome_name", values = input_coordinates$chr,start = input_coordinates$start, end = input_coordinates$end,
                                     to = "GRCh38")

print(converted_coordinates)
