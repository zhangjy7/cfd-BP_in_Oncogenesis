library(dplyr)
library(DT)
library(SummarizedExperiment)
library(mgsa)

Projects = read.table('/Volumes/InMoLab/TCGA data/Analysis/TCGA_projects.tsv', sep = '\t', header = T)

Projects_nms = Projects[,1]
Projects_nms = as.character(Projects_nms)


for (it in 1:11) {
  
  print((it*3-2):(it*3))
  query_TCGA <- GDCquery(project = Projects_nms[(it*3-2):(it*3)],
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         file.type = "htseq.counts.gz"
                         #barcode = listSamples,
                         #legacy = TRUE
  )
  
  GDCdownload(query_TCGA, directory = "/Volumes/InMoLab/TCGA data/Original Data/Cancers/TCGA_proj/GDCdata_counts")
}

for(it in 1:33){
  clinical <- GDCquery_clinic(project = Projects_nms[it], type = "clinical", save.csv = T)
}

for (it in 1:33){
  
  query_TCGA <- GDCquery(project = Projects_nms[it],
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         file.type = "htseq.counts.gz"
                         #barcode = listSamples,
                         #legacy = TRUE
  )
  info = getResults(query_TCGA)
  write.csv(info, paste0(Projects_nms[it], '_counts.csv'), quote = F, row.names = F)
}

# for (it in 1:33){
#   case = read.csv(paste0())
#   case = (getResults(query_TCGA)[,8])
#   CaseNA = (which(case==""))
#   query_lung$results[CaseNA, 8]
#   
#   case_clin = clinical[,67]
#   counts_id = (getResults(query_TCGA)[,8])
#   
#   }


