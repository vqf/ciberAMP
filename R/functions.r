.tumors_noN <- function(){
  result <- c("ACC", "BLCA", "BRCA", "CHOL", "COAD", "CESC", "DLBC", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "KICH", "LAML", "LIHC", "LGG","LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
  return(result)
}

.tumors_N <- function(){
  result <- c("KIRC", "LUSC", "LUAD","BRCA", "PRAD", "THCA", "COAD", "KIRP", "STAD", "LIHC", "HNSC", "ESCA", "KICH", "UCEC", "BLCA", "PAAD", "SARC", "READ", "CHOL", "THYM", "CESC", "GBM")
  return(result)
}

.setRowMatrix <- function(l){
  ncols <- length(l);
  result <- matrix(ncol = ncols, nrow = 0)
  colnames(result) <- l
  return(result)
}


.log2 <- function(df){
  result <- log2(df + 1)
  return(result)
}

.downloadExpression <- function(tumor){
  cohort <- paste("TCGA-", tumor, sep="")
  query <- TCGAbiolinks::GDCquery(project = cohort,
                                  legacy = TRUE,
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  platform = "Illumina HiSeq",
                                  file.type = "results",
                                  experimental.strategy = "RNA-Seq",
                                  sample.type = c("Primary Tumor", "Solid Tissue Normal"))
  TCGAbiolinks::GDCdownload(query)
  tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
  return(tumor.exp)
}

.filterExpression <- function(tumor, sign, object, cor.cut, norm.method, filt.method,
                              filt.qnt.cut, filt.var.func, filt.var.cutoff, filt.eta,
                              filt.FC){
  #The count matrix is filtered out of outliers using AAIC and setting a threshold for correlation minimum.

  dataPrep <- TCGAbiolinks::TCGAanalyze_Preprocessing(object = tumor.exp, cor.cut = cor.cut, filename = paste(tumor, "_AAIC_expression.png", sep=""))

  dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(tabDF = dataPrep,
                                                      geneInfo = TCGAbiolinks::geneInfo,
                                                      method = norm.method)

  dataFilt <- TCGAbiolinks::TCGAanalyze_Filtering(tabDF = dataNorm,
                                                  method = filt.method,
                                                  qnt.cut =  filt.qnt.cut,
                                                  var.func = filt.var.func,
                                                  var.cutoff = filt.var.cutoff,
                                                  eta = filt.eta,
                                                  foldChange = filt.FC)

  dataFilt <- as.data.frame(dataFilt)
  dataFilt <- subset(dataFilt, rownames(dataFilt) %in% sign) #<-----
  return(dataFilt)
}

.getDataDEGs <- function(dataFilt, FDR.DEA, FC){


  samplesNT <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                                   typesample = c("NT"))
  samplesTP <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                                   typesample = c("TP"))
  dataDEGs <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                            mat2 = dataFilt[,samplesTP],
                                            metadata = FALSE,
                                            Cond1type = "Normal",
                                            Cond2type = "Tumor",
                                            fdr.cut = filt.FDR.DEA,
                                            logFC.cut = filt.FC,
                                            method = dea.method)

  dataDEGs$Tumor <- rep(tumor, times = nrow(dataDEGs))
  return(dataDEGs)
}

.getDataDEGs.scnavsdip <- function(dataFilt, FDR.DEA, FC, pat.dip, pat.scna){


  
  dataDEGs.scnavsdip <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[, intersect(colnames(dataFilt), pat.dip)],
                                            mat2 = dataFilt[, intersect(colnames(dataFilt), pat.scna)],
                                            metadata = FALSE,
                                            Cond1type = "Diploid",
                                            Cond2type = "SCNA",
                                            fdr.cut = filt.FDR.DEA,
                                            logFC.cut = filt.FC,
                                            method = dea.method)

  dataDEGs.scnavsdip$Tumor <- rep(tumor, times = nrow(dataDEGs))
  return(dataDEGs)
}