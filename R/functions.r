.tumors_noN <- function(){
# This function creates a complete list of TCGA cohorts
.tumors_all <- function(){
  result <- c("ACC", "BLCA", "BRCA", "CHOL", "COAD", "CESC", "DLBC", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "KICH", "LAML", "LIHC", "LGG","LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
  return(result)
}

# This function creates a complete list of TCGA cohorts with normal samples as control
.tumors_N <- function(){
  result <- c("KIRC", "LUSC", "LUAD","BRCA", "PRAD", "THCA", "COAD", "KIRP", "STAD", "LIHC", "HNSC", "ESCA", "KICH", "UCEC", "BLCA", "PAAD", "SARC", "READ", "CHOL", "THYM", "CESC", "GBM")
  return(result)
}

# This function creates matrices with 0 rows and the specified number (l) of columns
.setRowMatrix <- function(l){
  ncols <- length(l);
  result <- matrix(ncol = ncols, nrow = 0)
  colnames(result) <- l
  return(result)
}

# This function downloads tumor expression from GDC and includes it in a SummarizedExperiment object knwon as "tumor.exp"
.downloadExpression <- function(tumor) {
  if(tumor %in% tumors.with.normal) {
    # If the required cohort has normal samples, then we have to indicate it in the TCGAbiolinks::GDCquery function argument "sample.type".
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
  }else if(tumor == c("SKCM") {
    # If tumor is SKCM cohort, then "sample.type" argument should be specified for only primary tumors --> in this cohort there are many metastatic samples which are going to be analyzed  <-- SHOULD WE GIVE THE PEOPLE OPTION TO exCLUDE THEM???
    cohort <- paste("TCGA-", tumor, sep="")
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c("Primary solid Tumor", "Metastatic"))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)

  }else if(tumor == c("LAML")) {
    # If tumor is SKCM cohort, then "sample.type" argument should be specified since this is a blood type malignancy.
    cohort <- paste("TCGA-", tumor, sep="")
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c("Primary Blood Derived Cancer - Peripheral Blood"))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)

  }else if(tumor %in% tumors) {
    # In any other case, if the user's "cohorts" input is valid (exists in TCGA), then this will be a non-normal sample cohort and this is specified in "sample.type" argument as well.
    cohort <- paste("TCGA-", tumor, sep="")
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c("Primary solid Tumor"))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)

  }
}

# This function filters lower expressed genes and normalizes it by "EDASeq::withinLaneNormalization" and "EDASeq::betweenLaneNormalization" functions (both built in TCGAbiolinks::TCGAanalyze_Normalization R package function)
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



  dataDEGs.scnavsdip <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[gene, intersect(colnames(dataFilt), pat.dip), drop = FALSE],
                                            mat2 = dataFilt[gene, intersect(colnames(dataFilt), pat.scna), drop = FALSE],
                                            metadata = FALSE,
                                            Cond1type = "Diploid",
                                            Cond2type = "SCNA",
                                            fdr.cut = filt.FDR.DEA,
                                            logFC.cut = filt.FC,
                                            method = dea.method)

  dataDEGs.scnavsdip$Tumor <- rep(tumor, times = nrow(dataDEGs))
  return(dataDEGs)
}
