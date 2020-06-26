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
.downloadExpression <- function(tumor, tumors.with.normal) {
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
  }
  else if(tumor == "SKCM") {
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

  }
  else if(tumor == "LAML") {
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

  }
  else if(tumor %in% tumors) {
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

  dataPrep <- TCGAbiolinks::TCGAanalyze_Preprocessing(object = object, cor.cut = cor.cut, filename = paste(tumor, "_AAIC_expression.png", sep=""))

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

# This function compares tumor vs normal expression differences and uses "edgeR::exactTest" function to assign significance to those comparisons: finds differentially expressed genes (DEGs)
.getDataDEGs <- function(tumor, dataFilt, FDR.DEA, FC, dea.method){


  samplesNT <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                                   typesample = c("NT"))
  samplesTP <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                                   typesample = c("TP"))
  dataDEGs <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                            mat2 = dataFilt[,samplesTP],
                                            metadata = FALSE,
                                            Cond1type = "Normal",
                                            Cond2type = "Tumor",
                                            fdr.cut = FDR.DEA,
                                            logFC.cut = FC,
                                            method = dea.method)

  dataDEGs$Tumor <- rep(tumor, times = nrow(dataDEGs))
  dataDEGs$Gene_Symbol <- rownames(dataDEGs)
  dataDEGs <- dataDEGs[,c("Gene_Symbol", "logFC", "logCPM", "PValue", "FDR", "Tumor")]
  return(dataDEGs)
}

# This function downloads or prepares somatic copy number alteration data matrices for further steps.
.getSCNAmatrix <- function(tumor) {
  # If user does not provide a suitable cna.matrix, then we use GDC to download thresholded resulting one for the under analysis tumor.
  gistic <- TCGAbiolinks::getGistic(tumor, type = "thresholded")
  rownames(gistic) <- gistic[,1]
  gistic <- gistic[,4:ncol(gistic)]
  colnames(gistic) <- substr(colnames(gistic), start = 1, stop = 12)
  return(gistic)

}

# This function splits cohort by a gene's copy number in three groups: del, amp and neutro (diploid)
.selectDel <- function(new, cna.thr) {
  if(cna.thr == "Deep") {
    group.del <- new[new[,2] == -2, ]
    return(group.del)
  }else if(cna.thr == "Shallow") {
    group.del <- new[new[,2] == -1, ]
    return(group.del)
  }else if(cna.thr == "Both"){
    group.del <- new[new[,2] <= -1, ]
    return(group.del)
  }
}

.selectAmp <- function(new, cna.thr) {
  if(cna.thr == "Deep") {
    group.amp <- new[new[,2] == 2, ]
    return(group.amp)
  }else if(cna.thr == "Shallow") {
    group.amp <- new[new[,2] == 1, ]
    return(group.amp)
  }else if(cna.thr == "Both"){
    group.amp <- new[new[,2] >= 1, ]
    return(group.amp)
  }
}

.selectNeutro <- function(new, cna.thr){
  group.neutro <- NULL
  if(cna.thr == "Deep") {
    group.neutro <- new[new[,2] == 0, ]
  }else if(cna.thr == "Shallow") {
    group.neutro <- new[new[,2] == 0, ]
  }else if(cna.thr == "Both"){
    group.neutro <- new[new[,2] == 0, ]
  }
  return(group.neutro)
}

.selectDiploid <- function(new, cna.thr) {
  if(cna.thr == "Deep") {
    group.neutro <- new[new[,2] == 0, ]
    return(group.neutro)
  }else if(cna.thr == "Shallow") {
    group.neutro <- new[new[,2] == 0, ]
    return(group.neutro)
  }else if(cna.thr == "Both"){
    group.neutro <- new[new[,2] == 0, ]
    return(group.neutro)
  }
}

# This function sets the minimum number of patients to be analyzed
.setMinPat <- function(new, pat.percentage) {
  min <- (nrow(new) * pat.percentage)/100
  if(min >= 3) {
    return(min)
  }else if(min < 3){
    minimum.patients <- 3
  }
  return(min)
}

# This function returns if analysis will continue by group.x <- group.amp or group.x <- group.del
.setGroupX <- function(group.del, group.amp, group.neutro, minimum.patients) {
  if(isTRUE(nrow(group.del) < minimum.patients)) {

    return(group.amp)

  }else if(isTRUE(nrow(group.amp) < minimum.patients)) {

    return(group.del)

  }else if(isTRUE(nrow(group.neutro) < minimum.patients)) {

    return(group.amp)

  }
}

# This function returns if analysis will continue by group.y <- group.neutro or group.y <- group.del
.setGroupY <- function(group.del, group.amp, group.neutro, minimum.patients) {
  if(isTRUE(nrow(group.del) < minimum.patients)) {

    return(group.neutro)

  }else if(isTRUE(nrow(group.amp) < minimum.patients)) {

    return(group.neutro)

  }else if(isTRUE(nrow(group.neutro) < minimum.patients)) {

    return(group.del)

  }
}

# This function evaluates if this is a AMP vs Diploid or DEL vs Diploi case.
.setCond <- function(group.del, group.amp, group.neutro, minimum.patients) {
  if(isTRUE(nrow(group.del) < minimum.patients)) {

    return(c("Group AMP vs Group DIPLOID"))

  }else if(isTRUE(nrow(group.amp) < minimum.patients)) {

    return(c("Group DEL vs Group DIPLOID"))

  }else if(isTRUE(nrow(group.neutro) < minimum.patients)) {

    return(c("Group AMP vs Group DEL"))

  }
}

# This function returns the gene del/amplified % of in patients
.setSizePat <- function(group.del, group.amp, group.neutro, minimum.patients, del.patients, amp.patients) {
  if(isTRUE(nrow(group.del) < minimum.patients)) {

    return(amp.patients)

  }else if(isTRUE(nrow(group.amp) < minimum.patients)) {

    return(del.patients)

  }else if(isTRUE(nrow(group.neutro) < minimum.patients)) {

    return(amp.patients + del.patients)

  }
}

# This function sets the samples IDs to be returned for those who want to consult exactly which barcode is del/amp
.setPatIDs <- function() {
  if(isTRUE(nrow(group.del) < minimum.patients)) {

    return(paste(rownames(group.amp), collapse = ","))

  }else if(isTRUE(nrow(group.amp) < minimum.patients)) {

    return(paste(rownames(group.del), collapse = ","))

  }else if(isTRUE(nrow(group.neutro) < minimum.patients)) {

    return(paste(c(rownames(group.amp), rownames(group.del)), collapse = ","))

  }
}

# This function compares if the genes is differentially expressed when SCN-altered in tumors.
.getDataDEGs_SCNA <- function(tumor, dataFilt, group.x, group.y, filt.FDR.DEA, filt.FC) {
  samplesNT <- rownames(group.y)
  samplesTP <- rownames(group.x)
  dataDEGs <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                            mat2 = dataFilt[,samplesTP],
                                            metadata = FALSE,
                                            Cond1type = "Diploid",
                                            Cond2type = "SCN-altered",
                                            fdr.cut = filt.FDR.DEA,
                                            logFC.cut = filt.FC,
                                            method = "exactTest")

  dataDEGs$Tumor <- rep(tumor, times = nrow(dataDEGs))
  return(dataDEGs)

}

#This function configures the new line for significant DEGs in SCN-altered samples
.newSCNAline <- function(dataDEGs.SCNA, cond, SCNA.prop.pat, pat.ids) {
  returns(c(rownames(dataDEGs.SCNA), dataDEGs.SCNA[1,2], dataDEGs.SCNA[1,3], dataDEGs.SCNA[1,4], dataDEGs.SCNA[1,5], cond, SCNA.prop.pat, pat.ids))
}

# This function converts the SCNA DE resulting matrix into df
.convertToDF <- function(SCNA.DEG.result) {
  SCNA.DEG.result <- as.data.frame(SCNA.DEG.result)
  SCNA.DEG.result$log2FC.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$log2FC.SCNAvsDip))
  SCNA.DEG.result$logCPM.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$logCPM.SCNAvsDip))
  SCNA.DEG.result$p.val.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$p.val.SCNAvsDip))
  SCNA.DEG.result$FDR.SCNAvsDip <- p.adjust(as.numeric(as.character(SCNA.DEG.result$p.val.SCNAvsDip)), method = "fdr")
  SCNA.DEG.result$Pat.percentage <- as.numeric(as.character(SCNA.DEG.result$Pat.percentage))
}

