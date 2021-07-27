realTP <- NULL
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
.setRowMatrix <- function(nrow, l){
  ncols <- length(l);
  result <- matrix(nrow = nrow, ncol = ncols)
  colnames(result) <- l
  return(result)
}

#' Crude hack to get Primary Tumor name
#'
#' @return Primary Tumor name
.primaryTumorName <- function(){
  result <- "Primary solid Tumor"
  if (is.null(realTP)){
    sink(file = tempfile())
    result <- tryCatch({
      suppressMessages(suppressWarnings({query <- TCGAbiolinks::GDCquery(project = "TCGA-ACC",
                                                                         data.category =  "Copy Number Variation",
                                                                         data.type = "Masked Copy Number Segment",
                                                                         sample.type = c(result))}))
      sink()
      realTP <- result
      return(result)
    }, error=function(e){
      if (grepl(pattern = 'sample.type', x = e)){
        return("Primary Tumor")
      }
    })
    sink()
    realTP <- result
    return(result)
  }
  else{
    return(realTP)
  }
}

# This function downloads tumor expression from GDC and includes it in a SummarizedExperiment object knwon as "tumor.exp"
.downloadExpression <- function(tumor, tumors, tumors.with.normal) {
  cohort <- paste("TCGA-", tumor, sep="")
  tp <- .primaryTumorName()
  if(tumor %in% tumors.with.normal) {
    # If the required cohort has normal samples, then we have to indicate it in the TCGAbiolinks::GDCquery function argument "sample.type".
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c(tp, "Solid Tissue Normal"))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)
  }
  else if(tumor == "SKCM") {
    # If tumor is SKCM cohort, then "sample.type" argument should be specified for only primary tumors --> in this cohort there are many metastatic samples which are going to be analyzed  <-- SHOULD WE GIVE THE PEOPLE OPTION TO exCLUDE THEM???
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c(tp, "Metastatic"))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)

  }
  else if(tumor == "LAML") {
    # If tumor is SKCM cohort, then "sample.type" argument should be specified since this is a blood type malignancy.
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
    query <- TCGAbiolinks::GDCquery(project = cohort,
                                    legacy = TRUE,
                                    data.category = "Gene expression",
                                    data.type = "Gene expression quantification",
                                    platform = "Illumina HiSeq",
                                    file.type = "results",
                                    experimental.strategy = "RNA-Seq",
                                    sample.type = c(tp))
    TCGAbiolinks::GDCdownload(query)

    tumor.exp <- TCGAbiolinks::GDCprepare(query = query)
    return(tumor.exp)

  }
}

# This function filters lower expressed genes and normalizes it by "EDASeq::withinLaneNormalization" and "EDASeq::betweenLaneNormalization" functions (both built in TCGAbiolinks::TCGAanalyze_Normalization R package function)
.filterExpression <- function(tumor, sign, object, pp.cor.cut, norm.method, filt.method,
                              filt.qnt.cut, filt.var.func, filt.var.cutoff, filt.eta,
                              filt.FC){

  dataPrep <- TCGAbiolinks::TCGAanalyze_Preprocessing(object = object, cor.cut = pp.cor.cut, filename = paste(tumor, "_AAIC_expression.png", sep=""))

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
.getDataDEGs <- function(tumor, dataFilt, filt.FDR.DEA, FC){

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
                                            logFC.cut = FC,
                                            method = "exactTest")

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
  if(min >= 2) {
    return(min)
  }else if(min < 2){
    return(2)
    }
   }

# This function compares if the genes is differentially expressed when SCN-altered in tumors.
.getDataDEGs_SCNA <- function(tumor, dataFilt, group.x, group.y, filt.FDR.DEA, filt.FC, gene) {
  samplesNT <- rownames(group.y)
  samplesTP <- rownames(group.x)
  dataDEGs <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = dataFilt[, samplesNT],
                                            mat2 = dataFilt[,samplesTP],
                                            metadata = FALSE,
                                            Cond1type = "Diploid",
                                            Cond2type = "SCN-altered",
                                            fdr.cut = filt.FDR.DEA,
                                            logFC.cut = filt.FC,
                                            method = "exactTest")

  dataDEGs$Tumor <- rep(tumor, times = nrow(dataDEGs))

  if(gene %in% rownames(dataDEGs)) {
    dataDEGs <- dataDEGs[gene, ]
    return(dataDEGs)
  }else{
    return(NULL)
  }
}

#This function configures the new line for significant DEGs in SCN-altered samples.
.newSCNAline <- function(dataDEGs.SCNA, cond, SCNA.prop.pat, pat.ids) {
  return(c(rownames(dataDEGs.SCNA), dataDEGs.SCNA[1,1], dataDEGs.SCNA[1,2], dataDEGs.SCNA[1,3], dataDEGs.SCNA[1,4], dataDEGs.SCNA[1,5], cond, SCNA.prop.pat, pat.ids))
}

# This function converts the SCNA DE resulting matrix into df
.convertToDF <- function(SCNA.DEG.result) {
  SCNA.DEG.result <- as.data.frame(SCNA.DEG.result)
  SCNA.DEG.result$log2FC.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$log2FC.SCNAvsDip))
  SCNA.DEG.result$logCPM.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$logCPM.SCNAvsDip))
  SCNA.DEG.result$p.val.SCNAvsDip <- as.numeric(as.character(SCNA.DEG.result$p.val.SCNAvsDip))
  SCNA.DEG.result$FDR.SCNAvsDip <- p.adjust(as.numeric(as.character(SCNA.DEG.result$p.val.SCNAvsDip)), method = "fdr")
  SCNA.DEG.result$Pat.percentage <- as.numeric(as.character(SCNA.DEG.result$Pat.percentage))
  return(SCNA.DEG.result)
}

# This function combines dataDEGs and SCNA.DEG.result: 1) evaluated if dataDEGs is null or not and 2) creates a new integrated matrix with both data (if dataDEGs is null -> NAs introduced)
.mergeDEGs <- function(dataDEGs, SCNA.DEG.result, pat.percentage, genes, cosmic.genes) {
  if(!is.null(dataDEGs) && nrow(dataDEGs[dataDEGs$Gene_Symbol %in% genes, ]) > 0) {
    s <-  merge(dataDEGs, SCNA.DEG.result, by = "Gene_Symbol", all = TRUE)
    s$Condition <- as.character(s$Condition)
    s$Pat.IDs <- as.character(s$Pat.IDs)
    s$Tumor <- as.character(s$Tumor)
    s$TCGA_Tumor <- as.character(s$TCGA_Tumor)

    for(i in 1:nrow(s)) {
      if(is.na(s$logFC[i])) {
        s$logFC[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$logCPM[i])) {
        s$logCPM[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$PValue[i])) {
        s$PValue[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$FDR[i])) {
        s$FDR[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Tumor[i])) {
        s$Tumor[i] <- as.character(s$TCGA_Tumor[i])
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$log2FC.SCNAvsDip[i])) {
        s$log2FC.SCNAvsDip[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$logCPM.SCNAvsDip[i])) {
        s$logCPM.SCNAvsDip[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$p.val.SCNAvsDip[i])) {
        s$p.val.SCNAvsDip[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$FDR.SCNAvsDip[i])) {
        s$FDR.SCNAvsDip[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$TCGA_Tumor[i])) {
        s$TCGA_Tumor[i] <- as.character(s$Tumor[i])
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Condition[i])) {
        s$Condition[i] <- c("Not SCN-associated DEG")
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Pat.percentage[i])) {
        s$Pat.percentage[i] <- 0
      }
    }

    return(s)

  }else if(is.null(dataDEGs) || nrow(dataDEGs[dataDEGs$Gene_Symbol %in% genes, ]) == 0) {

    d <- as.data.frame(.setRowMatrix(nrow(SCNA.DEG.result), c("Gene_Symbol", "logFC", "logCPM", "PValue", "FDR", "Tumor")))
    d$Gene_Symbol <- SCNA.DEG.result$Gene_Symbol
    s <-  merge(d, SCNA.DEG.result, by = "Gene_Symbol", all = TRUE)
    s$Condition <- as.character(s$Condition)
    s$Pat.IDs <- as.character(s$Pat.IDs)
    s$Tumor <- as.character(s$Tumor)
    s$TCGA_Tumor <- as.character(s$TCGA_Tumor)

    for(i in 1:nrow(s)) {
      if(is.na(s$logFC[i])) {
        s$logFC[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$logCPM[i])) {
        s$logCPM[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$PValue[i])) {
        s$PValue[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$FDR[i])) {
        s$FDR[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Tumor[i])) {
        s$Tumor[i] <- as.character(s$TCGA_Tumor[i])
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$log2FC.SCNAvsDip[i])) {
        s$log2FC.SCNAvsDip[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$logCPM.SCNAvsDip[i])) {
        s$logCPM.SCNAvsDip[i] <- 0
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$p.val.SCNAvsDip[i])) {
        s$p.val.SCNAvsDip[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$FDR.SCNAvsDip[i])) {
        s$FDR.SCNAvsDip[i] <- 1
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$TCGA_Tumor[i])) {
        s$TCGA_Tumor[i] <- as.character(s$Tumor[i])
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Condition[i])) {
        s$Condition[i] <- paste("Not SCN-altered in more than", pat.percentage, " ")
      }
    }

    for(i in 1:nrow(s)) {
      if(is.na(s$Pat.percentage[i])) {
        s$Pat.percentage[i] <- 0
      }
    }

    return(s)
  }
}

# This function returns for each copy number altered gene all those cosmic ones highly correlating: 1) CN-altered samples overlap > 70% and 2) DE correlation

.getSigPairs <- function(SCNA.DEG.result, cna, cna.thr, tumor) {
  scn.degs <- as.character(SCNA.DEG.result[SCNA.DEG.result$log2FC.SCNAvsDip != 0, ]$Gene_Symbol)
  mutMat <- cna[,scn.degs, drop = FALSE]
  mutMat <- as.matrix(mutMat)
  if(cna.thr == "Deep") {
    mutMat[mutMat == '-1'] <- '0'
    mutMat[mutMat == '1'] <- '0'
    mutMat[mutMat != '0'] <- '1'
  }else if(cna.thr == "Shallow") {
    mutMat[mutMat == '-2'] <- '0'
    mutMat[mutMat == '2'] <- '0'
    mutMat[mutMat != '0'] <- '1'
  }else{
    mutMat[mutMat != '0'] <- '1'
  }
  mutMat <- as.matrix(mutMat)
  print("Calculating the interactions...")
  interactions = sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat),
                                                           function(j) {
                                                             print(paste(colnames(mutMat)[i], colnames(mutMat)[j], sep="_"))
                                                             f <- try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
                                                             if (class(f) == "try-error")
                                                               NA
                                                             else ifelse(f$estimate > 1, -log10(f$p.val), log10(f$p.val))
                                                           }))
  print("Calculating the odds ratio...")
  oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat),
                                                                      function(j) {
                                                                        print(paste(colnames(mutMat)[i], colnames(mutMat)[j], sep="_"))
                                                                        f <- try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
                                                                        if (class(f) == "try-error")
                                                                          f = NA
                                                                        else f$estimate
                                                                      }))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  sigPairs = which(x = 10^-abs(interactions) < max(0.05), arr.ind = TRUE)
  sigPairsTbl = data.table::rbindlist(lapply(X = seq_along(1:nrow(sigPairs)),
                                             function(i) {
                                               x = sigPairs[i, ]
                                               g1 = rownames(interactions[x[1], x[2], drop = FALSE])
                                               g2 = colnames(interactions[x[1], x[2], drop = FALSE])
                                               tbl = as.data.frame(table(apply(X = mutMat[, c(g1,g2), drop = FALSE], 1, paste, collapse = "")))
                                               combn = data.frame(t(tbl$Freq))
                                               colnames(combn) = tbl$Var1
                                               pval = 10^-abs(interactions[x[1], x[2]])
                                               fest = oddsRatio[x[1], x[2]]
                                               d = data.table::data.table(gene1 = g1, gene2 = g2,pValue = pval, oddsRatio = fest)
                                               d = cbind(d, combn)
                                               d
                                             }), fill = TRUE)

  sigPairsTbl <- sigPairsTbl[sigPairsTbl$gene1 != sigPairsTbl$gene2, ]
  sigPairsTbl <- sigPairsTbl[sigPairsTbl$gene2 %in% all_cosmic_genes(), ]
  sigPairsTbl$Event <- ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurence", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[, c('gene1', 'gene2')], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTblSig <- sigPairsTbl[order(as.numeric(sigPairsTbl$pValue)) & !duplicated(sigPairsTbl$pair), ]
  sigPairs.cosmic <- sigPairsTblSig[sigPairsTblSig$gene1 %in% all_cosmic_genes() | sigPairsTblSig$gene2 %in% all_cosmic_genes(), ]
  sigPairs.cosmic$TCGA_Tumor <- rep(tumor, nrow(sigPairs.cosmic))
  colnames(sigPairs.cosmic)[1] <- "Gene_Symbol"
  colnames(sigPairs.cosmic)[2] <- "Gene_Symbol_COSMIC"
  return(sigPairs.cosmic)
}

# This function returns for each copy number altered gene all those cosmic ones highly correlating: 1) CN-altered samples overlap > 70% and 2) DE correlation
.getOverlapCOSMIC <- function(SCNA.DEG.result, genes, cosmic.genes) {
  int.matrix <- .setRowMatrix(0, colnames(SCNA.DEG.result))
  input <- SCNA.DEG.result[SCNA.DEG.result$Gene_Symbol %in% genes & SCNA.DEG.result$log2FC.SCNAvsDip != 0, ]
  cosmic <- SCNA.DEG.result[SCNA.DEG.result$Gene_Symbol %in% cosmic.genes & SCNA.DEG.result$log2FC.SCNAvsDip != 0, ]
  colnames(cosmic) <- paste(colnames(cosmic), "COSMIC", sep="_")

  if(nrow(input) > 0 & nrow(cosmic) > 0) {
    for(j in 1:nrow(input)) {
      for(k in 1:nrow(cosmic)) {
        if(input$TCGA_Tumor[j] == cosmic$TCGA_Tumor[k] & as.character(input$Condition[j]) == as.character(cosmic$Condition[k])) {
          a <- unlist(strsplit(as.character(input$Pat.IDs[j]), ","))
          b <- unlist(strsplit(as.character(cosmic$Pat.IDs[k]), ","))
          int <- intersect(a, b)
          a.l <- length(a)
          b.l <- length(b)
          int.l <- length(int)
          prop.gene.cosmic <- (int.l/a.l) * 100
          prop.cosmic.gene <- (int.l/b.l) * 100
          if(prop.gene.cosmic >= 70) {
            line <- cbind(as.vector(input[j,]), cosmic[k, ])
            line$PROP_GENE_COSMIC <- prop.gene.cosmic
            line$PROP_COSMIC_GENE <- prop.cosmic.gene
            int.matrix <- rbind(int.matrix, line)
          }
        }
      }
    }

    if(nrow(int.matrix) > 0) {
      return(int.matrix)
    }else{
      int.matrix <- .setRowMatrix(nrow = 0, c("Gene_Symbol", "log2FC.SCNAvsDip", "logCPM.SCNAvsDip", "p.val.SCNAvsDip", "FDR.SCNAvsDip", "TCGA_Tumor", "Condition", "Pat.percentage", "Pat.IDs", "Gene_Symbol_COSMIC", "log2FC.SCNAvsDip_COSMIC", "logCPM.SCNAvsDip_COSMIC", "p.val.SCNAvsDip_COSMIC", "FDR.SCNAvsDip_COSMIC", "TCGA_Tumor_COSMIC", "Condition_COSMIC", "Pat.percentage_COSMIC", "Pat.IDs_COSMIC", "PROP_GENE_COSMIC", "PROP_COSMIC_GENE"))
      return(int.matrix)
    }
  }
}

