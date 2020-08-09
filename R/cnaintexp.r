#' Perform joint expression and copy-number variation analyses on a set of genes using a set of tumor samples.
#'
#' @param genes List of official gene symbols for analysis. See function \code{all_genes} for a list of available genes.
#' @param cohorts List of tumors to study. See function \code{all_tumors} for a list of available tumors.
#' Defaults to all tumors with either an empty list or a list containing `ALL`.
#' @param writePath Path where results are stored. Defaults to the current folder.
#' @param pat.percentage Minimum percentage of patients per group.
#' @param pp.cor.cut Threshold to filter samples by AICC. Passed to `TCGAanalyze_Preprocessing`.
#' @param norm.method Method of normalization, such as `gcContent` or `geneLength` (default).
#' Passed to `TCGAanalyze_Normalization`.
#' @param filt.method Method of filtering, such as `quantile` (default), `varFilter`, `filter1`, `filter2`.
#' @param filt.qnt.cut Threshold selected as mean for filtering. Defaults to 0.25.
#' @param filt.var.func Filtering function. Defaults to `IQR`. See `genefilter` documentation for available methods.
#' @param filt.var.cutoff Threshold for `filt.var.funct`.
#' @param filt.eta Parameter for `filter1`. Defaults to 0.05.
#' @param filt.FDR.DEA Threshold to filter differentially expressed genes according their corrected p-value.
#' Passed to `TCGAanalyze_DEA`.
#' @param filt.FC Parameter for `filter2`. Defaults to 1.
#' @param p.val.thr Threshold for reported p values. Defaults to 1 (report all).
#' @param cna.thr Threshold level for copy-number variation analysis. Can be `Deep`, `Shallow` or `Both`.
#' `Deep`, consider likely homozygous deletions for loss and high-level amplifications (high copy-number, focal) for gain.
#' `Shallow`, consider likely heterozygous deletions for loss and low-level amplifications (lower copy-number, broad) for gain.
#' `Both`, consider all loss and gain events.
#' @param exp.mat Custom expression matrix, which will not be filtered. Defaults to `NULL`,
#' which means that the method considers TCGA expression matrices belonging to `cohorts`.
#' @param cna.mat Custom copy-number analysis matrix, which will not be filtered. Defaults to `NULL`,
#' which means that the method considers TCGA copy-number analysis matrices belonging to `cohorts`.
#'
#' @return List containing three data frames:
#' 1 - SCNA-mRNA correlations for queried genes
#' 2 - SCNA-mRNA correlations for Cancer Census genes (COSMIC) | Known oncodriver/TS genes.
#' 3 - SCNA overlapping between queried and CGC genes.
#' @export
#'
#' @examples
ciberAMP <- function(genes = c(),
                      cohorts = c(),
                      writePath = NULL,
                      pat.percentage = 10,
                      pp.cor.cut = 0.6,
                      norm.method = "geneLength",
                      filt.method = "quantile",
                      filt.qnt.cut = 0.25,
                      filt.var.func = "IQR",
                      filt.var.cutoff = 0.75,
                      filt.eta = 0.05,
                      filt.FDR.DEA = 0.01,
                      filt.FC = 1,
                      p.val.thr = 1,
                      cna.thr = "Deep",
                      exp.mat = NULL,
                      cna.mat = NULL) {

  list.of.packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "stringr", "RTCGAToolbox", "car", "EDASeq", "edgeR")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages) > 0) {BiocManager::install(new.packages)}

  cosmic.genes <- all_cosmic_genes()
  genes <- as.character(genes)
  sign <- c(genes, cosmic.genes)

  if (is.null(writePath)){
    writePath = getwd()
    print(paste("Defaulting write path to ", writePath, sep = ""))
  }

  setwd(writePath)

  tumors <- .tumors_all()
  tumors.with.normal <- .tumors_N()

  if(length(cohorts) == 0 || "ALL" %in% cohorts) {
    cohorts <- tumors
  }

  EXPintCNA.results <- NULL
  COSMIC.ov.result <- NULL

  for(i in 1:length(cohorts)) {
    tumor <- cohorts[i]
    if(tumor %in% tumors) {
      print(paste("Analyzing ", tumor, "...", sep=""))
    }else{
      print(paste("The introduced cohort: ", tumor, " was not found among the available TCGA cohorts. Bad-spelling is a common mistake, please check it or select any of the available ones.", sep=""))
      next
    }

    dataDEGs <- NULL
    if(is.null(exp.mat) && tumor %in% tumors.with.normal) {
    # If the user does not provide an expression matrix as indicated...
      tumor.exp <- .downloadExpression(tumor, tumors, tumors.with.normal)
      dataFilt <- .filterExpression(tumor, sign, tumor.exp, pp.cor.cut, norm.method,
                                    filt.method, filt.qnt.cut, filt.var.func,
                                    filt.var.cutoff, filt.eta, filt.FC)


      dataDEGs <- .getDataDEGs(tumor, dataFilt, filt.FDR.DEA, filt.FC)
      write.table(dataDEGs, file = paste("dataDEGs_", tumor, ".txt", sep=""), sep="\t", quote=FALSE)

    }else if(is.null(exp.mat) && tumor != tumors.with.normal){
      tumor.exp <- .downloadExpression(tumor, tumors, tumors.with.normal)
      dataFilt <- .filterExpression(tumor, sign, tumor.exp, pp.cor.cut, norm.method,
                                    filt.method, filt.qnt.cut, filt.var.func,
                                    filt.var.cutoff, filt.eta, filt.FC)

    }else{
      # If the user provides an expression matrix as indicated...
      dataFilt <- exp.mat
    }

    if(is.null(cna.mat)) {
      # If the user does not provide a SCNA matrix as indicated...
      gistic <- .getSCNAmatrix(tumor)

    }else if(!is.null(cna.mat)) {
      # If the user provides a SCNA matrix as indicated...
      gistic <- as.data.frame(cna.mat)

    }

    colnames(dataFilt) <- substr(colnames(dataFilt), 1, 12)

    dataFilt <- dataFilt[intersect(rownames(dataFilt), rownames(gistic)), intersect(colnames(dataFilt), colnames(gistic))]
    gistic <- gistic[intersect(rownames(dataFilt), rownames(gistic)), intersect(colnames(dataFilt), colnames(gistic))]

    exp <- as.data.frame(t(dataFilt[order(rownames(dataFilt)), ]))
    cna <- as.data.frame(t(gistic[order(rownames(gistic)), ]))

    exp <- exp[order(rownames(exp)), ]
    cna <- cna[order(rownames(cna)), ]

    save(exp, file = paste(tumor, "_exp_matrix.rda", sep=""))
    save(cna, file = paste(tumor, "_cna_matrix.rda", sep=""))

    SCNA.DEG.result <- .setRowMatrix(nrow = 0, c("Gene_Symbol", "log2FC.SCNAvsDip", "logCPM.SCNAvsDip", "p.val.SCNAvsDip", "FDR.SCNAvsDip", "TCGA_Tumor", "Condition", "Pat.percentage", "Pat.IDs"))

    for(j in 1:ncol(exp)) {
      gene <- colnames(exp)[j]
      new <- as.data.frame(.setRowMatrix(nrow(exp), c(paste(gene, "_exp", sep=""), paste(gene, "_cna", sep=""))))
      rownames(new) <- rownames(exp)
      new[,1] <- as.numeric(as.character(exp[,gene]))
      new[,2] <- as.numeric(as.character(cna[,gene]))

      group.del <- .selectDel(new, cna.thr)
      group.amp <- .selectAmp(new, cna.thr)
      group.neutro <- .selectDiploid(new, cna.thr)


      print(gene)
      print(paste("Deleted in ", nrow(group.del), " samples", sep=""))
      print(paste("Amplified in ", nrow(group.amp), " samples", sep=""))
      print(paste("Diploid in ", nrow(group.neutro), " samples", sep=""))
      print("------------------------")

      minimum.patients <- .setMinPat(new, pat.percentage)

      del.patients <- (nrow(group.del)/nrow(new)) * 100
      amp.patients <- (nrow(group.amp)/nrow(new)) * 100
      neutro.patients <- (nrow(group.neutro)/nrow(new)) * 100

      dataDEGs.SCNA <- NULL

      if(isTRUE(nrow(group.del) < minimum.patients & nrow(group.amp) < minimum.patients)) {

        next

      }else if(isTRUE(nrow(group.del) >= minimum.patients) & isTRUE(nrow(group.neutro) >= minimum.patients)) {

        group.x <- group.del
        group.y <- group.neutro

        cond <- c("Group DEL vs Group DIPLOID")

        SCNA.prop.pat <- del.patients

        pat.ids <- paste(rownames(group.del), collapse = ",")

        dataDEGs.SCNA <- .getDataDEGs_SCNA(tumor, dataFilt, group.x, group.y, filt.FDR.DEA, filt.FC, gene)

      }else if(isTRUE(nrow(group.amp) >= minimum.patients) & isTRUE(nrow(group.neutro) >= minimum.patients)) {

        group.x <- group.amp
        group.y <- group.neutro

        cond <- c("Group AMP vs Group DIPLOID")

        SCNA.prop.pat <- amp.patients

        pat.ids <- paste(rownames(group.amp), collapse = ",")

        dataDEGs.SCNA <- .getDataDEGs_SCNA(tumor, dataFilt, group.x, group.y, filt.FDR.DEA, filt.FC, gene)

      }

      if(!is.null(dataDEGs.SCNA)) {

        line <- .newSCNAline(dataDEGs.SCNA, cond, SCNA.prop.pat, pat.ids)
        SCNA.DEG.result <- rbind(SCNA.DEG.result, line)

      }
    }

    SCNA.DEG.result <- .convertToDF(SCNA.DEG.result)

    if(nrow(SCNA.DEG.result) > 0) {
      merge.dataDEGs <- .mergeDEGs(dataDEGs, SCNA.DEG.result, pat.percentage, genes, cosmic.genes)

      if(is.null(EXPintCNA.results)) {
        EXPintCNA.results <- merge.dataDEGs
      }else{
        EXPintCNA.results <- rbind(EXPintCNA.results, merge.dataDEGs)
      }
    }else{
      merge.dataDEGs <- .setRowMatrix(nrow = 0, c("Gene_Symbol", "logFC", "logCPM", "PValue", "FDR", "Tumor", "log2FC.SCNAvsDip", "logCPM.SCNAvsDip", "p.val.SCNAvsDip", "FDR.SCNAvsDip", "TCGA_Tumor", "Condition", "Pat.percentage", "Pat.IDs"))
      if(is.null(EXPintCNA.results)) {
        EXPintCNA.results <- merge.dataDEGs
      }else{
        EXPintCNA.results <- rbind(EXPintCNA.results, merge.dataDEGs)
      }
    }

    COSMIC.overlap <- .getOverlapCOSMIC(SCNA.DEG.result, genes, cosmic.genes)

    if(is.null(COSMIC.ov.result)) {
      COSMIC.ov.result <- COSMIC.overlap
    }else{
      COSMIC.ov.result <- rbind(COSMIC.ov.result, COSMIC.overlap)
    }

  }

  if(nrow(EXPintCNA.results) > 0 && genes %in% EXPintCNA.results$Gene_Symbol){
    write.table(EXPintCNA.results[EXPintCNA.results$Gene_Symbol %in% genes, ], "CiberAMP EXPintCNA results.txt", sep="\t", quote=FALSE)
    write.table(EXPintCNA.results[EXPintCNA.results$Gene_Symbol %in% cosmic.genes, ], "CiberAMP EXPintCNA COSMIC genes results.txt", sep="\t", quote=FALSE)
    write.table(COSMIC.ov.result, "CiberAMP COSMIC overlap results.txt", sep="\t", quote=FALSE)

    end <- list(EXPintCNA.results[EXPintCNA.results$Gene_Symbol %in% genes, ], EXPintCNA.results[EXPintCNA.results$Gene_Symbol %in% cosmic.genes, ], COSMIC.ov.result)
    return(end)
  }else{
    print("Any found")
  }
}

#' Get available genes
#'
#' @return List of available genes
#' @export
all_cosmic_genes <- function(){
  cosmic.genes <- c("A1CF","ACVR1","AKT1","AKT2","AKT3","AR","ARAF","ARHGAP5","BCL2L12","CACNA1D","CALR","CARD11","CCNE1","CCR4","CCR7","CD28","CD79A","CD79B","CDH17","CDK4","CHD4","CSF1R","CSF3R","CTNNA2","CTNND2","CXCR4","CYSLTR2","DDR2","DGCR8","EGFR","ERBB3","FGFR4","FLT3","FLT4","FOXA1","FUBP1","GATA2","GNA11","GNAQ","GNAS","GRM3","H3F3A","H3F3B","HIF1A","HIST1H3B","HRAS","IDH1","IDH2","IKBKB","IL6ST","IL7R","JAK3","JUN","KAT7","KCNJ5","KDR","KIT","KNSTRN","KRAS","MAP2K1","MAP2K2","MAPK1","MDM2","MDM4","MET","MITF","MPL","MTOR","MUC16","MUC4","MYCL","MYCN","MYD88","MYOD1","NRAS","NT5C2","PIK3CA","PIK3CB","PPM1D","PREX2","PRKACA","PTPN11","RAC1","REL","SALL4","SF3B1","SGK1","SKI","SMO","SOX2","SRC","SRSF2","STAT3","TNC","TRRAP","TSHR","U2AF1","USP8","WAS","XPO1","ZEB1","ABL1","ABL2","ACKR3","AFF3","AFF4","ALK","ATF1","BCL11A","BCL2","BCL3","BCL6","BCL9","BIRC6","BRAF","BRD3","BRD4","CCND1","CCND2","CCND3","CD74","CDK6","CHST11","CREB1","CREB3L2","CRLF2","CRTC1","CTNNB1","DDIT3","DDX5","DDX6","DEK","ELK4","ERBB2","ERG","ETV1","ETV4","ETV5","EWSR1","FCGR2B","FCRL4","FEV","FGFR1","FGFR2","FGFR3","FLI1","FOXP1","FOXR1","FSTL3","GLI1","HEY1","HIP1","HLF","HMGA1","HMGA2","HNRNPA2B1","HOXA13","HOXC11","HOXC13","HOXD11","HOXD13","JAK2","KAT6A","KDM5A","KMT2A","LCK","LMO1","LMO2","LPP","LYL1","MAF","MAFB","MALT1","MAML2","MECOM","MLLT10","MLLT4","MN1","MSI2","MTCP1","MYB","MYC","NCOA2","NFATC2","NPM1","NR4A3","NTRK3","NUP98","NUTM1","OLIG2","P2RY8","PAX3","PBX1","PDCD1LG2","PDGFB","PDGFRA","PDGFRB","PIM1","PLAG1","PLCG1","POU2AF1","POU5F1","PRDM16","PSIP1","RAF1","RAP1GDS1","RARA","RET","ROS1","RSPO3","SET","SETBP1","SH3GL1","SND1","SRSF3","SSX1","SSX2","SSX4","STAT6","STIL","SYK","TAF15","TAL1","TAL2","TCF7L2","TCL1A","TEC","TFE3","TFEB","TLX1","TLX3","TNFRSF17","TRIM27","USP6","WHSC1","WHSC1L1","WWTR1","ZNF521","APOBEC3B","ATP1A1","BCL9L","BCORL1","BMPR1A","BTK","CBLC","CDKN1A","CUX1","DAXX","DDB2","EPAS1","ERBB4","EZH2","FES","FOXL2","GATA1","GATA3","GPC3","IRS4","JAK1","KDM6A","KLF4","KMT2D","LEF1","MAP2K4","MAP3K1","MAP3K13","NFE2L2","NKX2-1","NOTCH2","PABPC1","POLQ","PTK6","QKI","RAD21","RECQL4","RHOA","TBX3","TERT","TP63","ARNT","BCL11B","BIRC3","CBL","CIC","CREBBP","ELF4","ESR1","FOXO1","FOXO3","FOXO4","HOXA11","HOXA9","IRF4","MALAT1","MKL1","NFKB2","NOTCH1","NTRK1","PAX5","PRKAR1A","RUNX1","RUNX1T1","STAT5B","SUZ12","TBL1XR1","TCF3","TET1","TP53","TRIM24","WT1","ACVR2A","AMER1","APC","ARID1B","ARID2","ASXL1","ASXL2","ATM","ATP2B3","ATR","ATRX","AXIN1","AXIN2","B2M","BAP1","BARD1","BLM","BRCA1","BRCA2","BRIP1","BUB1B","CASP8","CBLB","CCNC","CDC73","CDH1","CDH10","CDK12","CDKN1B","CDKN2A","CDKN2C","CEBPA","CHD2","CHEK2","CNOT3","CNTNAP2","CSMD3","CTCF","CUL3","CYLD","DDX3X","DICER1","DNM2","DNMT3A","DROSHA","EED","ELF3","ERCC2","ERCC3","ERCC4","ERCC5","ETNK1","EXT2","FAM46C","FANCA","FANCC","FANCD2","FANCE","FANCF","FANCG","FAS","FAT1","FAT4","FBLN2","FBXO11","FBXW7","FEN1","FH","FLCN","GRIN2A","HNF1A","ID3","KDM5C","KEAP1","KLF6","KMT2C","LARP4B","LRP1B","LZTR1","MAX","MED12","MEN1","MLH1","MSH2","MSH6","MUTYH","NBN","NCOR1","NCOR2","NF2","NFKBIE","NTHL1","PALB2","PBRM1","PHF6","PHOX2B","PIK3R1","PMS2","POLD1","POLE","POLG","POT1","PPP2R1A","PPP6C","PRDM1","PRDM2","PRF1","PTCH1","PTEN","PTPN13","PTPN6","PTPRB","PTPRC","PTPRT","RB1","RBM10","RNF43","ROBO2","RPL10","RPL5","SBDS","SDHA","SDHAF2","SDHB","SDHC","SDHD","SETD2","SFRP4","SH2B3","SIRPA","SMAD2","SMAD3","SMAD4","SMARCA4","SMARCB1","SMARCD1","SMARCE1","SMC1A","SOCS1","SPEN","SPOP","STAG1","STAG2","STK11","SUFU","TET2","TGFBR2","TMEM127","TNFAIP3","TNFRSF14","TRAF7","TSC1","TSC2","UBR5","VHL","WNK2","WRN","XPA","XPC","ZFHX3","ZMYM3","ZNRF3","ZRSR2","ABI1","ARHGAP26","ARHGEF12","ARID1A","BCL10","BCOR","BTG1","CAMTA1","CARS","CASC5","CBFA2T3","CBFB","CCDC6","CCNB1IP1","CD274","CDH11","CDX2","CIITA","CLTC","CLTCL1","CNBP","CREB3L1","DDX10","EBF1","EIF3E","ELL","EP300","EPS15","ETV6","EXT1","FHIT","FUS","IKZF1","KAT6B","LRIG3","MLF1","MYH9","NAB2","NCOA4","NDRG1","NF1","NRG1","PER1","PML","PPARG","PTPRK","RAD51B","RANBP2","RHOH","RMI2","RPL22","RSPO2","SFPQ")
  r <- c("SLC34A2","TPM3","TRIM33","WIF1","YWHAE","ZBTB16","ZNF278","ZNF331")
  cosmic.genes <- c(cosmic.genes, r)
  return(cosmic.genes)
}

#' Get available tumors
#'
#' @return List of available tumor samples
#' @export
all_tumors <- function(){
  non <- .tumors_noN()
  nor <- .tumors_N()
  result <- c(non, nor)
  return(result)
}

#' Plot CiberAMP results with ggplot
#' @param output List of results from CNAintEXP function.
#'
#' @return ggplot2 graph where Y axis = mRNA diff. expression between SCN-altered vs diploit tumors and X axis = mRNA diff. expression between Tumor and Normal tissue
#' @export
ggplot.CiberAMP <- function(df.exp){
  ggplot(df.exp, aes(x = logFC, y = log2FC.SCNAvsDip, color = TCGA_Tumor)) + geom_point(aes(size = Pat.percentage), shape = ifelse(df.exp$Condition %in% "Group AMP vs Group DIPLOID", 17, 16)) + scale_colour_manual(values = c("ACC" = "#ffcccc", "BLCA" = "#ffd9cc", "BRCA" = "#ffe6cc", "CHOL" = "#fff2cc", "COAD" = "#ffffcc", "CESC" = "#f2ffcc", "DLBC" = "#e6ffcc", "ESCA" = "#d9ffcc", "GBM" = "#ccffcc", "HNSC" = "#ccffd9", "KIRC" = "#ccffe6", "KIRP" = "#ccfff2", "KICH" = "#ccffff", "LAML" = "#ccf2ff", "LIHC" = "#cce6ff", "LGG" = "#ccd9ff", "LUAD" = "#ccccff", "LUSC" = "#d9ccff", "MESO" = "#e6ccff", "OV" = "#f2ccff", "PAAD" ="#ffccff", "PCPG" = "#ffccf2", "PRAD" = "#ffcce6", "READ" = "#ffccd9", "SARC" = "#ffcccc", "SKCM" = "#ff6666", "STAD" = "#b366ff",  "TGCT" = "#668cff", "THCA" = "#ff8c66", "THYM" = "#ff0000", "UCEC" = "#848785", "UCS" = "#767676", "UVM" = "#86C0C3")) + scale_size_continuous(range = c(1,10)) + theme_minimal() + xlab("mRNA diff. exp. tumor vs normal samples (log2(FC))") + ylab("mRNA diff. exp. SCN-altered vs diploid tumor samples (log2(FC))")
}



#' Interactive plot with ShinyR package
#' @param df First or second data frame from CNAintEXP list of results. Contains the correlation results for SCNA and mRNA differential expression for 1) queried or 2) COSMIC CGC genes.
#' @param int.df Third data frame from CNAintEXP list of results. Contains overlappings between SCN-altered queried and COSMIC CGC genes.
#'
#' @return It allows the user to directly interact with data using a shiny app
#' @export
int.plot.CiberAMP <- function(df, int.df){

  list.of.packages <- c("shiny", "plotly", "DT", "dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages) > 0) {install.packages(new.packages)}

  require(shiny)
  require(plotly)
  require(DT)
  require(dplyr)

  cohorts.list = list("BLCA" = "BLCA", "ACC" = "ACC", "BRCA" = "BRCA", "CESC" = "CESC", "CHOL" = "CHOL", "COAD" = "COAD", "DLBC" = "DLBC", "ESCA" = "ESCA", "GBM" = "GBM", "HNSC" = "HNSC", "KICH" = "KICH", "KIRC" = "KIRC", "KIRP" = "KIRP", "LAML" = "LAML", "LGG" = "LGG", "LIHC" = "LIHC", "LUAD" = "LUAD", "LUSC" = "LUSC", "MESO" = "MESO", "OV" = "OV", "PAAD" = "PAAD", "PRAD" = "PRAD", "READ" = "READ", "SARC" = "SARC", "SKCM" = "SKCM", "STAD" = "STAD", "TGCT" = "TGCT", "THCA" = "THCA", "THYM" = "THYM", "UCEC" = "UCEC", "UCS" = "UCS", "UVM" = "UVM" )
  tumors <- .tumors_all()
  tumors.with.normal <- .tumors_N()
  df$KEY <- paste(df$Gene_Symbol, df$TCGA_Tumor, sep="_")
  # Define UI for application that draws a histogram
  ui <- fluidPage(

    # Application title
    titlePanel("CiberAMP"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(
        textAreaInput("genes", "Enter your genes", width = "100%", height = "250px", value = "ALL"),
        actionButton("submit", "Submit"),
        sliderInput("pat.perc", "Select the minimum copy number altered samples", value = 10, min = 0, max = 100),
        numericInput("p.val.thr", "Select the adjusted p-value threshold.", value = 0.05, min = 0, max = 1),
        numericInput("tvsn.fc.thr", "Select the minimum fold-change threshold when compare tumor vs normal samples.", value = 0, min = 0, max = Inf),
        numericInput("cna.fc.thr", "Select the minimum fold-change threshold when compare CNA vs diploid samples.", value = 0, min = 0, max = Inf),
        width = 3
      ),

      # Show a plot of the generated distribution
      mainPanel(
        br(),
        plotlyOutput("plot", width = "100%", height = "100%"),
        br(),
        br(),
        conditionalPanel(condition = "!is.null(cosmic)", DT::dataTableOutput("cosmic")),
        br()
      )
    )
  )

  # Define server logic required to draw a histogram
  server <- function(input, output, session) {

    g <- eventReactive(input$submit, {
      unlist(strsplit(gsub("\n", ",", req(input$genes)), split = ","))
    })

    d <- reactive({
      if("ALL" %in% g()) {
        df[abs(df$logFC) >= input$tvsn.fc.thr & abs(df$log2FC.SCNAvsDip) >= input$cna.fc.thr & df$FDR.SCNAvsDip <= input$p.val.thr & df$Pat.percentage >= input$pat.perc, ]
      }else{
        df[df$Gene_Symbol %in% g() & abs(df$logFC) >= input$tvsn.fc.thr & abs(df$log2FC.SCNAvsDip) >= input$cna.fc.thr & df$FDR.SCNAvsDip <= input$p.val.thr & df$Pat.percentage >= input$pat.perc, ]
      }
    })

    output$plot <- plotly::renderPlotly({ plotly::plot_ly(d(), x = ~logFC, y = ~log2FC.SCNAvsDip, source = "plot", key = ~KEY, text = ~paste("Symbol:", Gene_Symbol, "<br>Condition:", Condition, "<br>Tumor:", TCGA_Tumor, "<br>CNA samples (%):", Pat.percentage), type = "scatter", mode = "markers", marker = list(size = ~Pat.percentage, sizes = c(1, 100), opacity = 0.8, line = list(color = "white", width = 2)), color = ~TCGA_Tumor) %>%
        layout(title = "CiberAMP plot for TCGA cohorts", xaxis = list(title = "mRNA DE tumor vs normal samples (log2(FC))", showgrid = TRUE), yaxis = list(title = "mRNA DE SCN-altered vs diploid tumor samples (log2(FC))", showgrid = TRUE))
    })

    output$cosmic <- DT::renderDataTable({
      a <- event_data("plotly_click", source = "plot")
      if(length(a)) {
        a <- a[["key"]]
        a <- strsplit(as.character(a), split = "_")
        a <- unlist(a)
        cosmic <- int.df[int.df$Gene_Symbol %in% a[1] & int.df$TCGA_Tumor %in% a[2], c("Gene_Symbol", "TCGA_Tumor", "Gene_Symbol_COSMIC", "PROP_GENE_COSMIC", "PROP_COSMIC_GENE")]
        if(nrow(cosmic) > 0) {
          DT::datatable(cosmic, options = list(scrollX = TRUE), rownames = FALSE) %>% formatRound(c("PROP_GENE_COSMIC", "PROP_COSMIC_GENE"), 2)
        }else{
          NULL
        }
      }else if(!length(a)) {
        NULL
      }
    })
  }

  shinyApp(ui = ui, server = server)

}

