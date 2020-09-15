
------------------------------------------------------------------------

# CiberAMP | Empowering somatic copy number-driven mRNA differential expression in TCGA tumors.

CiberAMP integrates somatic copy number alterations (SCNAs) and mRNA expression data from the TCGA cohorts. It consists in two differential expression analysis (DEA): 1) between tumor vs. normal samples (if available) and 2) SCN-altered and diploid tumor samples. Tumor samples are classified as deeply/shallowly amplified/deleted or diploid based on GISTIC2.0 thresholded-by-gene file.

Finally, CiberAMP calculates the % of SCN-altered samples between your queried genes and the COSMIC Cancer Gene Census (CGC) list of known oncodrivers.

### Installation from GitHub ###

```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("r-lib/devtools")

# Install ciberAMP by devtools
devtools::install_github("vqf/ciberAMP", dependencies = TRUE)
```

------------------------------------------------------------------------

### Usage ###

```r
# Load the library
library(ciberAMP)

# Write your function
x <- ciberAMP(genes = c(), cohorts = c(), pat.percentage = 0, writePath = "PATH_TO_FOLDER")
```

Where:

* *genes* Your genes of interest approved symbols. This should be a character vector.
* *cohorts* By default, CiberAMP will compute every TCGA cohort (the 32 available), but you can use this argument to only analyze a few. You can consult the offical codes at TCGA's web: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
* *writePath* The folder path to save the results. TIP: if you want to re-run CiberAMP, if you use the same folder where all data was stored, you will never need to re-download it again what saves a lot of time and disk space! But, be careful, your results will be overwritten in every run!
* *pat.percentage* The minimum % of SCN-altered patients that should be considered to look for SCNA-driven differential expression.
* *pp.cor.cut* Threshold to filter samples by array-array intensity correlation (AICC) analysis. Passed to `TCGAanalyze_Preprocessing`.
* *norm.method* Method of normalization, such as `gcContent` or `geneLength` (default). See TCGAbiolinks R package for help.
* *filt.method* Method of filtering, such as `quantile` (default), `varFilter`, `filter1`, `filter2`. See TCGAbiolinks R package for help.
* *filt.qnt.cut* Threshold selected as quantile for filtering. Defaults to 0.25 (first quantile).
* *filt.var.func* Filtering function. Defaults to `IQR`. See `genefilter` documentation for available methods. See TCGAbiolinks R package for help.
* *filt.var.cutoff* Threshold for `filt.var.funct`. See TCGAbiolinks R package for help.
* *filt.eta* Parameter for `filter1`. Defaults to 0.05. See TCGAbiolinks R package for help.
* *filt.FDR.DEA* Threshold to filter differentially expressed genes according their corrected p-value.
* *filt.FC* Minimum log2(FC) value to considered a gene as differentially expressed. Defaults to 1.
* *cna.thr* Threshold level for copy-number variation analysis. Can be `Deep`, `Shallow` or `Both`
* *exp.mat* Custom normalized RNAseq counts expression matrix of only tumors. Defaults to `NULL`.
* *cna.mat* Custom copy-number analysis matrix of only tumors. Defaults to `NULL`.

------------------------------------------------------------------------

# Looking into CiberAMP results

CiberAMP results into a list of 3 data frames that can be accessed by:

The x[[1]] data frame containing queried genes differential expression results:

* Column 1 -> queried gene approved symbol.
* Columns 2:4 -> queried genes tumor vs normal differential expression results.
* Column 5 -> TCGA cohort ID.
* Column 6:9 -> queried gene SCN-altered vs. diploid tumor samples differential expression results.
* Column 10  -> TCGA cohort.
* Column 11  -> SCN-altered vs. diploid tumor samples comparison: amplified vs. diploid or deleted vs. diploid.
* Column 12  -> % of samples SCN-altered.
* Column 14 -> SCN-altered TCGA sample barcodes.

The x[[2]] data frame containing CGC list differential expression results. Columns are the same as x[[1]].

Finally, the x[[3]] data frame containing the % of tumor samples in which queried gene's SCN-driven differential expression co-occurs with any CGC oncodrivers' SCNAs:

* Column 1 -> queried gene approved symbol.
* Column 2:5 -> queried gene copy number altered vs. diploid tumor samples DE results.
* Column 6 -> TCGA cohort ID.
* Column 7 -> SCN-altered vs. diploid tumor samples comparison type: amplified vs. diploid or deleted vs. diploid.
* Column 8 -> % of tumor SCN-altered samples for a queried gene.
* Column 9 -> queried gene SCN-altered TCGA tumor samples' barcodes.
* Column 10 -> COSMIC CGC gene approved symbol.
* Column 11:14 -> COSMIC CGC gene copy number altered vs. diploid tumor samples DE results.
* Column 15 -> TCGA cohort ID.
* Column 16 -> SCN-altered vs. diploid tumor samples comparison type: amplified vs diploid or deleted vs diploid.
* Column 17 -> % of tumor SCN-altered samples for a COSMIC CGC gene.
* Column 18 -> COSMIC CGC gene SCN-altered TCGA tumor samples' barcodes.
* Column 19 -> % of overlapping queried and COSMIC CGC genes SCN-altered tumor samples.
* Column 20 -> % of overlapping COSMIC CGC and queried genes SCN-altered tumor samples.

where x is the variable designated to recall CiberAMP outcomes.
