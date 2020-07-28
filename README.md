
------------------------------------------------------------------------

# CiberAMP | To discover copy number alteration-driven mRNA differential expression in TCGA tumors.

CiberAMP integrates somatic copy number alterations (SCNAs) and mRNA expression data from TCGA cohorts. It consists in two differential expression analysis (DEA): 1) mRNA differetial abundance between tumor and normal samples (if available) and 2) mRNA differential expression between copy number-altered and diploid tumor samples. First, samples are classified as tumors or normal and, secondly, tumor samples are classified as deeply/shallowly amplified/deleted based on GISTIC2.0 run on TCGA cohorts.

Finally, since highly copy number-altered oncodrivers can originate others genes alterations, CiberAMP calculated the % of SCN-altered samples between your genes and the COSMIC Cancer Gene Census (CGC) list. This returns two values: 1) your gene over COSMIC gene and 2) COSMIC gene over your gene.

### Installation from GitHub ###

```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("r-lib/devtools")

# Install ciberAMP by devtools
devtools::install_github("ciberAMP", dependencies = TRUE)
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

* *genes* Your genes of interest official symbols. This should be a character vector
* *cohorts* By default, CiberAMP will compute every TCGA cohort (the 32 available), but you can use this argument to only analyze a few. You can consult the offical codes at TCGA's web -> https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
* *writePath* The folder path to save the results. TIP: if you want to re-run CiberAMP, if you use the same folder where all data was stored, you will never need to re-download it again what saves a lot of time and disk space! But, be careful, every run your results will be overwritten!
* *pat.percentage* The minimum % of altered patients that should be considered to look for SCNA-driven differential expression. When in a cohort an event occurs recurrently it is more likely that such event may be under positive selection so being relevant for tumor fitness.
* *pp.cor.cut* Threshold to filter samples by AICC. Passed to `TCGAanalyze_Preprocessing`.
* *norm.method* Method of normalization, such as `gcContent` or `geneLength` (default).
* *filt.method* Method of filtering, such as `quantile` (default), `varFilter`, `filter1`, `filter2`.
* *filt.qnt.cut* Threshold selected as mean for filtering. Defaults to 0.25.
* *filt.var.func* Filtering function. Defaults to `IQR`. See `genefilter` documentation for available methods.
* *filt.var.cutoff* Threshold for `filt.var.funct`.
* *filt.eta* Parameter for `filter1`. Defaults to 0.05.
* *filt.FDR.DEA* Threshold to filter differentially expressed genes according their corrected p-value.
* *filt.FC* Parameter for `filter2`. Defaults to 1.
* *p.val.thr* Threshold for reported p values. Defaults to 1 (report all).
* *cna.thr* Threshold level for copy-number variation analysis. Can be `Deep`, `Shallow` or `Both`
* *exp.mat* Custom expression matrix, which will not be filtered. Defaults to `NULL`,
* *cna.mat* Custom copy-number analysis matrix, which will not be filtered. Defaults to `NULL`,

To better understand many of these values you can visit TCGAbilolinks web tutorial: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html

------------------------------------------------------------------------

# Looking into CiberAMP results

CiberAMP results into a list of 3 data frames that can be accessed by:

The x[[1]] data frame contains differential expression results for queried genes:

* Column 1 -> queried gene approved symbol.
* Columns 2:4 -> queried genes tumor vs normal differential expression results.
* Column 5 -> TCGA cohort.
* Column 6:9 -> queried gene copy number altered vs diploid tumor samples differential expression results.
* Column 10  -> TCGA cohort.
* Column 11  -> CN-altered vs diploid tumors comparison type: amplified vs diploid or deleted vs diploid.
* Column 12  -> % of samples CN-altered.
* Column 14 -> CN-altered TCGA sample barcodes.

The x[[2]] data frame contains differential expression results for COSMIC cancer census genes list. Column information is exactly the same as x[[1]].

Finally, the x[[3]] contains overlapping information between queried and CGC genes copy number overlap between samples. Genes may be CN-altered due to being lodged next to well-known and transcriptionally active in cancer cancer-related genes. In this data frame, you can consult if any of your queried CN-alterations would be highly correlated to any of known oncodrivers or tumor suppressors:

* Column 1 -> queried gene approved symbol.
* Column 2:5 -> queried gene copy number altered vs diploid tumor samples differential expression results.
* Column 6 -> TCGA cohort.
* Column 7 -> CN-altered vs diploid tumors comparison type: amplified vs diploid or deleted vs diploid.
* Column 8 -> % of samples CN-altered in queried gene.
* Column 9 -> CN-altered TCGA sample barcodes in queried gene.
* Column 10 -> COSMIC CGC gene approved symbol.
* Column 11:14 -> COSMIC CGC gene copy number altered vs diploid tumor samples differential expression results.
* Column 15 -> TCGA cohort.
* Column 16 -> CN-altered vs diploid tumors comparison type: amplified vs diploid or deleted vs diploid.
* Column 17 -> % of samples CN-altered in COSMIC CGC gene.
* Column 18 -> CN-altered TCGA sample barcodes in COSMIC CGC gene.
* Column 19 -> % of queried gene CN-altered samples that overlaps with COSMIC CGC gene CN-altered samples.
* Column 20 -> % of COSMIC CGC gene CN-altered samples that overlaps with queried gene CN-altered samples.

where x is the variable designated while calling ciberAMP function.
