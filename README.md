
------------------------------------------------------------------------

# CiberAMP | An R package to identify differential mRNA expression linked to somatic copy number variations in cancer datasets

CiberAMP is an R package that uses differential expression analyses to stablish accurate correlations between specific SCNVs and changes in expression in the genes affected by them. The algorithm has been designed to be an easy-to-access tool for the TCGA, the largest database in the world with genomic and transcriptomic data ofr more than 10,000 samples of 33 different human cancers.

Unlike other methods, CiberAMP can yield information on:
  (i) SCNV-DEGs (somatic copy number variations associated differentially expressed genes) in a cohort of TCGA tumor samples
  (ii) The type of copy number variation associated with each SCNV-DEG in terms of expression pattern and genomic context
  (iii) Insights on the potential functional relevance of each identified SCNV-DEG

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

* *genes* The list of genes of interest. It is a vector of gene official symbols according to the HGNC.
* *cohorts* The list of TCGA cohorts to be analyzed. By default, CiberAMP will be run on all TCGA cohorts. You can consult the official TCGA cohort IDs here: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations or in CiberAMP's manuscript (Table S1)
* *writePath* The path to the folder to save results TIP: if you want to re-run CiberAMP, if you use the same folder where all data was stored, you will not need re-download all data from the TCGA again. This can save you a lot of space in your disk, but be careful, results will be overwritten as well.
* *pat.percentage* The minimum % of copy number altered samples per gene that will be analyzed.
* *pp.cor.cut* Threshold to filter samples by array-array intensity correlation (AICC) analysis (0.6 by default). Passed to `TCGAanalyze_Preprocessing`.
* *norm.method* Method of normalization, such as `gcContent` or `geneLength` (default). See TCGAbiolinks R package for help.
* *filt.method* Method of filtering, such as `quantile` (default), `varFilter`, `filter1`, `filter2`. See TCGAbiolinks R package for help.
* *filt.qnt.cut* Threshold selected as quantile for filtering. Defaults to 0.25 (first quantile).
* *filt.var.func* Filtering function. Defaults to `IQR`. See `genefilter` documentation for available methods.
* *filt.var.cutoff* Threshold for `filt.var.funct`. See TCGAbiolinks R package for help.
* *filt.eta* Parameter for `filter1`. Defaults to 0.05. See TCGAbiolinks R package for help.
* *filt.FDR.DEA* Threshold to filter differentially expressed genes according their corrected p-value.
* *filt.FC* Minimum log2(FC) value to considered a gene as differentially expressed. Defaults to 0.58 (that corresponds to a differential expression of at least 50%).
* *cna.thr* Threshold level for copy-number variation analysis. Can be `Deep`, `Shallow` or `Both`
* *exp.mat* Custom normalized RNAseq counts expression matrix of only tumors. Defaults to `NULL`.
* *cna.mat* Custom copy-number analysis matrix of only tumors. Defaults to `NULL`.

------------------------------------------------------------------------

# Looking into CiberAMP results

CiberAMP returns a list of 3 data frames:

The first data frame contains all SCNV-DEGs and genes differentially expressed between tumor and normal samples exclusively. The secon data frame contains all the SCNV-DE known cancer drivers. These two data frames have the same format and in each column we can find:

* Column 1 -> Gene approved symbols
* Columns 2:4 -> Results from the differential expression analysis between tumor and healthy samples.
* Column 5 -> ID of the queried TCGA cohort.
* Column 6:9 -> Results from the differential expression analysis between copy number altered and diploid tumor samples.
* Column 10  -> ID of the queried TCGA cohort.
* Column 11  -> The type of comparison made: amplified vs. diploid or deleted vs. diploid.
* Column 12  -> Recurrence of gene amplifications or deletions in the cohort.
* Column 14 -> Barcodes of the samples harboring such SCNVs.

The third data frame contains the information about the significant co-occurring amplification or deletions between the SCNV-DEGs and known cancer drivers:

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

------------------------------------------------------------------------

# Looking into CiberAMP's logic classifier results

The logic classification algorithm integrated in CiberAMP's package allows the user to rate the potential candidates subdividing them into four subgroups.

First, the SCN-associated DEGs reported from the previous step are divided based on their significant genomic interactions with any COSMIC CGC oncogene in each cohort.
Secondly, these genes are further subdivided regarding their genomic location inside or outside enriched genomic regions. 
Finally, within each of the four resulting subgroups, genes are rated based on, first, their recurrency and, secondly, their SCN-associated FDR adjusted p-value.

```r
# Load the library
library(ciberAMP)

# Write your function
x <- CiberAMP.classifier(res1 = NULL, res3 = NULL, width.window = 6000000)
```
Where:
* *res1* The first data frame reported from the previous function
* *res3* The third data frame reported from the previous function
* *width.window* The window length in base pairs used for genomic enriched clusters calculation.

The outcomes of this function is a list of 4 data frames. The first data frame contains all the SCNV-DEGs that are not co-amplified or co-deleted with any known cancer driver gene and outside any enriched cluster. The second data frame conatins all SCNV-DEGs that are not co-amplified or co-deleted with any known cancer driver gene and located within an enriched cluster. The third data frame containes all SCNV-DEGs that are co-amplified or co-deleted with a known cancer driver gene and outside any enriched cluster. The fourth data frame contains all SCNV-DEGs that are co-amplified or co-deleted with a known cancer driver gene and within an enriched gene cluster.
