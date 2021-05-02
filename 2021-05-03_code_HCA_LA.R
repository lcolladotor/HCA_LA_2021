## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----"biocmanager", eval = FALSE---------------------------------------------------------------------------------------------------------
## ## Install BiocManager which is the utility for installing
## ## Bioconductor packages
## 
## ## For installing Bioconductor packages
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")


## ----"biocversion"-----------------------------------------------------------------------------------------------------------------------
BiocManager::version()


## ----"install_deps", eval = FALSE--------------------------------------------------------------------------------------------------------
## ## Install required packages
## BiocManager::install(
##     c(
##         'SingleCellExperiment',
##         # 'usethis',
##         # 'here',
##         'scran',
##         'scater',
##         'scRNAseq',
##         # 'org.Mm.eg.db',
##         'AnnotationHub',
##         'ExperimentHub',
##         # 'BiocFileCache',
##         # 'DropletUtils',
##         # 'EnsDb.Hsapiens.v86',
##         # 'TENxPBMCData',
##         # 'BiocSingular',
##         # 'batchelor',
##         'uwot',
##         # 'Rtsne',
##         # 'pheatmap',
##         # 'fossil',
##         # 'ggplot2',
##         # 'cowplot',
##         # 'RColorBrewer',
##         # 'plotly',
##         # 'iSEE',
##         'pryr',
##         'sessioninfo'
##     )
## )


## ----"data infraestructure code"---------------------------------------------------------------------------------------------------------
## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")

# Load the SingleCellExperiment package
library('SingleCellExperiment')
# Extract the count matrix from the 416b dataset
counts.416b <- counts(sce.416b)
# Construct a new SCE from the counts matrix
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Inspect the object we just created
sce

## How big is it?
pryr::object_size(sce)

# Access the counts matrix from the assays slot
# WARNING: This will flood RStudio with output!

# 1. The general method
assay(sce, "counts")[1:6, 1:3]
# 2. The special method for the assay named "counts"
counts(sce)[1:6, 1:3]

sce <- scater::logNormCounts(sce)
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# 1. The general method
assay(sce, "logcounts")[1:6, 1:3]
# 2. The special method for the assay named "logcounts"
logcounts(sce)[1:6, 1:3]

# assign a new entry to assays slot
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# List the assays in the object
assays(sce)
assayNames(sce)

## How big is it?
pryr::object_size(sce)

# Extract the sample metadata from the 416b dataset
colData.416b <- colData(sce.416b)
# Add some of the sample metadata to our SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Inspect the object we just updated
sce
# Access the sample metadata from our SCE
colData(sce)
# Access a specific column of sample metadata from our SCE
table(sce$block)

# Example of function that adds extra fields to colData
sce <- scater::addPerCellQC(sce.416b)
# Access the sample metadata from our updated SCE
colData(sce)

# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

## Add the lognorm counts again
sce <- scater::logNormCounts(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just wild type cells
# Remember, cells are columns of the SCE
sce[, sce$phenotype == "wild type phenotype"]

# Access the feature metadata from our SCE
# It's currently empty!
rowData(sce)

# Example of function that adds extra fields to rowData
sce <- scater::addPerFeatureQC(sce)
# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)


# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Annotate each gene with its chromosome location
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
rowData(sce)$chromosome <- chromosome

# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just genes on chromosome 3
# NOTE: which() needed to cope with NA chromosome names
sce[which(rowData(sce)$chromosome == "3"), ]

# Access the metadata from our SCE
# It's currently empty!
metadata(sce)

# The metadata slot is Vegas - anything goes
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
    analyst = c("Pete"))

# Access the metadata from our updated SCE
metadata(sce)

# E.g., add the PCA of logcounts
# NOTE: We'll learn more about PCA later
sce <- scater::runPCA(sce)
# Inspect the object we just updated
sce
# Access the PCA matrix from the reducedDims slot
reducedDim(sce, "PCA")[1:6, 1:3]

# E.g., add a t-SNE representation of logcounts
# NOTE: We'll learn more about t-SNE later
sce <- scater::runTSNE(sce)
# Inspect the object we just updated
sce
# Access the t-SNE matrix from the reducedDims slot
head(reducedDim(sce, "TSNE"))

# E.g., add a 'manual' UMAP representation of logcounts
# NOTE: We'll learn more about UMAP later and a
# 		  simpler way to compute it.
u <- uwot::umap(t(logcounts(sce)), n_components = 2)
# Add the UMAP matrix to the reducedDims slot
# Access the UMAP matrix from the reducedDims slot
reducedDim(sce, "UMAP") <- u

# List the dimensionality reduction results stored in # the object
reducedDims(sce)


## ----"quality_control"-------------------------------------------------------------------------------------------------------------------
## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
## Data
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
# Annotate each gene with its chromosome location
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SEQNAME"
)
# Identify the mitochondrial genes
is.mito <- which(location == "MT")

library('scater')
sce.416b <- addPerCellQC(sce.416b,
    subsets = list(Mito = is.mito))


## ----qc_metrics, cache=TRUE, dependson='all_code'--------------------------------------------------------------------
plotColData(sce.416b, x = "block", y = "detected")

plotColData(sce.416b, x = "block", y = "detected") +
    scale_y_log10()

plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)


## ----all_code_part2, cache = TRUE, dependson='all_code'--------------------------------------------------------------
# Example thresholds
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib),
    NExprs = sum(qc.nexprs),
    SpikeProp = sum(qc.spike),
    MitoProp = sum(qc.mito),
    Total = sum(discard)
)

qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(sce.416b$detected, log = TRUE,
    type = "lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extract the thresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
# Summarize the number of cells removed for each reason.
DataFrame(
    LibSize = sum(qc.lib2),
    NExprs = sum(qc.nexprs2),
    SpikeProp = sum(qc.spike2),
    MitoProp = sum(qc.mito2),
    Total = sum(discard2)
)

## More checks
plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)

batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
qc.lib3 <- isOutlier(sce.416b$sum,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.nexprs3 <- isOutlier(sce.416b$detected,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher",
    batch = batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher",
    batch = batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extract the thresholds
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib3),
    NExprs = sum(qc.nexprs3),
    SpikeProp = sum(qc.spike3),
    MitoProp = sum(qc.mito3),
    Total = sum(discard3)
)


## ----use_case, cache=TRUE, dependson= c('all_code', 'all_code_part2')------------------------------------------------
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")

discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor)
discard.ercc2 <- isOutlier(
    sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor,
    subset = sce.grun$donor %in% c("D17", "D2", "D7")
)

plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc)
)
plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc2)
)

# Add info about which cells are outliers
sce.416b$discard <- discard2

# Look at this plot for each QC metric
plotColData(
    sce.416b,
    x = "block",
    y = "sum",
    colour_by = "discard",
    other_fields = "phenotype"
) +
    facet_wrap( ~ phenotype) +
    scale_y_log10()

# Another useful diagnostic plot
plotColData(
    sce.416b,
    x = "sum",
    y = "subsets_Mito_percent",
    colour_by = "discard",
    other_fields = c("block", "phenotype")
) +
    facet_grid(block ~ phenotype)


## ----"session_info"----------------------------------------------------------------------------------------------------------------------
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


## ----"tangle", eval = FALSE--------------------------------------------------------------------------------------------------------------
## knitr::knit("2021-05-03_code_HCA_LA.Rmd", tangle = TRUE)

