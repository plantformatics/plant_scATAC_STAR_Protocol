###################################################################################################
#                             Socrates QC and clustering example                                  #
###################################################################################################


# System requirements ----------------------------------------------------
# requires R v4.0.0 or greater
# requires MACS2 in path


# load libraries ---------------------------------------------------------
library(Socrates)


# load data --------------------------------------------------------------
bed <- "Zm_B73_seedling.rep1.unique_Tn5.bed.gz" # 5 columns: chr start end barcode strand
ann <- "Zea_mays.AGPv4.36.gtf" # must end in either gtf or gff/gff3
ref <- "Zea_mays.AGPv4.36.fa.fai" # reference genome index (output from samtools faidx)



###################################################################################################
#                                    Data quality analysis                                        #
###################################################################################################


# create Socrates QC object ----------------------------------------------
obj <- loadBEDandGenomeData(bed, 
                            ann, 
                            ref, 
                            attribute="gene_id") # specify attribute to use for geneIDs in gff/gtf


# Identify accessible chromatin regions in bulk --------------------------
obj <- callACRs(obj, genomesize=1.6e9, 
                shift= -50, 
                extsize=100,
                fdr=0.05,
                output="Zm_B73_seedling.rep1.bulk_ACRs", 
                tempdir="./Zm_B73_seedling.rep1.bulk_ACRs_MACS2", 
                verbose=T)


# build metadata ---------------------------------------------------------
obj <- buildMetaData(obj, 
                     tss.window=2000, # users may wish to adjust window size based on the amount of intergenic space
                     verbose=TRUE)


# filter cells -----------------------------------------------------------
obj <- findCells(obj, 
                 doplot=T) # plot the results


# generate sparse matrix -------------------------------------------------
obj <- generateMatrix(obj, 
                      filtered=T, # use the filtered cells from `findCells`
                      windows=1000, # bin size (1-kb)
                      peaks=F, # set to TRUE to use the bulk-level ACRs to build the sparse matrix
                      verbose=T)


# convert to Socrates format for downstream analysis ---------------------
soc.obj <- convertSparseData(obj, 
                             verbose=T)


# save QC object ---------------------------------------------------------
saveRDS(obj, file="Zm_B73_seedling.rep1.QC_results.rds")



###################################################################################################
#                                    Clustering analysis                                          #
###################################################################################################


# get per cell feature counts --------------------------------------------
cell.counts <- Matrix::colSums(soc.obj$counts) # count number of features with Tn5 insertions per cell
cell.counts.z <- as.numeric(scale(cell.counts)) # convert features counts into Z-scores
cell.counts.threshold <- max(c(cell.counts[cell.counts.z < -1], 1000)) # minimum feature counts (greater of 1 std or 1000)


# clean sparse counts matrix ---------------------------------------------
soc.obj <- cleanData(soc.obj, 
                     min.c=cell.counts.threshold, # minimum number of accessible features per cell
                     min.t=0.01, # minimum feature frequency across cells
                     max.t=0.001, # maximum feature frequency across cells
                     verbose=T)


# normalize with regularized quasibinomial logistic regression -----------
soc.obj <- regModel(soc.obj, 
                    verbose=T)


# denoise with SVD -------------------------------------------------------
soc.obj <- reduceDims(soc.obj, 
                      n.pcs=20, # number of components to retain
                      cor.max=0.7, # filter components correlated above this value with per cell feature counts
                      verbose=T)


# reduce to 2-dimensions with UMAP ---------------------------------------
soc.obj <- projectUMAP(soc.obj, 
                       verbose=T)


# identify clusters using neighborhood graph -----------------------------
soc.obj <- callClusters(soc.obj, 
                        res=0.8, # resolution of clusters (lower = less clusters)
                        verbose=T)


# plot cluster membership on UMAP embedding ------------------------------
plotUMAP(soc.obj)


# save results -----------------------------------------------------------
saveRDS(soc.obj, file="Zm_B73_seedling.rep1.clustering_results.rds")



##################################
## explanation of soc.obj slots ##
##################################

## soc.obj$clusters - contains metadata for each barcode and the cluster membership
## soc.obj