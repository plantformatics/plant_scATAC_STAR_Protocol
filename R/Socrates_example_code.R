######################################
# Socrates QC and clustering example #
######################################

# requires R v4.0.0 or greater
# requires MACS2 in path

# load libraries
library(Socrates)

# load data 
bed <- "Zm_B73_seedling.rep1.unique_Tn5.bed.gz" # 5 columns: chr start end barcode strand
ann <- "Zea_mays.AGPv4.36.gtf" # must end in either gtf or gff/gff3
ref <- "Zea_mays.AGPv4.36.fa.fai" # reference genome index (output from samtools faidx)


###############
# QC analysis #
###############

# create Socrates QC object
obj <- loadBEDandGenomeData(bed, 
                            ann, 
                            ref, 
                            attribute="gene_id") # specify attribute to use for geneIDs in gff/gtf

# Identify accessible chromatin regions in bulk (this step may take a way for large BED files)
obj <- callACRs(obj, genomesize=1.6e9, 
                shift= -50, 
                extsize=100,
                fdr=0.05,
                output="Zm_B73_seedling.rep1.bulk_ACRs", 
                tempdir="./Zm_B73_seedling.rep1.bulk_ACRs_MACS2", 
                verbose=T)

# build metadata
obj <- buildMetaData(obj, 
                     tss.window=2000, # users may wish to adjust window size based on the amount of intergenic space
                     verbose=TRUE)

# filter cells
obj <- findCells(obj, 
                 doplot=T) # plot the results

# generate sparse matrix
obj <- generateMatrix(obj, 
                      filtered=T, # use the filtered cells from `findCells`
                      windows=1000, # bin size (1-kb)
                      peaks=F, # set to TRUE to use the bulk-level ACRs to build the sparse matrix
                      verbose=T)

# convert to Socrates format for downstream analysis. 
soc.obj <- convertSparseData(obj, 
                             verbose=T)

# save QC object
saveRDS(obj, file="Zm_B73_seedling.rep1.QC_results.rds")


#######################
# Clustering analysis #
#######################

