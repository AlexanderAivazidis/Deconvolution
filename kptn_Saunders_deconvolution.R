### Kptn bulk data deconvolution using Saunders et al. data

setwd('/nfs/users/nfs_a/aa16/Deconvolution/')
#setwd('/home/jovyan/Deconvolution/')

source('deconvolution.R')
import::from(myUtils, mapIdsMouse)

dataDirectory = '/lustre/scratch117/cellgen/team283/brainData'
#dataDirectory = '../data/'
resultsDirectory = 'results/'
figuresDirectory = 'figures/'

regionVector = c('Cerebellum', 'Frontal Cortex', 'Hippocampus', 'Striatum')

for ( i in 1:length(regionVector)){
  region = regionVector[i]
  # Prepare data:
  counts_sc = readRDS(paste(dataDirectory, 'Saunders/Saunders_Mouse_', region, '_counts.rds', sep = ''))
  coldata = readRDS(paste(dataDirectory, 'Saunders/Saunders_Mouse_', region, '_coldata.rds', sep = ''))
  rowdata = readRDS(paste(dataDirectory, 'Saunders/Saunders_Mouse_', region, '_rowdata.rds', sep = ''))
  reads_bulk = as.matrix(read.table(paste(dataDirectory, 'KptnMouse/', region, '_Kptn_4col_ReadsPerGene.txt', sep = '')))
  rownames(counts_sc) = mapIdsMouse(rownames(counts_sc), 'SYMBOL', 'ENSEMBL')
  
  # Get celltypes
  celltypes = as.character(coldata[,2])
  general_celltypes = unlist(lapply(celltypes, function(x) strsplit(x, split = '\\.')[[1]][1]))
  firstMarker = unlist(lapply(celltypes, function(x) strsplit(x, split = '\\.')[[1]][2]))
  general_celltypes[substring(firstMarker,1,5) == 'Slc17'] = paste('Glutamatergic', general_celltypes[substring(firstMarker,1,5) == 'Slc17'], sep = '')
  general_celltypes[substring(firstMarker,1,3) == 'Gad'] = paste('GABAergic', general_celltypes[substring(firstMarker,1,3) == 'Gad'], sep = '')
  
  # Run MuSiC:
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = general_celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'Deconvolution.rds', sep = ''))
  plotDeconvolution(results, save = FALSE)
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionSummary.csv', sep = ''))
}


