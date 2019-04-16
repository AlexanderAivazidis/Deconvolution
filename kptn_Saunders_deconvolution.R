### Kptn bulk data deconvolution using Saunders et al. data

farm = FALSE
import::from(myUtils, mapIdsMouse)

if (farm == TRUE){
  setwd('/nfs/users/nfs_a/aa16/Deconvolution/')
}else{
  setwd('/home/jovyan/Deconvolution/')
}

source('deconvolution.R')

if (farm == TRUE){
  dataDirectory = '/lustre/scratch117/cellgen/team283/brainData/' 
}else{
  dataDirectory = '../data/'
}

resultsDirectory = 'results/'
figuresDirectory = 'figures/'

regionVector = c('Cerebellum', 'FrontalCortex', 'Hippocampus', 'Striatum')

for ( i in 1:length(regionVector)){
  region = regionVector[i]
  print(region)
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
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionSpecificCelltypes.rds', sep = ''))
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionSpecificCelltypesSummary.csv', sep = ''))
  plotDeconvolution(results, save = TRUE, file = paste(figuresDirectory, region, 'DeconvolutionSpecificCelltypes.pdf', sep = ''),
                               width = 10, height = 10)
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = general_celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionGeneralCelltypes.rds', sep = ''))
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionGeneralCelltypesSummary.csv', sep = ''))
  plotDeconvolution(results, save = TRUE, file = paste(figuresDirectory, region, 'DeconvolutionGeneralCelltypes.pdf', sep = ''),
                               width = 10, height = 10)
}


