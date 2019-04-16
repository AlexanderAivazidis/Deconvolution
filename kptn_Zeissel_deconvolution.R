### Kptn bulk data deconvolution using Zeissel et al. data

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
  counts_sc = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_counts.rds', sep = ''))
  coldata = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_coldata.rds', sep = ''))
  rowdata = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_rowdata.rds', sep = ''))
  reads_bulk = as.matrix(read.table(paste(dataDirectory, 'KptnMouse/', region, '_Kptn_4col_ReadsPerGene.txt', sep = '')))
  rownames(counts_sc) = mapIdsMouse(rownames(counts_sc), 'SYMBOL', 'ENSEMBL')
  
  # Get celltypes
  celltypes = as.character(coldata[,2])
  general_celltypes = unlist(lapply(celltypes, function(x) strsplit(x, split = '\\.')[[1]][1]))

  # Run MuSiC:
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselSpecificCelltypes.rds', sep = ''))
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselSpecificCelltypesSummary.csv', sep = ''))
  plotDeconvolution(results, save = TRUE, file = paste(figuresDirectory, region, 'DeconvolutionZeisselSpecificCelltypes.pdf', sep = ''),
                    width = 10, height = 10)
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = general_celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselGeneralCelltypes.rds', sep = ''))
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselGeneralCelltypesSummary.csv', sep = ''))
  plotDeconvolution(results, save = TRUE, file = paste(figuresDirectory, region, 'DeconvolutionZeisselGeneralCelltypes.pdf', sep = ''),
                    width = 10, height = 10)
}


