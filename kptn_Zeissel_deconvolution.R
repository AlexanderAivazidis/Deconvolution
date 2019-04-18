### Kptn bulk data deconvolution using Zeissel et al. data

### Kptn bulk data deconvolution using Saunders et al. data

farm = FALSE

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

for (i in 1:length(regionVector)){
  region = regionVector[i]
  print(region)
  # Prepare data:
  counts_sc = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_counts.rds', sep = ''))
  coldata = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_coldata.rds', sep = ''))
  rowdata = readRDS(paste(dataDirectory, 'Zeissel/Zeissel_Mouse_', region, '_rowdata.rds', sep = ''))
  reads_bulk = as.matrix(read.table(paste(dataDirectory, 'KptnMouse/', region, '_Kptn_4col_ReadsPerGene.txt', sep = '')))
  rownames(counts_sc) = rowdata[,'accession']
  colnames(counts_sc) = paste(coldata[,'CellID'], coldata[,'MitoRiboRatio'])
  
  # Get celltypes
  celltypes = as.character(coldata[,'CombinedClusterName'])

  # Run MuSiC:
  results = getDeconvolution(counts_bulk = reads_bulk , counts_sc = counts_sc, clusters = celltypes, return = TRUE,
                             save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselSpecificCelltypes.rds', sep = ''))
  summarizeDeconvolution(results, return = FALSE, save = TRUE, file = paste(resultsDirectory, region, 'DeconvolutionZeisselSpecificCelltypesSummary.csv', sep = ''))
  plotDeconvolution(results, save = TRUE, file = paste(figuresDirectory, region, 'DeconvolutionZeisselSpecificCelltypes.pdf', sep = ''),
                    width = 10, height = 10)
}

