### Deconvolution Pipeline Functions ###

require(MuSiC)
require(Biobase)
library(xbioc)
library(ggpubr)

# count matrices need to be in matrix format and include genes in rownames and
# sample names in colnames

getDeconvolution = function(counts_bulk, counts_sc, clusters,  return = TRUE, save = FALSE, file = 'Deconvolution.rds')
{
  # Prepare single cell data:
  row_data = data.frame(genes = rownames(counts_sc))
  column_data = data.frame(samples = colnames(counts_sc), clusters = clusters)
  rownames(column_data) = colnames(counts_sc)
  phenoData = new('AnnotatedDataFrame', data = column_data)
  scSet = ExpressionSet(assayData = counts_sc, phenoData = phenoData)
  # Prepare bulk data:
  bulkSet = ExpressionSet(assayData = counts_bulk)
  results = music_prop(bulkSet, scSet, clusters = 'clusters', samples='samples')
  if (save){
    saveRDS(results, file = file)
  }
  if (return){
    return(results) 
  }
}

plotDeconvolution = function(results, save = TRUE, file = 'deconvolution.pdf', width = 10, height = 10){
  tab = results[[1]]
  binary = (colSums(tab) > 0)
  tab = tab[,binary]
  newTab = data.frame(unlist(lapply(1:dim(tab)[2], function(x) rep(colnames(tab)[x],dim(tab)[1]))),
                      rep(rownames(tab),dim(tab)[2]),
                      as.vector(tab))
  colnames(newTab) = c('cell_cluster', 'genotype', 'proportion')
  newTab[,2] = substring(newTab[,2], 1,2)
  p <- ggboxplot(as.data.frame(newTab), x = "genotype", y = "proportion",
                 color = "genotype", palette = "jco",
                 add = "jitter",
                 facet.by = "cell_cluster", short.panel.labs = FALSE)
  if (save){
    pdf(file = file, width = width, height = height)
    p + stat_compare_means(aes(group = genotype), label = "p.signif")#, method = 'wilcox.test')
    dev.off()
  }else{
    p + stat_compare_means(aes(group = genotype), label = "p.signif")#, method = 'wilcox.test')
  }
}

# Write Deconvolutionr results in a table with p-value
summarizeDeconvolution = function(results, return = TRUE, save = FALSE, file = 'DeconvolutionSummary.csv'){
  tab = results[[1]]
  pValue = unlist(lapply(1:dim(tab)[2], function(x) t.test(results[[1]][1:5,x], results[[1]][6:11,x])$p.value))
  mean_wt = unlist(lapply(1:dim(tab)[2], function(x) mean(results[[1]][1:5,x])))
  mean_hom = unlist(lapply(1:dim(tab)[2], function(x) mean(results[[1]][6:11,x])))
  tab = rbind(mean_wt, mean_hom, pValue, tab)
  # Sort by p-value:
  tab = tab[, order(tab[3,])]
  tab = round(tab,4)
  write.csv(tab, file = file, quote = FALSE, row.names = TRUE)
}




