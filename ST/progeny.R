library('progeny')

samps <- c('IR0h', 'IR6h', 'IR24h', 'IR72h', 'IR30D')

for (samp in samps){
  progeny_scores = progeny::progeny(expr = as.matrix(slide[['Spatial']]@data),
                                    scale=TRUE,
                                    organism='Mouse', top=500, perm=1, verbose=T)
  slide[['progeny']] = CreateAssayObject(counts = t(progeny_scores))
  saveRDS(paste0('data/', samp, '_progeny.rds'))
}
