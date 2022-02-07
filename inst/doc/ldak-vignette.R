## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = 'png',
  eval = F
)

## ----plink_qc-----------------------------------------------------------------
#  for(i in 1:22) {
#         bigsnpr::snp_plinkQC(plink.path = '/bin/plink', prefix.in = paste0('euro_only/chr', i),
#                            prefix.out = paste0('euro_only_qc/chr', i),
#                            geno = 0, maf = 0.01, hwe = 1e-10)
#  }
#  bigsnpr::snp_plinkQC(plink.path = '/bin/plink', prefix.in = 'euro_only/chrX',
#                          prefix.out = 'euro_only_qc/chrX',
#                          geno = 0, maf = 0.01, hwe = 1e-10)
#  
#  bigsnpr::snp_plinkQC(plink.path = '/bin/plink', prefix.in = 'euro_only/chrY',
#                          prefix.out = 'euro_only_qc/chrY',
#                          geno = 0, maf = 0.01, hwe = 1e-10)

## ----gwas_bim_join------------------------------------------------------------
#  library(data.table)
#  system("mkdir plinkRanges")
#  system("mkdir filtered")
#  system("mkdir ldak")
#  
#  # We assume this contains columns SNPID, CHR19 and BP19 (so-called because of hg19), REF, and ALT
#  # note that gwas_dat needs to be a data.table
#  gwas_dat <- fread('gwas_sum_stats.tsv.gz', sep = '\t', header = T)
#  
#  # Iterating over chromosomal bim files
#  for(i in 1:22) {
#    # bim files have no header
#    bim_dat <- fread(sprintf('euro_only_qc/chr%d.bim', i), sep = '\t', header = F, col.names = c('Chr', 'ID', 'Cm', 'BP19', 'A1', 'A2') )
#  
#    bim_join <- merge(bim_dat, gwas_dat[CHR19 == i], by.x = 'BP19', by.y = 'BP19')
#  
#    # Make sure alleles match, although for two-sided association p-values we don't care whether ref/alt is reversed
#    bim_join <- bim_join[(REF == A1 & ALT == A2) | (REF == A2 & ALT == A1)]
#  
#    bim_join <- bim_join[, .(Chr, BP19, BP19, SNPID)]
#  
#    if(any(duplicated(bim_join, by='BP19'))) {
#      warning(sprintf('%d duplicates removed from output', sum(duplicated(bim_join, by = 'BP19'))))
#    }
#  
#    # Remove duplicates
#    bim_join <- unique(bim_join, by='BP19')
#  
#    fwrite(bim_join, file = sprintf('plinkRanges/chr%d.tsv', i), row.names = F, sep = '\t', col.names = F, quote = F)
#  }

## ----merge_bim_combined_weights-----------------------------------------------
#  library(data.table)
#  
#  bim_dat <- fread('filtered/chr_all.bim', sep = '\t', header = F, col.names = c('Chr', 'ID', 'Cm', 'BP19', 'A1', 'A2'))
#  
#  weights_dat <- fread('ldak/combined_weights.all', sep = ' ', header = T)
#  
#  # Drop rows with a missing ID or weight value
#  weights_dat <- na.omit(weights_dat, cols = c('Predictor', 'Weight'))
#  
#  join_dat <- merge(weights_dat, bim_dat[, .(ID, Chr, BP19, A1, A2)], all.x = T, by.x = c('Predictor', 'Chr'), by.y = c('ID', 'Chr'), sort = F)
#  
#  fwrite(join_dat, file = 'ldak/combined_weights_meta.all', sep = ' ', col.names = T, row.names = F, quote = F)

## -----------------------------------------------------------------------------
#  weights_dat <- fread('ldak/combined_weights_meta.all', sep = ' ', header = T, select = c('Predictor', 'Weight', 'Chr', 'BP19', 'A1', 'A2'))
#  
#  # We assume this contains columns SNPID, CHR19 and BP19 (so-called because of hg19), REF, and ALT
#  gwas_dat <- fread('gwas_sum_stats.tsv.gz', sep = '\t', header = T)
#  
#  gwas_dat <- merge(gwas_dat, weights_dat, by.x = c('CHR19', 'BP19'), by.y = c('Chr', 'BP19'), sort = F)
#  
#  # Drop rows where the ref/alt allele pairing differs from that already present
#  gwas_dat <- gwas_dat[((REF == A1 & ALT == A2) | (REF == A2 & ALT == A1))]
#  
#  # Drop now-redundant allele columns
#  gwas_dat[, c('A1', 'A2') := NULL ]

