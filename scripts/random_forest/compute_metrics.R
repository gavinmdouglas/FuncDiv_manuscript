# Compute contributional alpha diversity metrics for all three datasets in order to run random forest.
# Also get community-wide function relative abundance to include as well.
# Note that tables containing the relative pathway levels contributed by each taxon were already output
# (which are a distinct measure from the taxon's relative abundance!)

rm(list = ls(all.names = TRUE))

library(FuncDiv)

metrics_to_ignore <- c("menhinicks_richness", "mcintoshs_evenness", "mcintoshs_dominance",
                       "margalefs_richness", "fishers_alpha", "faiths_pd")

metrics_to_run <- names(FuncDiv_alpha_metrics)[which(! names(FuncDiv_alpha_metrics) %in% metrics_to_ignore)]


CRC_pathway_genus_abun <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_pathway_genera_abun.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

CRC_pathway_genus_alpha <- alpha_div_contrib(metrics = metrics_to_run,
                                             contrib_tab = CRC_pathway_genus_abun,
                                             replace_NA = TRUE,
                                             ncores = 10,
                                             samp_colname = "sample",
                                             func_colname = "func",
                                             abun_colname = "genus_relabun",
                                             taxon_colname = "genus")

for (m in names(CRC_pathway_genus_alpha)) {
  
  outfile <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_pathway_alpha_contrib/", m, ".tsv", sep = "")
  
  write.table(x = CRC_pathway_genus_alpha[[m]],
              file = outfile,
              quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
  
}

CRC_pathway_summed_abun <- CRC_pathway_genus_abun
CRC_pathway_summed_abun <- CRC_pathway_summed_abun[, -which(colnames(CRC_pathway_summed_abun) == "genus")]
CRC_pathway_summed_abun <- aggregate(formula = genus_relabun ~ ., data = CRC_pathway_summed_abun, FUN = sum)
CRC_pathway_summed_abun <- reshape2::dcast(data = CRC_pathway_summed_abun, formula = func ~ sample, value.var = "genus_relabun")
rownames(CRC_pathway_summed_abun) <- CRC_pathway_summed_abun$func
CRC_pathway_summed_abun <- CRC_pathway_summed_abun[, -which(colnames(CRC_pathway_summed_abun) == "func")]

write.table(x = CRC_pathway_summed_abun,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_summed_genus_rel_levels_by_pathway.tsv",
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)




IBD_pathway_genus_abun <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_pathway_genera_abun.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

IBD_pathway_genus_alpha <- alpha_div_contrib(metrics = metrics_to_run,
                                             contrib_tab = IBD_pathway_genus_abun,
                                             replace_NA = TRUE,
                                             ncores = 10,
                                             samp_colname = "sample",
                                             func_colname = "func",
                                             abun_colname = "genus_relabun",
                                             taxon_colname = "genus")

for (m in names(IBD_pathway_genus_alpha)) {
  
  outfile <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_pathway_alpha_contrib/", m, ".tsv", sep = "")
  
  write.table(x = IBD_pathway_genus_alpha[[m]],
              file = outfile,
              quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
  
}

IBD_pathway_summed_abun <- IBD_pathway_genus_abun
IBD_pathway_summed_abun <- IBD_pathway_summed_abun[, -which(colnames(IBD_pathway_summed_abun) == "genus")]
IBD_pathway_summed_abun <- aggregate(formula = genus_relabun ~ ., data = IBD_pathway_summed_abun, FUN = sum)
IBD_pathway_summed_abun <- reshape2::dcast(data = IBD_pathway_summed_abun, formula = func ~ sample, value.var = "genus_relabun")
rownames(IBD_pathway_summed_abun) <- IBD_pathway_summed_abun$func
IBD_pathway_summed_abun <- IBD_pathway_summed_abun[, -which(colnames(IBD_pathway_summed_abun) == "func")]

write.table(x = IBD_pathway_summed_abun,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_summed_genus_rel_levels_by_pathway.tsv",
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)




STH_pathway_genus_abun <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_pathway_genera_abun.tsv.gz",
                                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

STH_pathway_genus_alpha <- alpha_div_contrib(metrics = metrics_to_run,
                                             contrib_tab = STH_pathway_genus_abun,
                                             replace_NA = TRUE,
                                             ncores = 10,
                                             samp_colname = "sample",
                                             func_colname = "func",
                                             abun_colname = "genus_relabun",
                                             taxon_colname = "genus")

for (m in names(STH_pathway_genus_alpha)) {
  
  outfile <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_pathway_alpha_contrib/", m, ".tsv", sep = "")
  
  write.table(x = STH_pathway_genus_alpha[[m]],
              file = outfile,
              quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
  
}

STH_pathway_summed_abun <- STH_pathway_genus_abun
STH_pathway_summed_abun <- STH_pathway_summed_abun[, -which(colnames(STH_pathway_summed_abun) == "genus")]
STH_pathway_summed_abun <- aggregate(formula = genus_relabun ~ ., data = STH_pathway_summed_abun, FUN = sum)
STH_pathway_summed_abun <- reshape2::dcast(data = STH_pathway_summed_abun, formula = func ~ sample, value.var = "genus_relabun")
rownames(STH_pathway_summed_abun) <- STH_pathway_summed_abun$func
STH_pathway_summed_abun <- STH_pathway_summed_abun[, -which(colnames(STH_pathway_summed_abun) == "func")]

write.table(x = STH_pathway_summed_abun,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_summed_genus_rel_levels_by_pathway.tsv",
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
