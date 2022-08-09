rm(list = ls(all.names = TRUE))

oob_out <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/STH_oob_accuracy.tsv",
                      header = TRUE, sep = "\t", row.names = 1)

genera_based_acc <- oob_out["contributing_genus_clr", 1]


metrics_to_ignore <- c("menhinicks_richness", "mcintoshs_evenness", "mcintoshs_dominance",
                       "margalefs_richness", "fishers_alpha", "faiths_pd")
metrics <- names(FuncDiv_alpha_metrics)[which(! names(FuncDiv_alpha_metrics) %in% metrics_to_ignore)]

max(oob_out[metrics, 1])
mean(oob_out[metrics, 1])
sd(oob_out[metrics, 1])

mean(oob_out[paste(metrics, ".w_relabun", sep = ""), 1])
sd(oob_out[paste(metrics, ".w_relabun", sep = ""), 1])
