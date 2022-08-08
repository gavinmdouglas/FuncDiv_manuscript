rm(list = ls(all.names = TRUE))

library(ComplexUpset)
library(ggplot2)


setwd("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/varImp/")

metrics_to_ignore <- c("menhinicks_richness", "mcintoshs_evenness", "mcintoshs_dominance",
                       "margalefs_richness", "fishers_alpha", "faiths_pd")
metrics <- names(FuncDiv_alpha_metrics)[which(! names(FuncDiv_alpha_metrics) %in% metrics_to_ignore)]

metrics <- c(metrics, "contributing_genus_clr")

disease <- "STH"

disease_raw_varImp <- list()

for (m in metrics) {
  
  in_varImp <- read.table(paste(disease, m, "tsv", sep = "."), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  colnames(in_varImp) <- m
  disease_raw_varImp[[m]] <- in_varImp
  
}

disease_varImp <- do.call(cbind, disease_raw_varImp)

disease_varImp_rank <- data.frame(apply(disease_varImp, 2, function(x) { rank(plyr::desc(x)) }), check.names = FALSE)

# Classify features as in top 20 or not.
disease_varImp_rank_top20 <- disease_varImp_rank
disease_varImp_rank_top20[disease_varImp_rank_top20 > 20] <- 0
disease_varImp_rank_top20[disease_varImp_rank_top20 > 0] <- 1

# Remove features not found in top 20 of any metric.
disease_varImp_rank_top20 <- disease_varImp_rank_top20[which(rowSums(disease_varImp_rank_top20) > 0), ]

upset_out <- upset(disease_varImp_rank_top20,
                    intersect = colnames(disease_varImp_rank_top20),
                    min_size = 2,
                    set_sizes = FALSE)

upset_outfile <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/raw_upset_plots/",
                       disease,
                       "_top20_upset_raw.pdf",
                       sep = "")

ggsave(filename = upset_outfile,
       device = "pdf",
       plot = upset_out,
       width = 8, height = 5)

