# Quick commands to get #'s to report in text.

tmp <- read.table("../../input/prepped/STH_pathway_genera_abun.tsv.gz",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

tmp2 <- read.table("../../input/prepped/STH_pathway_community_rel_levels.tsv.gz",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

tmp3 <- read.table("../../input/prepped/STH_metadata.tsv.gz",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(tmp3) <- tmp3$sample_id

intersecting_samples <- colnames(tmp2)[which(colnames(tmp2) %in% rownames(tmp3))]

tmp3 <- tmp3[intersecting_samples, ]
