# Swap in genus abundances into pathway tables (so use that rather than the inferred pathway abundances).
# Also save the community-wide pathway abundances (i.e., the unstrat tables)

rm(list = ls(all.names = TRUE))

prep_tables <- function(RDS_path) {

  output_tables <- readRDS(file = RDS_path)
  
  relabun_genus <- output_tables$relabun
  
  relabun_genus$genus <- rownames(relabun_genus)
  
  relabun_genus$genus <- gsub("\\|s__.*$", "", relabun_genus$genus)
  relabun_genus$genus <- gsub("^.*\\|g__", "g__", relabun_genus$genus)
  
  relabun_genus <- aggregate(formula = . ~ genus,
                             data = relabun_genus,
                             FUN = sum)
  
  rownames(relabun_genus) <- relabun_genus$genus
  relabun_genus <- relabun_genus[, -which(colnames(relabun_genus) == "genus")]
  relabun_genus <- data.frame(sweep(x = relabun_genus, MARGIN = 2, STATS = colSums(relabun_genus), FUN = '/')) * 100
  
  # Regroup to genera and remove unclassified taxa.
  pathway_strat_prepped <- output_tables$strat
  pathway_strat_prepped$sample <- as.character(pathway_strat_prepped$sample)
  pathway_strat_prepped <- pathway_strat_prepped[-which(pathway_strat_prepped$taxon == "unclassified"), ]
  pathway_strat_prepped$genus <- pathway_strat_prepped$taxon
  pathway_strat_prepped$genus <- gsub(".s__.*$", "", pathway_strat_prepped$genus)
  
  # Sanity check that all genera ids in pathway table are actually present.
  missing_genera <- pathway_strat_prepped$genus[which(! pathway_strat_prepped$genus %in% rownames(relabun_genus))]
  if (length(missing_genera) > 0) {
    stop("Some genera names in pathway table not found in relative abundance table.")
  }
  
  # Remove taxon and relabun columns.
  pathway_strat_prepped <- pathway_strat_prepped[, -which(colnames(pathway_strat_prepped) %in% c("taxon", "relabun"))]
  
  # Simplify function ids.
  pathway_strat_prepped$func <- gsub(pattern = ":.*$", replacement = "", x = pathway_strat_prepped$func)
  
  # Remove duplicate lines.
  pathway_strat_prepped <- pathway_strat_prepped[-which(duplicated(pathway_strat_prepped)), ]
  
  # Add in genera abundances.
  # Do this in a roundabout way to different lists, because otherwise it took too much memory and/or time.
  new_subsets <- list()
  for (s in unique(pathway_strat_prepped$sample)) {
    new_subsets[[s]] <- pathway_strat_prepped[which(pathway_strat_prepped$sample == s), ]
    new_subsets[[s]]$genus_relabun <- relabun_genus[new_subsets[[s]]$genus, s]
  }
  
  pathway_strat_prepped <- do.call(rbind, new_subsets)
  
  rownames(pathway_strat_prepped) <- NULL
  
  pathabun_unstrat <- output_tables$unstrat
  rownames(pathabun_unstrat) <- gsub(pattern = ":.*$", replacement = "", x = rownames(pathabun_unstrat))
  pathabun_unstrat <- data.frame(sweep(x = pathabun_unstrat, MARGIN = 2, STATS = colSums(pathabun_unstrat), FUN = '/')) * 100
  
  return(list(path_by_genus.abun_strat = pathway_strat_prepped,
              pathabun_unstrat = pathabun_unstrat,
              metadata = output_tables$metadata))
  
}

CRC_prepped <- prep_tables("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/CRC_output_tables.rds")
IBD_prepped <- prep_tables("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/IBD_output_tables.rds")
STH_prepped <- prep_tables("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/STH_output_tables.rds")

write.table(x = CRC_prepped$path_by_genus.abun_strat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_pathway_genera_abun.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = CRC_prepped$pathabun_unstrat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_pathway_community_rel_levels.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = CRC_prepped$metadata,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/CRC_metadata.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



write.table(x = IBD_prepped$path_by_genus.abun_strat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_pathway_genera_abun.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = IBD_prepped$pathabun_unstrat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_pathway_community_rel_levels.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = IBD_prepped$metadata,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/IBD_metadata.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



write.table(x = STH_prepped$path_by_genus.abun_strat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_pathway_genera_abun.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(x = STH_prepped$pathabun_unstrat,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_pathway_community_rel_levels.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = STH_prepped$metadata,
            file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/STH_metadata.tsv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

