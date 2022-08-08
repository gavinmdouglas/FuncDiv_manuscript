# Run RF models on each input table.
# Also run alternative versions of RF models with relative abundance variables included too.

rm(list = ls(all.names = TRUE))

library(FuncDiv)
library(ranger)

clr_transform_by_col <- function(in_df) {
  # Perform CLR transformation on table assuming samples are columns.
  return(data.frame(apply(in_df, 2, function(x){ log(x) - mean(log(x)) }), check.names = FALSE))
}

metrics_to_ignore <- c("menhinicks_richness", "mcintoshs_evenness", "mcintoshs_dominance",
                       "margalefs_richness", "fishers_alpha", "faiths_pd")
metrics <- names(FuncDiv_alpha_metrics)[which(! names(FuncDiv_alpha_metrics) %in% metrics_to_ignore)]

for (disease in c("CRC", "IBD", "STH")) {

  metadata_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/", disease, "_metadata.tsv.gz", sep = "")
  metadata <- read.table(file = metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(metadata) <- metadata$sample_id
  
  tables <- list()
  for (m in metrics) {
    
    contrib_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/",
                          disease, "_pathway_alpha_contrib/", m, ".tsv.gz", sep = "")
    
    metric_raw <- read.table(file = contrib_path, header = TRUE, row.names = 1)
    
    for (s in colnames(metric_raw)) {
      metric_raw[which(abs(metric_raw[, s]) >= "inf"), s] <- min(metric_raw[, s])
    }
    
    tables[[m]] <- data.frame(scale(metric_raw, center = TRUE, scale = TRUE),
                              check.names = FALSE)
  }
  
  pathway_genera_abun_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/",
                                    disease, "_summed_genus_rel_levels_by_pathway.tsv.gz", sep = "")
  
  tables[["contributing_genus_clr"]] <- read.table(file = pathway_genera_abun_path, header = TRUE, row.names = 1)
  tables[["contributing_genus_clr"]] <- data.frame(sweep(x = tables[["contributing_genus_clr"]],
                                                         MARGIN = 2,
                                                         STATS = colSums(tables[["contributing_genus_clr"]]),
                                                         FUN = '/')) * 100
  tables[["contributing_genus_clr"]] <- tables[["contributing_genus_clr"]] + 0.1
  tables[["contributing_genus_clr"]] <- clr_transform_by_col(tables[["contributing_genus_clr"]])
  
  
  pathway_levels_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/prepped/",
                          disease, "_pathway_community_rel_levels.tsv.gz", sep = "")
  
  tables[["pathway_levels_clr"]] <- read.table(file = pathway_levels_path, header = TRUE, row.names = 1)
  tables[["pathway_levels_clr"]] <- data.frame(sweep(x = tables[["pathway_levels_clr"]],
                                                         MARGIN = 2,
                                                         STATS = colSums(tables[["pathway_levels_clr"]]),
                                                         FUN = '/')) * 100
  tables[["pathway_levels_clr"]] <- tables[["pathway_levels_clr"]] + 0.1
  tables[["pathway_levels_clr"]] <- clr_transform_by_col(tables[["pathway_levels_clr"]])

  
  # Need to transpose tables and add "disease_state" column.
  RF_input <- list()
  
  intersecting_samples <- colnames(tables$richness)[which(colnames(tables$richness) %in% rownames(metadata))]
  
  for (category in names(tables)) {
    RF_input[[category]] <- data.frame(t(tables[[category]]), check.names = FALSE)
    
    missing_samples <- which(! intersecting_samples %in% rownames(RF_input[[category]]) )
    
    if (length(missing_samples) > 0) { stop("Some samples missing") }
    
    RF_input[[category]] <- RF_input[[category]][intersecting_samples, ]
    RF_input[[category]]$disease_state <- as.factor(metadata[intersecting_samples, "disease"])
  }
  

  # Also make additional RF input tables for contrib div alpha metrics with relabun variables added too.
  genus_table <- RF_input$contributing_genus_clr
  genus_table <- genus_table[, -which(colnames(genus_table) == "disease_state")]
  colnames(genus_table) <- paste(colnames(genus_table), ".relabun", sep = "")
  for (m in metrics) {
    category_name <- paste(m, "w_relabun", sep = ".")
    RF_input[[category_name]] <- cbind(RF_input[[m]], genus_table)
  }
  
  
  # Run RF on all of these input tables.
  major_class_prop <- max(table(RF_input$richness$disease_state)) / nrow(RF_input$richness)
  oob_accuracy <- c(major_class_prop)
  
  for (category in names(RF_input)) {
    
    RF_out <- ranger(importance = "permutation",
                     data = RF_input[[category]],
                     num.trees = 10000,
                     dependent.variable.name = "disease_state")
    
    oob_accuracy <- c(oob_accuracy, 1 - RF_out$prediction.error)
    
    
    RDS_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/RF_RDS/",
                      disease, ".", category, "_RF.rds", sep = "")
    saveRDS(object = RF_out, file = RDS_path)
    
    
    varImp_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/varImp/",
                         disease, ".", category, ".tsv", sep = "")
    varImp_df <- data.frame(pathway = names(RF_out$variable.importance),
                            permutation_varImp = RF_out$variable.importance,
                            stringsAsFactors = FALSE)
    write.table(x = varImp_df, file = varImp_path, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
    
  }
  
  oob_path <- paste("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/output/",
                       disease, "_oob_accuracy.tsv", sep = "")
  oob_accuracy_df <- data.frame(category = c("major_class_prop", names(RF_input)), oob_accuracy = oob_accuracy, stringsAsFactors = FALSE)
  write.table(x = oob_accuracy_df, file = oob_path, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

}
