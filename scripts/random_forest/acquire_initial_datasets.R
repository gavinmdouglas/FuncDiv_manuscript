rm(list = ls(all.names = TRUE))

library(curatedMetagenomicData)
library(mia)
library(phyloseq)
library(stringr)
library(reshape2)


prep_curatedMetagenomicData_pathways_and_relabun <- function(pathway_string, relabun_string) {
  
  pathway_raw <- curatedMetagenomicData(pathway_string, dryrun = FALSE)
  relabun_raw <- curatedMetagenomicData(relabun_string, dryrun = FALSE)
  
  pathway_raw_phyloseq <- mia::makePhyloseqFromTreeSummarizedExperiment(pathway_raw[[pathway_string]],
                                                                        abund_values = "pathway_abundance")
  
  relabun_raw_phyloseq <- mia::makePhyloseqFromTreeSummarizedExperiment(relabun_raw[[relabun_string]],
                                                                        abund_values = "relative_abundance")
  
  relabun_tab <- data.frame(phyloseq::otu_table(relabun_raw_phyloseq))
  
  strat_pathway_tab <- data.frame(phyloseq::otu_table(pathway_raw_phyloseq))
  
  strat_pathway_tab <- strat_pathway_tab[-which(rownames(strat_pathway_tab) == "UNMAPPED"), ]
  strat_pathway_tab <- strat_pathway_tab[-grep("UNINTEGRATED", rownames(strat_pathway_tab)), ]
  
  unstrat_pathway_tab <- strat_pathway_tab[-grep("\\|", rownames(strat_pathway_tab)), ]
  strat_pathway_tab <- strat_pathway_tab[grep("\\|", rownames(strat_pathway_tab)), ]
  
  pathway_and_taxa <- str_split_fixed(rownames(strat_pathway_tab), "\\|", 2)
  
  strat_pathway_tab$func <- pathway_and_taxa[, 1]
  strat_pathway_tab$taxon <- pathway_and_taxa[, 2]
  
  strat_pathway_tab_melt <- melt(data = strat_pathway_tab,
                                 ids = c("func", "taxon"),
                                 variable.name = "sample",
                                 value.name = "relabun")
  
  return(list(strat=strat_pathway_tab_melt,
              unstrat=unstrat_pathway_tab,
              relabun=relabun_tab))
}


# CRC dataset
YachidaS_2019_subset <- sampleMetadata[which(sampleMetadata$study_name == "YachidaS_2019"), ]
YachidaS_2019_subset <- YachidaS_2019_subset[which(YachidaS_2019_subset$disease %in% c("CRC", "healthy")), ]

### This command returns all versions of the pathway abundance table that are available.  
curatedMetagenomicData("YachidaS_2019.pathway_abundance", dryrun = TRUE)

curatedMetagenomicData("YachidaS_2019.relative_abundance", dryrun = TRUE)



YachidaS_2019_output_tables <- prep_curatedMetagenomicData_pathways_and_relabun(pathway_string = "2021-10-14.YachidaS_2019.pathway_abundance", 
                                                                                relabun_string = "2021-10-14.YachidaS_2019.relative_abundance")

YachidaS_2019_output_tables[["metadata"]] = YachidaS_2019_subset


# IBD dataset
NielsenHB_2014_subset <- sampleMetadata[which(sampleMetadata$study_name == "NielsenHB_2014"), ]
NielsenHB_2014_subset <- NielsenHB_2014_subset[which(NielsenHB_2014_subset$days_from_first_collection == 0), ]

curatedMetagenomicData("NielsenHB_2014.pathway_abundance", dryrun = TRUE)

NielsenHB_2014_output_tables <- prep_curatedMetagenomicData_pathways_and_relabun(pathway_string = "2021-03-31.NielsenHB_2014.pathway_abundance",
                                                                                 relabun_string = "2021-03-31.NielsenHB_2014.relative_abundance")

NielsenHB_2014_output_tables[["metadata"]] = NielsenHB_2014_subset


# Soil-transmitted helminth (STH) infection
RubelMA_2020_subset <- sampleMetadata[which(sampleMetadata$study_name == "RubelMA_2020"), ]

curatedMetagenomicData("RubelMA_2020.pathway_abundance", dryrun = TRUE)

RubelMA_2020_output_tables <- prep_curatedMetagenomicData_pathways_and_relabun(pathway_string = "2021-10-14.RubelMA_2020.pathway_abundance",
                                                                               relabun_string = "2021-10-14.RubelMA_2020.relative_abundance")

RubelMA_2020_output_tables[["metadata"]] = RubelMA_2020_subset


saveRDS(object = YachidaS_2019_output_tables,
        file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/CRC_output_tables.rds")

saveRDS(object = NielsenHB_2014_output_tables,
        file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/IBD_output_tables.rds")

saveRDS(object = RubelMA_2020_output_tables,
        file = "/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/random_forest/input/raw_download/STH_output_tables.rds")
