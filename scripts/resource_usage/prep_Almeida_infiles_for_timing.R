# Prep infiles to be used for assessing the runtime and memory usage.
# This will include filtering out rare features and also possible sub-sampling tables 
# to see how the pipelines work with fewer features.

rm(list = ls(all.names = TRUE))

library(FuncDiv)

full_abun <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/MAG_abun/bwa_depth_min25coverage.tsv.gz",
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

full_ko <- read.table("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/MAG_annot/kegg_ortholog_annot.tsv.gz",
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

full_subsetted_tables <- subset_func_and_abun_tables(func_table = full_ko,
                                                     abun_table = full_abun)

# Subsample to ~500 samples and ~2000 functions.
set.seed(14514)
test_ko <- full_subsetted_tables$func[sample(1:nrow(full_subsetted_tables$func), size = 2050), ]
test_abun <- full_subsetted_tables$abun[, sample(1:ncol(full_subsetted_tables$abun), size = 525)]

test_subsetted_tables <- subset_func_and_abun_tables(func_table = test_ko,
                                                     abun_table = test_abun)


halved_func_tab <- test_subsetted_tables$func[sample(1:nrow(test_subsetted_tables$func), size = round(0.5 * nrow(test_subsetted_tables$func))), ]
quartered_func_tab <- halved_func_tab[sample(1:nrow(halved_func_tab), size = round(0.5 * nrow(halved_func_tab))), ]

halved_abun_tab <- test_subsetted_tables$abun[, sample(1:ncol(test_subsetted_tables$abun), size = round(0.5 * ncol(test_subsetted_tables$abun)))]
quartered_abun_tab <- halved_abun_tab[, sample(1:ncol(halved_abun_tab), size = round(0.5 * ncol(halved_abun_tab)))]

halved_func_prepped_files <- subset_func_and_abun_tables(func_table = halved_func_tab, abun_table = test_subsetted_tables$abun)
halved_abun_prepped_files <- subset_func_and_abun_tables(func_table = test_subsetted_tables$func, abun_table = halved_abun_tab)

quartered_func_prepped_files <- subset_func_and_abun_tables(func_table = quartered_func_tab, abun_table = test_subsetted_tables$abun)
quartered_abun_prepped_files <- subset_func_and_abun_tables(func_table = test_subsetted_tables$func, abun_table = quartered_abun_tab)

setwd("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/prepped_infiles")

write.table(x = test_subsetted_tables$func, file = "complete_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = test_subsetted_tables$abun, file = "complete_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = halved_func_prepped_files$func, file = "halved.func_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = halved_func_prepped_files$abun, file = "halved.func_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = halved_abun_prepped_files$func, file = "halved.abun_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = halved_abun_prepped_files$abun, file = "halved.abun_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = quartered_func_prepped_files$func, file = "quartered.func_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = quartered_func_prepped_files$abun, file = "quartered.func_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = quartered_abun_prepped_files$func, file = "quartered.abun_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = quartered_abun_prepped_files$abun, file = "quartered.abun_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
