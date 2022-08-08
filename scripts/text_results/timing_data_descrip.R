# Descriptions of files used for resource usage tests

rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/prepped_infiles/")

complete_test.func <- read.table(file = "complete_test.func.tsv.gz",
                                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

complete_test.abun <- read.table(file = "complete_test.abun.tsv.gz",
                                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

dim(complete_test.func)
dim(complete_test.abun)

write.table(x = halved_func_prepped_files$func, file = "halved.func_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = halved_func_prepped_files$abun, file = "halved.func_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = halved_abun_prepped_files$func, file = "halved.abun_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = halved_abun_prepped_files$abun, file = "halved.abun_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = quartered_func_prepped_files$func, file = "quartered.func_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = quartered_func_prepped_files$abun, file = "quartered.func_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = quartered_abun_prepped_files$func, file = "quartered.abun_test.func.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = quartered_abun_prepped_files$abun, file = "quartered.abun_test.abun.tsv", col.names = NA, row.names = TRUE, sep = "\t", quote = FALSE)
