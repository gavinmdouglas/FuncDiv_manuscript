rm(list = ls(all.names = TRUE))

setwd("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/prepped_infiles")

complete_test.func <- read.table(file = "complete_test.func.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
complete_test.abun <- read.table(file = "complete_test.abun.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
dim(complete_test.func)
dim(complete_test.abun)

halved.func_test.func <- read.table(file = "halved.func_test.func.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
halved.func_test.abun <- read.table(file = "halved.func_test.abun.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
dim(halved.func_test.func)
dim(halved.func_test.abun)

halved.abun_test.func <- read.table(file = "halved.abun_test.func.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
halved.abun_test.abun <- read.table(file = "halved.abun_test.abun.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
dim(halved.abun_test.func)
dim(halved.abun_test.abun)


quartered.func_test.func <- read.table(file = "quartered.func_test.func.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
quartered.func_test.abun <- read.table(file = "quartered.func_test.abun.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
dim(quartered.func_test.func)
dim(quartered.func_test.abun)

quartered.abun_test.func <- read.table(file = "quartered.abun_test.func.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
quartered.abun_test.abun <- read.table(file = "quartered.abun_test.abun.tsv.gz", header = TRUE, row.names = 1, sep = "\t")
dim(quartered.abun_test.func)
dim(quartered.abun_test.abun)
