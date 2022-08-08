# Script to run entire alpha diversity metric workflow on specified function copy number and taxonomic abundances tables.
# Will run three metrics across all input functions / samples.

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(FuncDiv))

parser <- ArgumentParser()

parser$add_argument('func_abun', metavar = 'FUNC_TABLE', type = "character", nargs = 1,
                    help = 'Path to table containing function copy numbers across taxa.')

parser$add_argument('taxa_abun', metavar = 'ABUN_TABLE', type = "character", nargs = 1,
                    help = 'Path to table containing taxa abundances across samples.')

parser$add_argument('tree', metavar = 'TREE', type = "character", nargs = 1,
                    help = 'Path to tree, which is needed to compute Faith\'s PD.')

parser$add_argument('outfile', metavar = 'OUTPUT', type = "character", nargs = 1,
                    help = 'Path to output RDS file.')

parser$add_argument('--num_cores', metavar = 'Number of cores', type = "integer", nargs = 1,
                    default = 1,
                    help = 'Number of cores to run parallelizable steps.')

args <- parser$parse_args()

metrics_to_run <- c("richness", "gini_simpson_index", "faiths_pd")

func_tab <- read.table(file = args$func_abun, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
abun_tab <- read.table(file = args$taxa_abun, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
in_tree <- ape::read.tree(args$tree)

output <- alpha_div_contrib(metrics = metrics_to_run,
                            func_tab = func_tab,
                            abun_tab = abun_tab,
                            in_tree = in_tree,
                            ncores = args$num_cores)

saveRDS(object = output, file = args$outfile)
