# Script to run beta diversity metric workflow on specified function copy number and taxonomic abundances tables.
# Will run two metrics across all input functions / samples.
# Note that the

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(FuncDiv))

parser <- ArgumentParser()

parser$add_argument('func_abun', metavar = 'FUNC_TABLE', type = "character", nargs = 1,
                    help = 'Path to table containing function copy numbers across taxa.')

parser$add_argument('taxa_abun', metavar = 'ABUN_TABLE', type = "character", nargs = 1,
                    help = 'Path to table containing taxa abundances across samples.')

parser$add_argument('in_tree', metavar = 'TREE', type = "character", nargs = 1,
                    help = 'Path to newick treefile.')

parser$add_argument('outdir', metavar = 'OUTPUT', type = "character", nargs = 1,
                    help = 'Path to output directory to make')

parser$add_argument('--num_cores', metavar = 'Number of cores', type = "integer", nargs = 1,
                    default = 1,
                    help = 'Number of cores to run parallelizable steps.')

args <- parser$parse_args()

metrics_to_run <- c("weighted_unifrac")

func_tab <- read.table(file = args$func_abun, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
abun_tab <- read.table(file = args$taxa_abun, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

in_tree <- ape::read.tree(args$in_tree)

null_out <- beta_div_contrib(metrics = metrics_to_run,
                            func_tab = func_tab,
                            abun_tab = abun_tab,
                            in_tree = in_tree,
                            ncores = args$num_cores,
                            return_objects = FALSE,
                            write_outfiles = TRUE,
                            outdir = args$outdir)
