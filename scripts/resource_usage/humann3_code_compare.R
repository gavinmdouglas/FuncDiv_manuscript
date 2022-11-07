# A quick comparison of the efficiency of running HUMAnN3's example R code from their tutorial vs. FuncDiv on the same test dataset.
# Note that this is based on the instructions here (accessed Nov. 3, 2022):
# https://github.com/biobakery/biobakery/wiki/Metagenomic-Visualizations#compute-alpha-diversity-for-each-pathway-in-each-sample

rm(list = ls(all.names = TRUE))

library(FuncDiv)
library(tidyverse)
library(vegan)

setwd("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/Almeida_2019_dataset/prepped_infiles")

MGX_pwys_raw = read_tsv("/data1/gdouglas/projects/contrib_div/FuncDiv_manuscript/data/HUMAnN3_tutorial_table/mgx_pwys.tsv.gz")

MGX_pwys = MGX_pwys_raw |> 
  separate(Pathway, 
           into = c("pwy", "taxa"),
           sep = "\\|") |> 
  filter(!(pwy %in% c("UNMAPPED", "UNINTEGRATED"))) |>
  filter(!(is.na(taxa) | taxa == "unclassified"))


alpha_div = function(pwy_data) {
  
  relative_pwy = pwy_data |>
    mutate(across(c(-pwy, -taxa), ~ .x / sum(.x))) |>
    mutate(across(c(-pwy, -taxa), ~ replace_na(.x, 0))) |>
    select(-taxa) |>
    group_by(pwy)
  
  if (nrow(relative_pwy) == 1) {
    res = relative_pwy |>
      summarise(across(.fns = ~ if_else(sum(.x) > 0,
                                        0,
                                        -1)))
  } else {
    res = relative_pwy |>
      summarise(across(.fns = ~ if_else(sum(.x) > 0,
                                        diversity(.x, index = 'simpson'),
                                        -1)))
  }
  
  sample_values = res |> select(-pwy) |> unlist()
  
  all_null = all(sample_values == -1)
  
  if (all_null) {
    return(NULL)
  } else {
    return(res)
  }
  
}

# Ran each command 10 times and then averaged to get proper estimate of runtime.


humann3_code_runtime <- as.numeric()
funcdiv_code_runtime <- as.numeric()
for (i in 1:10) {
  start_time <- Sys.time()
  MGX_a_div = MGX_pwys |> 
    group_split(pwy) |> 
    map(alpha_div) |> 
    bind_rows()
  end_time <- Sys.time()
  humann3_code_runtime <- c(humann3_code_runtime, end_time - start_time)

  start_time <- Sys.time()
  MGX_pwys_contib <- reshape2::melt(data = MGX_pwys, ids = c("pwy", "taxa"), value.name = "abun", variable.name = "sample")
  MGX_pwys_contib <- MGX_pwys_contib[-which(MGX_pwys_contib$abun == 0), ]
  MGX_a_div_FuncDiv <- alpha_div_contrib(metrics = "gini_simpson_index",
                                         ncores = 1,
                                         contrib_tab = MGX_pwys_contib,
                                         samp_colname = "sample",
                                         func_colname = "pwy",
                                         taxon_colname = "taxa",
                                         abun_colname = "abun")
  end_time <- Sys.time()
  funcdiv_code_runtime <- c(funcdiv_code_runtime, end_time - start_time)
}

mean(humann3_code_runtime)
sd(humann3_code_runtime)

mean(funcdiv_code_runtime)
sd(funcdiv_code_runtime)
