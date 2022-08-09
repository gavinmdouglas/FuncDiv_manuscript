# FuncDiv_manuscript

Code used for showcasing the [`FuncDiv` R package](https://github.com/gavinmdouglas/FuncDiv), which is for an in-progress manuscript.

* `data/`
  * `Almeida_2019_dataset/` - input data used for assessing `FuncDiv` resource usage.
  * `random_forest/` - input and output files of random forest analysis.
  
* `scripts/` - code used to run manuscript analyses, split up by code for random forest and resource usage analyses, and minor code for results reported in text.


### Additional analysis details

Due to space restrictions, several details were left out of our description of the example application (i.e., the random forest analysis). These details are provided below.

The metagenomic dataset we analyszed was previously pre-processed as part of the curatedMetagenomicData ([Pasolli et al., 2017](https://www.nature.com/articles/nmeth.4468)). We downloaded the 2021/10/14 version of this processed dataset, which included the relative abundances of genera across samples, and a breakdown of which genera contribute each pathway per sample. Prior to building the models, we converted the genera abundances to relative abundances across all pathways and added a pseudocount of 0.1%. We then performed a centred log-ratio transformation of this table by sample. We also mean-centred and scaled by sample each contributional alpha diversity table. The random forest models were built with the ranger R package v0.14.1 ([Wright and Ziegler, 2017](https://www.jstatsoft.org/article/view/v077i01)).



