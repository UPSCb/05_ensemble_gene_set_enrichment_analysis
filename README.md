# RnaSeqTutorial05

Ensemble Gene Set Enrichment Analysis tutorial

## Installation

In R do the following to install the dependencies

```R
if(!require("BiocManager")) {install.packages("BiocManager")}
BiocManager::install(c("dplyr","edgeR","EGSEA","EGSEAdata","here","learnr","limma","org.Hs.eg.db","readr","tibble","UpSetR"))
```

Then install the package using _e.g._ `pak`

 ```R
 if(!require("pak")){BiocManager::install("pak")}
 pak::pkg_install("UPSCb/Tutorial_05_ensemble_gene_set_enrichment_analysis")
 ```

### Troubleshooting

If you get an error message stating: `could not find tools necessary to compile a package`, then:

* on **Windows**, you will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), select the right version for you R installation.

* on **Mac OS**, make sure that xcode is installed. In the terminal, run `xcode-select --install`.

If the above do not resolve the installation issue, then run this instead:

```R
if(!require("pak")){BiocManager::install("devtools")}
devtools::install_github("UPSCb/Tutorial_05_ensemble_gene_set_enrichment_analysis")
```

## Getting started

This is the first in a series of tutorials.

1. 01_data_science_intro: an introduction to data science using the `palmerpenguins` package
2. 02_exploratory_data_analysis: the exploratory data analysis of a public RNA-Seq dataset
3. 03_differential_gene_expression: the differential expression DE analysis of the same dataset
4. 04_gene_set_enrichment_analysis: gene set enrichment analysis (GSEA) on the DE gene sets
5. 05_ensemble_gene_set_enrichment_analysis: ensemble GSEA of the DE gene sets and comparison to GSEA
6. 06_differential_transcript_usage: the differential transcript usage analysis of the same dataset
7. 07_non_model_organism_exploratory_data_analysis: the exploratory data analysis of a Norway spruce RNA-Seq dataset
8. 08_non_model_organism_differential_gene_expression: the differential expression DE analysis of the Norway spruce RNA-Seq dataset

To start the tutorial run:

```{r tutorial}
learnr::run_tutorial("05_ensemble_gene_set_enrichment_analysis", package = "RnaSeqTutorial05")
```

## Accessing the code

### Tutorials

The tutorial is available as `.Rmd` files in the `tutorials` directory:

```{r tutorial list}
dir(system.file(package="RnaSeqTutorial05","tutorials"))
```

```{r tutorial paths}
list.files(path=system.file(package="RnaSeqTutorial05","tutorials"),
           pattern="*.Rmd",
           recursive=TRUE)
```

### Scripts

There may be script(s) available in the `scripts` directory:

```{r script paths}
list.files(path=system.file(package="RnaSeqTutorial05","scripts"),
           pattern="*.R",
           recursive=TRUE)
```
