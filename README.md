# PanvaR

An R package for prioritizing candidate genes from GWA studies. 


## Installation 

Ensure remotes package is installed. 
``` r
install.packages("remotes")
```

Install panvaR:
``` r
remotes::install_github("danforthcenter/panvaR", build_vignettes = TRUE)
```

To get started, view the vignette:
```r
browseVignettes("panvaR")
```

We have also implemented a GUI for visualizing results. To use the GUI:
```r
panvar_gui()
```
# Reproduction

Please visit the following Zenodo repository for result reproduction: {Placeholder for Collin to add Zenodo link}

# Dependencies

While PanvaR uses R to visualize and combine data, it relies on some utilities that are not available on all operating systems. Consequently, we are not supporting Windows operating systems at this time. 

## Non-R Dependencies

`panvaR` is an R package that also relies on several external, command-line bioinformatics tools to function correctly. These tools must be installed and available in your system's `PATH`.

The required external dependencies are:

*   **PLINK 2 (`plink2`)**: (Tested with Version 2.0.0-a.6.9LM) The primary workhorse for genetic data manipulation. It is used for:
    *   Converting VCF files to PLINK's binary format (`.bed`, `.bim`, `.fam`).
    *   Quality control filtering (e.g., Minor Allele Frequency, missingness).
    *   Calculating Linkage Disequilibrium (LD) to find SNPs correlated with a tag SNP.
    *   Extracting specific genomic regions (windowing).
    *   Generating reports on missing data.

*   **BCFtools (`bcftools`)**: (Tested with version 1.21) Primarily used for indexing and retrieving data from compressed VCF/BCF files (`.vcf.gz` or `.bcf`), creating `.tbi` index files using the `tabix` command.  It may also be used to list chromosomes, though this is a secondary function.

*   **VCFtools (`vcftools`)**: (Tested with Version 0.1.16) Used for filtering and manipulating VCF files. `panvaR` uses it to create new VCF files based on SNP lists (e.g., those identified through LD analysis) and potentially for other filtering operations. It's worth noting that some functionality might overlap with BCFtools, and future development could explore consolidating these dependencies.

*   **tabix**: (Tested with version 1.21)A generic tool for indexing and retrieving data from tab-delimited files. `panvaR` uses it to create and check for index files (`.tbi`) for compressed VCF files (`.vcf.gz`), which is essential for fast data access.

*   **Java Runtime Environment (JRE)**: (Version 17 or later) Required to run SnpSift, which is used for annotation extraction.

*   **SnpSift**: (Version 4.3t or later) This tool is part of the SnpEff suite and is used to parse and extract functional annotation information (e.g., variant effect, impact, gene name) from the `INFO` field of a VCF file. A `snpSift.jar` file is included within the `panvaR` package, but a system-wide Java installation is still necessary to execute it.

## R Dependencies

TLDR: So long as you are on the most up-to-date versions of the following packages, we do not expect many breaking-level issues from the R dependencies. That being said, the workflow has been tested with the following packages and listed are their versions.

The requested R package information is as follows:

| Package Name | Version Details |
|---|---|
| bigsnpr | 1.12.18 (2024-11-26) CRAN (R 4.4.3) |
| bigstatsr | 1.6.1 (2024-09-09) CRAN (R 4.4.3) |
| data.table | 1.17.0 (2025-02-22) CRAN (R 4.4.3) |
| DT | 0.33 (2024-04-04) CRAN (R 4.4.3) |
| dplyr | 1.1.4 (2023-11-17) CRAN (R 4.4.3) |
| forcats | 1.0.0 (2023-01-29) CRAN (R 4.4.3) |
| ggplot2 | 3.5.1 (2024-04-23) CRAN (R 4.4.3) |
| knitr | 1.49 (2024-11-08) CRAN (R 4.4.3) |
| lubridate | 1.9.4 (2024-12-08) CRAN (R 4.4.3) |
| modelr | 0.1.11 (2023-03-22) CRAN (R 4.4.3) |
| plotly | 4.10.4 (2024-01-13) CRAN (R 4.4.3) |
| purrr | 1.0.4 (2025-02-05) CRAN (R 4.4.3) |
| readr | 2.1.5 (2024-01-10) CRAN (R 4.4.3) |
| rmarkdown | Not found in session info |
| shiny | 1.10.0 (2024-12-14) CRAN (R 4.4.3) |
| shinyBS | 0.61.1 (2022-04-17) CRAN (R 4.4.3) |
| shinyFiles | 0.9.3 (2022-08-19) CRAN (R 4.4.3) |
| shinyjs | 2.1.0 (2021-12-23) CRAN (R 4.4.3) |
| stringr | 1.5.1 (2023-11-14) CRAN (R 4.4.3) |
| sys | 3.4.3 (2024-10-04) CRAN (R 4.4.3) |
| tibble | 3.2.1 (2023-03-20) CRAN (R 4.4.3) |
| tidyr | 1.3.1 (2024-01-24) CRAN (R 4.4.3) |
| tidyverse | 2.0.0 (2023-02-22) CRAN (R 4.4.3) |
| usethis | 3.1.0 (2024-11-26) CRAN (R 4.4.3) |
| parallel | Base R package (not listed in session info as separate entry) |

# Operating systems support

For now, `panvaR` has been tested and runs successfully on POSIX compliant operating systems (MacOS, various Linux distributions). We have **not** been able to test the package on WSL2 or WSL1 - while the WSL2 and WSL1 are POSIX compliant, because we have not been able to test the package properly in the said environments we cannot yet comment on whether the package will work there.

As of this commit `PanvaR` definitely does not work on Windows due the workflows reliance on `HTSLib` based tools `BCFtools`, `VCFtools`, and `tabix` which do not have binaries of Windows. It is possible that Windows could be supported in the future as we explore the removal of `HTSLib` based dependencies.
