# Panvar

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

To use the GUI:
```r
panvar_gui()
```
## System Dependencies

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