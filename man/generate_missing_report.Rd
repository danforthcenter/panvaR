% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_generate_missing_reports.r
\name{generate_missing_report}
\alias{generate_missing_report}
\title{generate_missing_reports}
\usage{
generate_missing_report(path_to_vcf_file)
}
\arguments{
\item{path_to_vcf_file}{Path to the vcf file.}
}
\value{
A list containing the following elements:
\item{table_for_snps}{A table that shows the missing values for the SNPS in the given VCF file.}
\item{table_for_lines}{A table that shows the missing values for the lines in the given VCF file.}
}
\description{
This function takes a vcf file and generate the missing reports for both lines and SNPs
}
\examples{
generate_missing_report("/path/to/your/vcf/file.vcf")
}
