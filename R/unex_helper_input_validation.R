# The tag snp should be supplied in the as
# "chrom:bp" format
tag_snp_check <- function(tag_snp) {
    if (grepl(":", tag_snp)) {
        return(TRUE)
    } else {
        message("The tag snp should be supplied in the as
             \"chrom:bp\" format")
        return(FALSE)
    }
}
