genes_in_range <- function(file_name, chrom, start, stop){

    # open the file as a data.table
    dt <- data.table::fread(file_name)

    # check if the data.table has the column chrom and Ext_start
    if(!"chrom" %in% colnames(dt) || !"Ext_start" %in% colnames(dt)){
        stop("The gene location file does not have the column chrom and Ext_start - these are need with impeccable spelling!")
    }

    # subset the data.table for for the given start and stop values on the column `Ext_start`
    subset = dt[chrom == chrom & Ext_start >= start & Ext_start <= stop]

    # Deduplcate subset on the column gene
    return(unique(subset$gene))

}


list_of_genes <- function(file_name, chrom, start, stop){

    # check if the file exists
    if(!file.exists(file_name)){
        stop("The file path that you supplied does not exist")
    }

    # check if start and stop are numeric and non-negative
    if(!is.numeric(start) || !is.numeric(stop) || start < 0 || stop < 0){
        stop("start and stop must be numeric and non-negative.")
    }

    # check to make sure that stop is greater than start
    if(stop < start){
        stop("stop must be greater than start.")
    }

    # Now that all the params have been checked pass them to the function that can safely work under the assumption that everything is as it should be
    return(genes_in_range(file_name, chrom, start, stop))
}