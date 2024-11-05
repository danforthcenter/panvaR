# The general goal of the functions in this script is to verify  -
# that the inputs are correct.


# Check if the gwas table has the required data fields -
# if not then stop, if yes then return the right subset

check_gwas_table <- function(current_gwas_table){

    needed_call_names <- c("CHROM","BP","Pvalues")

    # check if the supplied table is a path or a table
    if(is.character(current_gwas_table)){
        current_gwas_table <- fread(current_gwas_table)
    }

    # are all the fields present in the current table
    test_fields <- all(needed_call_names %in% colnames(current_gwas_table))

    if(!test_fields){
        print("You need to have these fields :- CHROM, BP and Pvalues.")
        stop("The GWAS table does not have the right data fields or has them mis-spelled.")
		stop()
    } else{
        # If true select for the three fields that we need
        current_gwas_table <- current_gwas_table %>%
            select(CHROM, BP, Pvalues)
        return(current_gwas_table)
    }
}