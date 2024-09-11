# The general goal of the functions in this script is to verify  -
# that the inputs are correct.

# Check if the gwas table has the required data fields
check_gwas_table <- function(current_gwas_table){

    needed_call_names <- c("CHROM","BP","Pvalues")

    # are all the fields present in the current table
    test_fields <- all(needed_call_names %in% colnames(current_gwas_table))

    if(!test_fields){
        print("Your GWAS table either does not have the right data fields or has them mis-spelled.")

        print("You need to have these fields :- CHROM, BP and Pvalues.")

        return(FALSE)
    } else{
        return(TRUE)
    }
}