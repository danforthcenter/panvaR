# The general goal of the functions in this script is to verify  -
# that the inputs are correct.


# Check if the gwas table has the required data fields -
# if not then stop, if yes then return the right subset

check_gwas_table <- function(current_gwas_table){
  
  needed_call_names <- c("CHROM","BP","Pvalues")
  
  # check if the supplied table is a path or a table
  if(is.character(current_gwas_table)){
    if (!file.exists(current_gwas_table)) {
        stop("The GWAS table file path provided does not exist: ", current_gwas_table)
    }
    current_gwas_table <- data.table::fread(current_gwas_table)
  } else if (!is.data.frame(current_gwas_table)) {
    stop("Input to check_gwas_table must be a file path or a data.frame/data.table.")
  }
  
  # are all the fields present in the current table
  test_fields <- all(needed_call_names %in% colnames(current_gwas_table))
  
  if(!test_fields){
    missing_cols <- setdiff(needed_call_names, colnames(current_gwas_table))
    print(paste("You need to have these fields :-", paste(needed_call_names, collapse=", ")))
    stop(paste("The GWAS table does not have the right data fields. Missing:", paste(missing_cols, collapse=", ")))
  } else{
    # If true select for the three fields that we need
    current_gwas_table <- as.data.table(current_gwas_table) %>%
      select(all_of(needed_call_names))
    return(current_gwas_table)
  }
}