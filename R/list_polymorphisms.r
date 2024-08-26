# The goal of this script is to translate the bash code into R
# This really only needs the process_gene() function translated

does_gene_exist <- function(gene_name, gene_location){

    grep_call <- "grep"
	
    grep_args = c(
    "--F", gene_name,
    gene_location 
    )

	tryCatch(
		{
			try <- exec_wait(
				grep_call,
				args = grep_args,
				std_out = TRUE,
				std_err = TRUE
			)
		},
		error = function(e){
            # Custom error message
			print("There was an error with running grep on your input.")
            print(paste("Execution attempt produced error:-", e$message))
			stop()
            1 # Return 1 on error
        }
	)
}

tabix_subset <- function(start, stop, subset_file, vcf_file){

    tabix_call <- "tabix"
	
    the_args = c(
    "-h", paste0(vcf_file,":",start, "-", stop),
    ">",
    subset_file
    )

	tryCatch(
		{
			try <- exec_wait(
				grep_call,
				args = the_args,
				std_out = TRUE,
				std_err = TRUE
			)
		},
		error = function(e){
            # Custom error message
			print("There was an error with running grep on your input.")
            print(paste("Execution attempt produced error:-", e$message))
			stop()
            1 # Return 1 on error
        }
	)

}

process_gene_function <- function(gene_name, vcf_file, gene_location){


    # open the gene location file as a data.table
    gene_location <- fread(gene_location)

    # use the does_gene_exist() function to check if the gene exists
    if(does_gene_exist(gene_name, gene_location) == 0){
        stop(paste0("The gene does not exist in the gene location file, might be time to check your inputs for ", gene_name))
    }

    # For this specific gene get the start and stop position from the data table 
    # that has the gene_location data
    current_gene <- gene_location %>% 
        filter(Gene == "Sevir.J000601")

    gene_start <- current_gene$Ext_Start

    gene_stop <- current_gene$Ext_Start

    # using tabix subset the vcf file

    # first create a tempfile that will hold the tabix subset
    current_tabix_subset <- tempfile()


}

process_gene <- function(gene_name, vcf_file, gene_location){

    # check if the gene_location files exists
    if(!file.exists(gene_location)){
        stop("The gene location file does not exist")
    }

    # check if the vcf_file exists
    if(!file.exists(vcf_file)){
        stop("The vcf file does not exist")
    }

    # Using the proper_tbi function defined in R/general_functions to check if the vcf file is indexed
    proper_tbi(vcf_file)

    # Now that all the tests have passed pass the params to the checker the function
    # that can safely assume that the input is most likely properly formatted
    process_gene_function(gene_name, vcf_file, gene_location)

}