execute_snpsift <- function(file_path) {
    
    # Create a temporary file for the output
    snp_eff_table <- temp_file()

    # Set up the path to the SnpSift jar file
    jar_path <- system.file("java", "snpSift.jar", package = "panvaR")

    # Set up the arguments for the Java call
    binary_args <- c(
        "-jar", jar_path, "extractFields", file_path, 
        "CHROM", "POS", "ANN[*].FEATUREID", 
        "REF", "ALT", "ANN[*].EFFECT", 
        "ANN[*].AA", "ANN[*].IMPACT"
    )
    
    # Try to execute the command and catch any errors
    tryCatch(
        {
            error_message <- temp_file()  # Create a temporary file for error messages
            
            # Execute the system call to run SnpSift
            exec_wait(
                "java",  # Main binary call is Java
                args = binary_args,  # Pass the arguments
                std_out = snp_eff_table,  # Capture standard output in the table
                std_err = error_message  # Capture error messages
            )
            
            # snp_eff_table is the path -
            # let's return a data.table as well
            snp_eff_datatable <- fread(snp_eff_table)

            # Change the names of the table
            colnames(snp_eff_datatable) <- c(
                "CHROM",
                "BP",
                "GENE",
                "REF",
                "ALT",
                "EFFECT",
                "AA",
                "IMPACT"
            )

            # Make a list and return it

            return(
                list(
                    path = snp_eff_table,
                    table = snp_eff_datatable
                )
            )
            
        },
        error = function(e) {
            print(paste("Execution attempt produced error:", e$message))
            return(1)  # Return 1 on error
        }
    )
}