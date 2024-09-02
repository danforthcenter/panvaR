# A function to return just the basename of a file
# This is something I have been using often - as such it is useful to just have around
base_name_func <- function(path_to_file,super_name = FALSE, include_dir = FALSE){

    # use this code to get the basename of the file without the extension
    name_string_items <- strsplit(basename(path_to_file),"\\.")[[1]]

    #this gets everything other than the last item
    # this is also a good fail-safe against supplying a file that has no extension
    name_items <- name_string_items[-length(name_string_items)]

    # do you need the directory name of the file
    dir_name <- dirname(path_to_file)

    if (super_name){
        the_base = paste(name_items, collapse = ".")
    } else {
        the_base = name_items[1]
    }

    if (include_dir){
        output_name = paste0(dir_name,"/",the_base,sep = "")
    } else {
        output_name = the_base
    }

    return(output_name)
}


# Is the binary you want to use available in the shell?
locate_bin_on_shell <- function(command){

    the_command <- paste0(command)

    result <- tryCatch(
        {
            test <- sys::exec_wait(command,
            args = "-v",
            std_out = TRUE,
            std_err = FALSE
            )
            0 # Return 0 on success
        },
        error = function(e) {
            # Custom error message
            print(paste0(command, " is either not installed or not available in this shell."))
            print(paste("Execution attempt produced error:-", e$message))
            1 # Return 1 on error
        }
    )

    if (result == 0){
        return(TRUE)
    }
	else{
		return(FALSE)
	}
}

# Count and provide the number of cells that should be used
good_core_count <- function(){
	
	core_count = nb_cores() - 1 
	if(core_count < 1){
	    core_count = 1
	}

	return(core_count)
}

# A function to get all the BPs to keep
snps_to_keep <- function(ld_table){

    # This functions only takes the very striclty defined
    # and shaped table as made by the function `ld_filtered_snp_list`
    # in the file `ld_filtered_snp_list.R`

    # This function is need because the way
    # Plink2 calculates LD, it leaves out the tag snp
    # and keeps everything else
    # that is all that we are accomodating for.
    
    # Get other snps
    other_snps <- ld_table %>%
        select(Subject_snp_chrom,Subject_snp_bp)  %>%
        rename(Chrom = Subject_snp_chrom,bp = Subject_snp_bp)

    # Get tag snps
    tag_snp <- ld_table %>%
        select(Tag_snp_chrom, Tag_snp_bp) %>%
        distinct() %>%
        rename(Chrom = Tag_snp_chrom,bp = Tag_snp_bp)
    
    keep_list <- rbind(other_snps,tag_snp)

    return(keep_list)
}

# Sanitize the tables to keep
keep_table_sanitizer <- function(table) {

    keep_table_path <- tempfile()
    
    if(is.character(table)){

        keep_table <- fread(table)

        keep_table %>%
            fwrite(keep_table_path,sep = "\t",col.names = FALSE)
    } else{
        keep_table <- as.data.table(table)

        keep_table %>%
            fwrite(keep_table_path,sep = "\t",col.names = FALSE)
    }

    return(keep_table_path)
    
}
