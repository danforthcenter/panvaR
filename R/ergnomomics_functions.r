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
