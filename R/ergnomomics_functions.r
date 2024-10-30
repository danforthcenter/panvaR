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

# This function get the extention of the file supplied
extention_func <- function(path_to_file,super_name = FALSE, include_dir = FALSE){

    # use this code to get the basename of the file without the extension
    name_string_items <- strsplit(basename(path_to_file),"\\.")[[1]]

    # Just get the last item in the list
    extension <- name_string_items[length(name_string_items)]

    return(extension)
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
	
	core_count = detectCores() - 1 
	if(core_count < 1){
	    core_count = 1
	}

	return(core_count)
}

# A function to get all the BPs to keep
snps_to_keep <- function(ld_table){

    # This function only takes the very striclty defined
    # and shaped table as made by the function `ld_filtered_snp_list`
    # in the file `ld_filtered_snp_list.R`

    # This function is needed because the way
    # Plink2 calculates LD, it leaves out the tag snp
    # and keeps everything else
    # that is all that we are accomodating for.
    
    # Get other snps
    other_snps <- ld_table %>%
        select(Subject_snp_chrom,Subject_snp_bp)  %>%
        rename(CHROM = Subject_snp_chrom,BP = Subject_snp_bp)

    # Get tag snps
    tag_snp <- ld_table %>%
        select(Tag_snp_chrom, Tag_snp_bp) %>%
        distinct() %>%
        rename(CHROM = Tag_snp_chrom,BP = Tag_snp_bp)
    
    keep_list <- rbind(other_snps,tag_snp)

    return(keep_list)
}

# Sanitize the tables to keep
keep_table_sanitizer <- function(table) {

    keep_table_path <- temp_file()
    
    if(is.character(table)){

        keep_table <- fread(table)

		# TODO:- How to make sure the Chrs -
		# Come before the loci in the the -
		# keep tables?

        keep_table %>%
            fwrite(keep_table_path,sep = "\t",col.names = FALSE)
    } else{
        keep_table <- as.data.table(table,col.names=TRUE)

		keep_table <- keep_table %>%
            select(CHROM,BP)

        keep_table %>%
            fwrite(keep_table_path,sep = "\t",col.names = FALSE)
    }

    return(keep_table_path)
    
}

# A function to get the tag snp out of the gwas results table
tag_snp_func <- function(gwas_results){

    # This function assumes that the supplied tables has the following fields
    # 1. 'the_bp`
    # 2. 'the_chromosomes'
    # 3. 'the_pvalues' -> -log10() of the pvalues 

    # Just double checking to make sure that the table is arranged properly
    
    current_table <- gwas_results %>% 
        arrange(desc(as.numeric(Pvalues)))

    tag_snp_row <- current_table %>%
        slice(1)

    tag_snp <- list(
        tag_snp_bp = as.integer(tag_snp_row$BP), 
        tag_snp_chromosome = tag_snp_row$CHROM

    )

    return(tag_snp)
}

# The user might supply the tag snp formatted as - 
# "<chr>:<bp>". This function provides the convienience -
# to make the split easily.

tag_snp_splitter <- function(tag_snp){
    
    tag_snp_items <- strsplit(tag_snp,":")[[1]]
    
    tag_snp_chrom <- tag_snp_items[1]
    
    tag_snp_bp <- as.double(tag_snp_items[2])
    
    return_list <- list(
        chrom = tag_snp_chrom,
        bp = tag_snp_bp
    )
    
    return(return_list)
}

# Make the table of LD that retains the tag SNP
ld_table_maker <- function(ld_table){

    ld_table_subject <- ld_table %>% 
        select(
            Subject_snp_chrom,
            Subject_snp_bp,
            Phased_r2
        ) %>%
        rename(CHROM = Subject_snp_chrom,BP = Subject_snp_bp, LD = Phased_r2)

    ld_table_tag <- ld_table %>%
        select(Tag_snp_chrom, Tag_snp_bp) %>%
        distinct() %>%
        rename(CHROM = Tag_snp_chrom, BP = Tag_snp_bp)
    
    ld_table <- rbind(ld_table_subject,ld_table_tag, fill = TRUE)

    return(ld_table)
    
}


# Function:- `window_unit_func` 
# Context:- The user can either supply a numeric and granular window in numbers or bp -
# or supply a kilobase window as characters "500kb".
# Goal:- Take an input from the user and return the window value.

window_unit_func <- function(window_value) {

    # Check if the input is of type character
    if (typeof(window_value) == 'character') {

        # Convert the user-supplied input to numeric and multiply it by 1000
        kb_numeric = as.numeric(gsub("KB", "", window_value, ignore.case = TRUE))

        # Check if conversion was successful
        if (is.na(kb_numeric)) {
            stop("Error: Input could not be converted to numeric.")
        }

        kb_window = kb_numeric * 1000

        return(kb_window)
    }

    # If the input is already numeric, return it
    if (typeof(window_value) == 'double') {
        return(window_value)
    }

    stop("Error: Unsupported input type.")
}