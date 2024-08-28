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