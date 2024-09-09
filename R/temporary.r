# Write an R function that checks if .panvar exits 
# if it does delete the files in it and return the path
# if the directory does not exits then create it and return the path

temporary_directory <- function(delete_files = FALSE) {
  dir_path <- ".panvar"
  
  if (dir.exists(dir_path) && delete_files) {
    # Delete the contents of the directory
    file.remove(list.files(dir_path, full.names = TRUE))
  }
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path)
  }
  
  return(dir_path)
}

# A function to roll my own tempfiles 
# Might as well for the simplicity of it
# This file just returns a path to a tempfile in 
# the .panvar directory

temp_file <- function(create_file = FALSE) {

  random_name <- paste(sample(LETTERS, 5, replace = TRUE), collapse = "") # A random letter generator

  if (create_file) {
    # Create the file
    file.create(paste0(".panvar/", random_name))
  }

  return(paste0(".panvar", "/", random_name))
  
}