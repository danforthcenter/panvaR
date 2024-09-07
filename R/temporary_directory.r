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