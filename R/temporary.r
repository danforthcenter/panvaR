# R's temporary directory and temp file functions are limited.
# This is why we will have to write out own.

# User have expressed the desire to have the output directory not be -
# a hidden directory.

# Function to manage the temporary directory
# Checks if the directory exists; if not, creates it. Deletes contents if specified.

temporary_directory <- function(delete_files = FALSE, working_directory = "panvar") {
  
  # Use provided working directory or default to "panvar"
  dir_path <- ifelse(is.null(working_directory), "panvar", working_directory)

  # Delete the contents of the directory if it exists and delete_files is TRUE
  if (dir.exists(dir_path) && delete_files) {
    file.remove(list.files(dir_path, full.names = TRUE))
  }

  # Create the directory if it does not exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path)
  }

  return(dir_path)
}

# Function to generate a temporary file path in the specified directory
# Optionally creates the file

temp_file <- function(create_file = FALSE, working_directory = "panvar") {
  
  # Use provided working directory or default to "panvar"
  dir_path <- ifelse(is.null(working_directory), "panvar", working_directory)
  
  # Generate a random file name
  random_name <- paste0(sample(LETTERS, 5, replace = TRUE), collapse = "")

  # Create the final path
  final_path <- file.path(dir_path, random_name)

  # Create the file if specified
  if (create_file) {
    file.create(final_path)
  }

  return(final_path)
}
