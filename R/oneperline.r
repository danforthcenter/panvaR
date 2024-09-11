split_vcf_eff <- function(input_file = NULL, output_file = NULL) {
  # Function to process a single line
  process_line <- function(line) {
    if (grepl("^#", line)) {
      return(line)
    }
    
    fields <- strsplit(trimws(line), "\t")[[1]]
    info <- fields[INFO_FIELD_NUM]
    info_fields <- strsplit(info, ";")[[1]]
    
    eff_field <- grep("^(EFF|ANN)=", info_fields, value = TRUE)
    if (length(eff_field) == 0) {
      return(line)
    }
    
    field_name <- sub("=.*", "", eff_field)
    effs <- strsplit(sub("^[^=]+=", "", eff_field), ",")[[1]]
    
    other_info <- info_fields[!grepl("^(EFF|ANN)=", info_fields)]
    
    pre <- paste(fields[1:(INFO_FIELD_NUM-1)], collapse = "\t")
    post <- paste(fields[(INFO_FIELD_NUM+1):length(fields)], collapse = "\t")
    
    output_lines <- character(length(effs))
    for (i in seq_along(effs)) {
      new_info <- if (length(other_info) > 0) {
        paste(c(other_info, paste0(field_name, "=", effs[i])), collapse = ";")
      } else {
        paste0(field_name, "=", effs[i])
      }
      output_lines[i] <- paste(pre, new_info, post, sep = "\t")
    }
    return(output_lines)
  }

  INFO_FIELD_NUM <- 8  # R uses 1-based indexing

  # Main processing
  if (is.null(input_file)) {
    con_in <- stdin()
  } else {
    con_in <- file(input_file, "r")
  }

  # Handle output
  if (is.null(output_file)) {
    output_file <- temp_file()
  }
  con_out <- file(output_file, "w")

  # Main processing
  tryCatch({
    while (TRUE) {
      line <- readLines(con_in, n = 1)
      if (length(line) == 0) break
      processed_lines <- process_line(line)
      writeLines(processed_lines, con_out)
    }
  }, finally = {
    if (!is.null(input_file)) close(con_in)
    close(con_out)
  })

  # Return the output file name
  return(output_file)
}

