# # Dependencies
library(dplyr)
library(tidyr)
library(data.table)
# Load the argparse package
library(argparse)

parser <- ArgumentParser()

# Add arguments
parser$add_argument("--geno", help="Path to genotype file")
parser$add_argument("--pheno", help="Path to phenotype file")
parser$add_argument("--alpha", help="Significance level (default: 0.05)", default=0.05, type="double")
parser$add_argument("--output", help="Path to output directory (default: current working directory)", default=getwd(), type="character")

# Parse the command line arguments
args <- parser$parse_args()

# Check if files exist
if (!file.exists(args$geno)) {
  stop("Genotype file not found")
}

if (!file.exists(args$pheno)) {
  stop("Phenotype file not found")
}

# Check if output directory exists, and create it if it doesn't
if (!dir.exists(args$output)) {
  dir.create(args$output)
}

# Get base names of files
pheno_base <- basename(args$pheno)
geno_base <- basename(args$geno)

# Create output file name
output_name <- file.path(args$output, paste0(pheno_base, "_", geno_base))

# # Functions

#Rijan: Our genotype data comes with an ID field that has extraneous characters at the end
#that need to be cleaned up. This function takes the geno-type dataframe and cleans up the 
#ID field of extraneous characters.

replace_version_substring <- function(dt, col_name) {
    
    # Rijan: Not sure why but the current syntax in the function will modify the table in place, which is not ideal, hence a copy need to be made
    dt_copy <- copy(dt)

    dt_copy[, (col_name) := sub("\\.v[0-9]\\.\\d$", "", get(col_name))]

    return(dt_copy)

}

# Check if the column is a lineage column or a metadata column
# Lineage columns will only have 4 values 0,1,2 or NA - which 
# means we can make a heuristic to filter lineages out of the 
# genotype data and leave metadata behind
# Note: Is heuristic - might have false-positives.

is_lineage <- function(x) {
  all(x %in% c(0, 1, 2, NA))
}

# This function considers 1s and 2s to be mutations and documents there rates
# so that the rates can be used to filter for SNPs with high enough mutation 
# rates.

mutation_rates <- function(dt){
    
    # Calculate row sums for integer columns
    mutation_counts <- rowSums(dt == 1 | dt == 2, na.rm = TRUE)
    
    total_lines <- ncol(dt)
    
    mutation_percentages <- (mutation_counts / total_lines) * 100
    
    return(mutation_percentages)
}

# function to convert integerts to strings

integer_to_type <- function(x) {
  case_when(
    x == 0 ~ "wild",
    x == 1 ~ "hybrid",
    x == 2 ~ "mut",
    TRUE ~ as.character(x)
  )
}

adjusted_p_values <- function(genotype_data, phenotype_data) {


	# this will be used down the road because Tukey is not going to produce pairwise data all the time.
	filler_tukey_data <- data.table(
	    'mut-hybrid' = c(NA),
	    'wild-mut' = c(NA),
	    'wild-hybrid' = c(NA)
	)

	tukey_data <- data.table(
	    'mut-hybrid' = numeric(),
	    'wild-hybrid' = numeric(),
	    'wild-mut' = numeric()
	)

	for (col_num in 1:ncol(genotype_data)){
	
	    geno_data = genotype_data[[col_num]]
	
	    current_anova_table = cbind(phenotype_data, geno_data)
	
	    current_anova_table <- current_anova_table %>% as.data.table()
	
	    current_anova_table <- current_anova_table %>% drop_na()
	
	    # check if the table as at least 2 groups
	
	    # ANOVA cannot be done if there is only one group
	
	    if(length(unique(current_anova_table[[2]])) >= 2 ){
		
	        # get ANOVA
	        anova_output <- aov(current_anova_table[[1]] ~ current_anova_table[[2]])
	
	        # get Tukey
	        Tukey_output <- TukeyHSD(anova_output)
	
	        # make sure that we got tukey
	
	        if (!is.null(Tukey_output)){
			
	            # the `TukeyHSD` function formats its output somewhat awkwardly so we need to re-format it
	            result <- as.data.table(Tukey_output[1],keep.rownames = TRUE)
	
	            setnames(result, old = names(result), new = c("group", "diff", "lwr", "upr", "HSD")) # this is somewhat of a filler step and can be "hopped over" but makes the code more readeable
	
	            # make a valid tukey output
	
	            #tukey_output <- result[,c("group","HSD")] # these are the only fields that are worthwhile to us
	
	            names <- result %>% pull(1)
	
	            values <- result %>% pull(5)
	
	            current_hsd <- 
	                rbind(names,values) %>% 
	                as.data.table() %>% 
	                setnames(names) %>% 
	                slice(-1)
	
	
	            tukey_data <- rbindlist(list(tukey_data,current_hsd),fill = TRUE)
	
	            # if the Tukey test does not data add an empty row to the data.table
	
	        } else{
	            tukey_data <- rbindlist(list(tukey_data, filler_tukey_data),fill = TRUE)
	        }
	
	    # if the SNP does not have more than one type of group then add an empty row to the data.table
	    } else{
		
	        tukey_data <- rbindlist(list(tukey_data, filler_tukey_data),fill = TRUE)
	    }
	
	}

	return(tukey_data)
}

## Inputs

phenotype_data <- fread(args$pheno)

geno_base <- fread(args$geno)

# # Wrangling

# use the `replace_version_substring` function to remove extraneous chars from the ID field
impactful_SNPs <- replace_version_substring(geno_base,"ID")

genotype_data_effect_size_filter <- impactful_SNPs %>% 
    filter(EFFECT %in% c("MODERATE","HIGH"))

rates_of_mutation <- mutation_rates(genotype_data_effect_size_filter)

rates_of_mutation <- genotype_data_effect_size_filter %>% 
    select(is_lineage) %>% 
    mutation_rates()


genotype_data_effect_size_filter_mutation_rate_filter <- genotype_data_effect_size_filter[rates_of_mutation >= 1,]

genotype_data <- genotype_data_effect_size_filter_mutation_rate_filter %>% 
    select(is_lineage) %>%  # select columns that are lineages
    mutate(across(everything(), integer_to_type)) %>% # convert the types of genotype data from integer based to character based
    t() %>% # transpose the data to prepare for a left join
    as.data.table(keep.rownames = TRUE) %>% # turn things into a data.table for type reasons
    setnames("rn","Lineage") # change names for the sake of the leftjoin


names(phenotype_data)[1] = "Lineage" # The column is named to "PlantID" for the sake of standerdization

phenotype_genotype_data <- left_join(phenotype_data, genotype_data, by = "Lineage")

phenotype_anova_input <- phenotype_genotype_data[[2]]

genotype_anova_input <- phenotype_genotype_data %>%
    select(c(-1,-2))


adjusted_pvalues <- adjusted_p_values(genotype_anova_input, phenotype_anova_input)

tally_table <- t(genotype_anova_input)

filler_count_data <- data.table(
	    'mut' = numeric(),
	    'wild' = numeric(),
	    'hybrid' = numeric()
	)


for(i in 1:nrow(tally_table)){
    
    counts <- table(tally_table[i,]) %>% 
        t() %>% 
        as.list()
        
    
    filler_count_data <- rbindlist(list(filler_count_data, counts), fill = TRUE)
    
}

metadata_fields <- genotype_data_effect_size_filter_mutation_rate_filter %>% 
    select(!is_lineage)


summary_table <- cbind(metadata_fields,adjusted_pvalues, filler_count_data)

summary_table$bonferroni_value <- args$alpha / nrow(summary_table)

summary_table <- summary_table %>% 
    mutate(Significance = 'mut-hybrid' < bonferroni_value | 'wild-hybrid' < bonferroni_value | 'wild-mut' < bonferroni_value)


fwrite(summary_table,output_name,sep="\t")
