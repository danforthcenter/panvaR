# %% [markdown]
# # Dependencies

# %% [code] {"_execution_state":"idle"}
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(gtable)
library(patchwork)
library(stringr)
library(feather)
library(data.table)

# %% [markdown]
# # Global constants

# %% [code]
transcript_trailer_to_remove <- "\\.v3.1$"

##Transcript Identifer - For grabbing individual genes from transcript IDs ## '$' denotes trailing number
##Setaria = "\\.[1-9]$"

individual_transcript_identifier <- "\\.[1-9]$"

##Plot specific names
population_name <- "Sorghum BAP"

# %% [code]
gene_db_path = "/kaggle/input/panvar-phenotype/Sbicolor_313_v3.1.gene.gff.feather"

# %% [code]
annotation_db_path = "/kaggle/input/panvar-phenotype/Sbicolor_454_v3.1.1.annotation.feather"

# %% [code]
PAV_subpops_path = "/kaggle/input/panvar-phenotype/BAP_pop_assignments_q8.csv"

# %% [code]
PAV_calls_path = "/kaggle/input/panvar-phenotype/BAP_PAV_calls.feather"

# %% [markdown]
# # Input data

# %% [code]
gene_db <- as.data.table(read_feather(gene_db_path))

# %% [code]
annotation_db <- as.data.table(read_feather(annotation_db_path))

# %% [code]
PAV_subpops <- fread(PAV_subpops_path)

# %% [code]
PAV_calls <- as.data.table(read_feather(PAV_calls_path))

# %% [markdown]
# # Interacting with phenotype data

# %% [code]
phenotype = read.table("/kaggle/input/panvar-phenotype/BAP_WSC_pheno.txt",sep ="\t",header = TRUE)

# %% [markdown]
# `snp_targets` holds the information to a list of impactful genes using something similar to `glob` in python
# 
# - [ ] `snp_targets` needs to be abstracted in such a way that it is a list of genes and not some files that become a list

# %% [code]
snp_targets = list.files("/kaggle/input/panvar-phenotype/individual_files_impactful",pattern="impactful.txt",recursive = T, full.names = F)

# %% [markdown]
# - go through the `snp_targets` 
# 
#     - Open each `impactful_file`
#     
#         - Remove unwanted suffix from the `ID` field
#         
# - [ ] why bother with `size > 0` in `mutation_db` when you can have the same impact on the same data field upstream with `alt_portion`
# 
# - [ ] The way levels are being used in this codebase really irks me

# %% [markdown]
# change log
# 
# | Old      | new          |
# |----------|--------------|
# | `$PER_ALT` | `$alt_portion` |

# %% [markdown]
# ## major tables
# 
# 1. pre_mutation_db -> An iterator that will hold one gene's polymorphism data
# 
# 2. mutation_db -> A subset of `pre_mutation_db` with out the lineages
# 
# 3. pheno_geno_table_sig -> Iterator used to hold pre_mutation_db filtered for "HIGH","MODERATE" values in the `EFFECT` and for `gene` in the `ID` field. `gene` is an iterator that goes through `gene_of_interest`

# %% [code]
for(single_snp_target in snp_targets){
    
    path_to_target <- paste("/kaggle/input/panvar-phenotype/individual_files_impactful",single_snp_target, sep="/")
    single_snp_target <- gsub(single_snp_target, pattern = ".txt", replacement = "")
    
    
    pre_mutation_db <- fread(path_to_target, sep = "\t", header = T)
    
    pre_mutation_db$ID <- gsub(x = pre_mutation_db$ID, pattern = "\\.v3.1$", replacement = "")
    
    number_genos <- length(pre_mutation_db) - 9 # Rijan this assumes that the first 9 fields are the ones with the the meta data.
    
    pre_mutation_db$alt_portion = (rowSums(pre_mutation_db[,9:length(pre_mutation_db)] == 1 | pre_mutation_db[,9:length(pre_mutation_db)] == 2, na.rm = TRUE) / number_genos) * 100 # Rijan: this is now fixed to get the right percentage of the alts in each polymorphism

    ##Rijan: Filter the `alt_portion` for the value above 1 percent
    
    pre_mutation_db <- pre_mutation_db %>% filter(alt_portion >= 1)
    
    mutation_db <- cbind(pre_mutation_db[,1:8],alt_portion = pre_mutation_db$alt_portion)
    
    mutation_db$size <- mutation_db$alt_portion / 10
    
    # Rijan: Map the variant to size. Off the top of my head, I do not remember this being used downstream but let's wait and see
    shape_map <- c(
        "missense_variant" = 21,
        "synonymous_variant" = 22,
        "disruptive_inframe_deletion" = 23
    )
    
    mutation_db$shape <- shape_map[mutation_db$TYPE]

    mutation_db$shape[is.na(mutation_db$shape)] <- 25
    
    # Rijan: Not sure what the goal here is but it is being done. This can most likely be removed.
    mutation_db_filt <- droplevels(mutation_db %>% filter(size > 0))
    mutation_db_filt$ID <- as.factor(mutation_db_filt$ID)
    
    # Rijan: This might be used downstream but for now it is a place holder
    
    output_db <- pre_mutation_db
    
    # Rijan: This use of the levels syntax is really, I mean really bugging me. 
    genes_of_interest <- unique(output_db$ID)
    
    # A placeholder for the variant db that weill be used down the road
    
    variant_sig_db <- data.table()
    
    # Rijan: `genes_of_interst` holds the unique values in output_db's ID field
    for(gene in genes_of_interest){
        
        # Rijan: No other way around this maybe?
        effects_of_interest <- c("HIGH","MODERATE")
        
        # filter pre_mutation_db on gene and effects_of_interst field
        pheno_geno_table_sig <- pre_mutation_db %>% filter(EFFECT %in% effects_of_interest, ID == gene)
        
        # the pheno_geno_table_sig will most likely be used to calculate something and it probably needs to be in a specific shape for that to work, which might be why `alt_portion` is being set to NULL
        pheno_geno_table_sig$alt_portion <- NULL
        
        if(nrow(pheno_geno_table_sig) >= 1){
            
            genos_table_sig <- pheno_geno_table_sig %>% 
                select(9:ncol(pheno_geno_table_sig))
            
            # Rijan transpose the table
            
            genos_table_sig <- t(genos_table_sig)
            
            genos_table_sig <- as.data.table(genos_table_sig,keep.rownames = TRUE)
            
            genos_table_sig2 <- genos_table_sig %>%
                mutate(across(everything(), ~ case_when(
                . == " " ~ "",
                . == "0" ~ "Ref/Ref",
                . == "1" ~ "Ref/Alt",
                . == "2" ~ "Alt/Alt",
                TRUE ~ as.character(.)
              )))
            
            # Rijan: the following three lines are just to assign names to the the `genos_table_sig2` in a where the column names reflect the genotype of the lines haplotypes
            
            pheno_geno_table_sig$names <- paste(pheno_geno_table_sig$CHROM, pheno_geno_table_sig$POS, pheno_geno_table_sig$TYPE, sep = "_")
            
            colnames(genos_table_sig2) <- c("PLANTID",t(pheno_geno_table_sig$names))
            
            traits_to_graph <- colnames(genos_table_sig2)[2:length(genos_table_sig2)]
      
            ##Merge variants with phenotypes, drop PlantIDs
            geno_pheno_full_table_sig <- merge(phenotype, genos_table_sig2, by.x = 1, by.y = 1, all=F)
            
            geno_pheno_full_table_sig[,1] <- NULL

            #geno_pheno_full_table_sigWhat does the P adjusted filed of this TukeyHSD(aov(anova_table[,1] ~ anova_table[,2]), ordered = TRUE) say?



            # assume 'dt' is your data.table
            for (col in names(geno_pheno_full_table_sig)[2:ncol(geno_pheno_full_table_sig)]) {
                
                pheno_data_array = geno_pheno_full_table_sig[,1]
                
                current_geno_data = geno_pheno_full_table_sig[,col]
                
                anova_table = cbind(pheno_data_array, current_geno_data)
                
                anova_table = anova_table[complete.cases(anova_table),]
                
                initial_anova = aov(anova_table[,1] ~ anova_table[,2])
                
                # Rijan: I am not sure why this needs to be modified any further than it already is but I might be wrong
                n_sample_df <- as.data.frame(table(anova_table[,2]))
                n_sample_df_m <- t(n_sample_df[,2])
                colnames(n_sample_df_m) <- paste("n_",t(n_sample_df[,1]), sep="")
                
                Tukey_table <- TukeyHSD(initial_anova,ordered = TRUE)
                
                pg_dt_tukey <- as.data.frame(Tukey_table$'anova_table[, 2]')
                temp_sig_pheno_db <- pg_dt_tukey
                temp_sig_pheno_db[,1:3] <- NULL
                current_sig_pheno_db <- as.data.frame(cbind(gene, col, t(temp_sig_pheno_db), n_sample_df_m))
                
                # accumulate the data present in current_sig_pheno_db
                # This is need to formalize the counting
                names(current_sig_pheno_db) <- gsub(x = names(current_sig_pheno_db), pattern = "Alt/Alt-Ref/Ref", replacement = "Ref/Ref-Alt/Alt")  
                names(current_sig_pheno_db) <- gsub(x = names(current_sig_pheno_db), pattern = "Ref/Alt-Ref/Ref", replacement = "Ref/Ref-Ref/Alt")  
                names(current_sig_pheno_db) <- gsub(x = names(current_sig_pheno_db), pattern = "Alt/Alt-Ref/Alt", replacement = "Ref/Alt-Alt/Alt")
                
                variant_sig_db <- plyr::rbind.fill(variant_sig_db, current_sig_pheno_db)
            }
            
        }
    }
    
    # Rijan: Placeholder for accumulating data down the road
    odat <- data.frame()
    
    # Rijan: This loop will tally and accumlate the EFFECT and TYPE fields of the output_db table, which itself holds the upstream data for the gene.

    for( ID_sample in genes_of_interest){
      
      sub_db <- output_db %>% filter(ID == ID_sample)
      gene <- t(as.data.frame(table(sub_db$EFFECT)))
      colnames(gene) <- gene[1,]
      gene <- gene[-c(1),]
      gene <- t(as.data.frame(gene))
      
      types <- t(as.data.frame(table(sub_db$TYPE)))
      colnames(types) <- types[1,]
      types <- types[-c(1),]
      types <- t(as.data.frame(types))
      
      temp_df <- cbind(ID_sample, gene, types)
      odat <- rbind(odat, temp_df)
    }
    
    # Rijan: Another loop to wrangle the data for regex 
    colnames(variant_sig_db)[2] <- "SNP_target"
    variant_sig_db$`Ref/Ref-Alt/Alt` <- as.numeric(as.character(variant_sig_db$`Ref/Ref-Alt/Alt`))
    variant_sig_db$`Ref/Ref-Ref/Alt` <- as.numeric(as.character(variant_sig_db$`Ref/Ref-Ref/Alt`))
    variant_sig_db$`Ref/Alt-Alt/Alt` <- as.numeric(as.character(variant_sig_db$`Ref/Alt-Alt/Alt`))
    
    odat_sig_test_all <- data.frame()
    
    # Rijan: combine `odat`s to make a final `odat_sig_test_all`
    
    for(gene_id in genes_of_interest){
    per_gene_variant_sig_db <- subset(variant_sig_db, subset = gene == gene_id)
    temp_min_sig_variant <- variant_sig_db[which.min(variant_sig_db$`Ref/Ref-Alt/Alt`),]
    odat_sig_test_all <- rbind(odat_sig_test_all, temp_min_sig_variant)
    }
     
    odat <- merge(odat, odat_sig_test_all, by.x = 1, by.y = 1, all.x = T, all.y = F)
    odat <- unique.data.frame(odat)
    
}

# %% [code]

