# still haven't fixed the non numeric chromsome ids in the example data so this is the id
# "Chr_05-6857045"
# Setaria_shattering_example_pruned

set_plink_path("/home/cluebbert/bin/plink2")

# test tag snp
out <- 
get_ld_in_window(tag.snp = "Chr_05-6857045",
                 window = 500,
                 in.dir = "inst/extdata/",
                 out.dir = "~/scratch/panvar_test/",
                 geno.bed = "Setaria_shattering_example_pruned")

# test qtl.df

