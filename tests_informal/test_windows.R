# test on windows 



# set global options 
set_plink_path("C:/Users/cluebbert/Documents/bin/plink2.exe")
set_panvar_prefix("SetShattering")
set_out_dir("C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/Panvar_TS/testing_data/")


# ------------------------------------------------------------------------\
# test sub modules --------
# ------------------------------------------------------------------------\

# ld
out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "Chr_05-6857045",
                           window = 50,
                           in.dir = "./inst/extdata",
                           out.dir = options()$panvar_outdir,
                           geno.bed = "Setaria_shattering_example_pruned")

# ------------------------------------------------------------------------\
# test windows plink --------
# ------------------------------------------------------------------------\

plink2 --bfile 'C:\Users\cluebbert\OneDrive - DDPSC\~Rprojects~\panvaR\inst\extdata\Setaria_shattering_example_pruned.bed' `
--r2-unphased --ld-snp 5-6857045 --ld-window-kb 50 --ld-window-r2 -1 --out 'C:/Users/cluebbert/OneDrive - DDPSC/~Rprojects~/Panvar_TS/testing_data/ld_out_temp'


# need to some how quote file paths because windows allows spaces
test.dir <- "some/directory/toquote"
x <- paste0("\'", test.dir, "\'")
x
