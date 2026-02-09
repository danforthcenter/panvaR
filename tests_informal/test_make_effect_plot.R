# set global options 
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")

# ld
out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "5-6857045",
                           window = 500,
                           in.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                           out.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full",
                           geno.bed = "SetShattering_PlinkQC_maf0.05_missing0.1")

gwas.df <- data.table::fread(file.path(options()$panvar_outdir, "SetShattering_GLM_GWASresults.csv"))

# effect
make_effect_plot(gwas.res = gwas.df, 
                 pvals.in.log = F,
                 ld.list = out_ld,
                 window = 10,
                 sig.line = 6)

# manhattan
make_panvar_manhattan(gwas.res = gwas.df,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 10,
                      sig.line = 6, orient = "V")