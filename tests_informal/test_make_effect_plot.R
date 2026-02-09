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
make_effect_plot(
  panvar.table.list = tables,
  # gwas.res = gwas.df,
  # ld.list = out_ld,
  pvals.in.log = F,
  window = 8,
  sig.line = 6
)



make_panvar_manhattan(
  panvar.table.list = tables,
  pvals.in.log = F,
  window = 8, 
  sig.line = 6,
  quantitative.annotation = "LD"
)

# manhattan
make_panvar_manhattan(gwas.res = gwas.df,
                      pvals.in.log = F,
                      ld.list = out_ld,
                      window = 10,
                      sig.line = 6, orient = "V")


# ------------------------------------------------------------------------\
# why different fill scales man/eff --------
# ------------------------------------------------------------------------\

# saved some dataframes to see if they were identical
# turned out I was filtering the plot.df in the effect but using the plot limits 
# to filter in the manhattan. This was leading to different fill scales being used.
# I made it so now both filter for only snps in the window when plotting 
effplot.df <- read.csv("./effect_test.csv")
manplot.df <- read.csv("./man_test.csv") %>% 
  filter(!is.na(plot.R2))
