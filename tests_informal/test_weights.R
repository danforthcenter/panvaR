
set_plink_path("/home/cluebbert/bin/plink2")
set_panvar_prefix("SetShattering")
set_out_dir("~/scratch/panvar_test/setaria_shatter_full")

out_ld <- 
  panvaR::get_ld_in_window(tag.snp = "5-6857045",
                           window = 500,
                           in.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full/",
                           out.dir = "/home/cluebbert/scratch/panvar_test/setaria_shatter_full",
                           geno.bed = "SetShattering_PlinkQC_maf0.05_missing0.1")

gwas.res <- read.csv("~/scratch/panvar_test/setaria_shatter_full/Shatter_GLM_res_with_snpeff_and_cad.csv")
ld.list <- out_ld
middle.snp <- out_ld$key.snp
window <- 25

# start function testing
gwas.sub <- gwas.res %>%
  as.data.frame() %>% 
  mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
  left_join(ld.list$table, by = "marker.ID") %>%
  filter(!is.na(.data$R2)) %>% 
  rename("LD" = "R2")
# add middlesnp back in
gwas.sub.mid.snp <- gwas.res %>% 
  filter(marker.ID == middle.snp) %>% 
  mutate(LD = 1)
gwas.sub <- bind_rows(gwas.sub, gwas.sub.mid.snp)

if(!pvals.in.log){
  gwas.sub <- gwas.sub %>%
    mutate(LOGPVAL = -log10(.data$PVAL))
} else {
  gwas.sub <- gwas.sub %>% 
    rename("LOGPVAL" = "PVAL")
}

score.in <- gwas.sub %>% 
  mutate(key.snp.pos = get_bp_from_id(middle.snp)) %>% 
  mutate(DIST = abs(POS - key.snp.pos)) %>% 
  # only need this for testing
  select(marker.ID, POS, DIST, LOGPVAL, LD, zero_shot)
  
# make a function that takes an arbitrary number of variables and their directions
# and outputs the min/max scaled 

minmaxnormal <- function(x) {
  return((x- min(x)) /(max(x)-min(x)))
}

cols <- c("DIST", "LOGPVAL", "LD", "zero_shot")
input.df <- score.in
directions <- c(-1, 1, 1, -1)
weights <- c(0, 1/3, 1/3, 1/3)
make_weights <- function(input.df, directions, weights = NULL, cols){}

if(is.null(weights)){
  weights <- rep.int(1, length(cols)) / length(cols)
}
dirs.df <- data.frame(var = cols, direction = directions, weight = weights)
out <- input.df %>% 
  select(marker.ID, all_of(cols)) %>% 
  tidyr::pivot_longer(cols = -marker.ID, names_to = "var", values_to = "value") %>% 
  left_join(dirs.df, by = "var") %>% 
  rowwise() %>% 
  mutate(value.direction = value * direction) %>% 
  group_by(var) %>% 
  mutate(value.normalized = minmaxnormal(value.direction)) %>% 
  group_by(marker.ID) %>% 
  mutate(weighted.score = weighted.mean(value.normalized, weight)) %>%
  # stop pipe here to test weighting
  ungroup() %>% 
  select(marker.ID, weighted.score) %>% 
  distinct()


# did the direction make sense?
ggplot(aes(x = value, y = value.normalized), data = x) +
  geom_point() +
  facet_wrap(~var, scales = "free")


# test weighting a little, see where to break pipe above
test <- x %>% 
  filter(marker.ID == "5-6357125") 

mean(test$value.normalized)

test2 <- x %>% 
  filter(marker.ID == "5-6357125") %>% 
  filter(var != "DIST")
mean(test2$value.normalized)