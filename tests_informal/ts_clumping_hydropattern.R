setwd("~/Projects/Hydropattern_finemap/")

hits <- read.csv("data/20260301_Hydropatterning_GWAS_sigSNP.csv")

hits <- hits %>% 
  mutate(marker.ID = paste0(Chromosome, "-", Position_B73v5)) %>% 
  rename(PVAL = p.value,
         CHR = Chromosome,
         POS = Position_B73v5) %>% 
  arrange(CHR, POS) %>% 
  mutate(hit.ID.num = 1:n()) %>% 
  relocate(hit.ID.num)

# x <- bigsnpr::snp_readBed("data/Schnable2023_V5_allchrs_hydropatterngenos.bed")
# geno <- bigsnpr::snp_attach(x)
snps <- data.table::fread("data/Schnable2023_V5_allchrs_hydropatterngenos.bim", data.table = F)
names(snps) <- c("chromosome", "marker.ID", "genetic.dist", "physical.pos", "allele1", "allele2")     


# get proxy snps 
for(i in 1:nrow(hits)){
  this.chr <- hits$CHR[i]
  this.pos <- hits$POS[i]
  
  this.df <- snps %>% 
    filter(chromosome == this.chr) %>% 
    # filter(physical.pos != this.pos) %>% 
    mutate(dist = abs(physical.pos - this.pos))
  
  out.marker <- this.df$marker.ID[which.min(this.df$dist)]
  
  this.out <- data.frame(hit.ID.num = hits$hit.ID.num[i],
                         closest.dist = this.df$dist[which.min(this.df$dist)],
                         closest.marker = out.marker)
  
  dist <- min(abs(this.df$POS - this.pos))
  if(i == 1){
    all.out <- this.out
  } else {
    all.out <- bind_rows(this.out, all.out)
  }
}


hits.proxy <- hits %>% 
  left_join(all.out, by = "hit.ID.num")

df.clump.input <- hits.proxy %>% 
  select(closest.marker, PVAL) %>% 
  mutate(CHR = get_chrom_from_id(closest.marker),
         POS = get_bp_from_id(closest.marker)) %>% 
  rename(marker.ID = closest.marker)


clumps <-
  panvaR::snp_make_clumps(geno.bed.filename = "Schnable2023_V5_allchrs_hydropatterngenos",
                          geno.bed.dir = "data/",
                          gwas.res = df.clump.input,
                          pvals.in.log = F,
                          window = 50,
                          ld.thresh = .6,
                          plink.path = normalizePath("~/bin/plink2"))

# problem snp: 1-43628780
# this doesn't report any non-na ld to any snps
# I think that this particular snp is of poor quality (no variance)
# I changed the flag "ld-r2" to -1 instead of 0 based on this:
# https://groups.google.com/g/plink2-users/c/zD2aME9Kb3o/m/3_yUsDQhAAAJ
# 

# plink2 --bfile Schnable2023_V5_allchrs_hydropatterngenos --r2-unphased --ld-snp 1-43628780 --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0 --out TEMPLD
# plink2 --bfile Schnable2023_V5_allchrs_hydropatterngenos --r2-unphased --ld-snp 3-112056458 --ld-window-kb 50 --ld-window 99999 --ld-window-r2 0 --out TEMPLD

3-112056458

current_args <- c(
  "--silent", 
  "--bfile",
  bfile,
  "--r2-unphased",
  "--ld-snp",
  snp.name,
  "--ld-window-kb",
  window,
  "--ld-window",
  "99999",
  "--ld-window-r2",
  "0",
  "--out",
  outfile
)