set_out_dir("~/scratch/mvp_ts/")

panvaR::vcf_to_plink2("/home/cluebbert/scratch/setaria_annotated_vcf.vcf.vcf")


rMVP::MVP.Data(fileBed = "~/scratch/mvp_ts/setaria_annotated_vcf",
               filePC = TRUE,
               pcs.keep = 10,
               filePhe = pheno,
               out = "~/scratch/mvp_ts/SetariaOrigChrom")

pheno <- read.delim("inst/extdata/Setaria_shattering_example_phenotype.tsv")
pheno <- read.table("~/scratch/mvp_ts/SetariaOrigChrom.phe", header = T)
geno <- attach.big.matrix("~/scratch/mvp_ts/SetariaOrigChrom.geno.desc")
map <- read.table("~/scratch/mvp_ts/SetariaOrigChrom.geno.map", header = T)

MVP(phe = pheno,
    geno = geno,
    map = map,
    method = "GLM",
    nPC.GLM = 3,
    maf = NULL,
    file.output = c("pmap", "plot"), # consider setting this to false and parsing the output a little more myself 
    outpath = "~/scratch/mvp_ts/",
    memo = "SetariaOrigChrom")
