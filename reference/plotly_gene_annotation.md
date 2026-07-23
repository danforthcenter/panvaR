# Plot genes and their locations

Use a table of genes, descriptions and their start and end points to
plot their locations and annotations along the genome.

## Usage

``` r
plotly_gene_annotation(
  panvar.table.list = NULL,
  annotation.table = NULL,
  middle.snp = NULL,
  window,
  gene.color = "blue",
  point.color = NULL,
  point.fill.scale = NULL
)
```

## Arguments

- panvar.table.list:

  list, output from
  [make_panvar_tables](https://danforthcenter.github.io/panvaR/reference/make_panvar_tables.md).
  Provide either this list or both annotation.table and middle.snp.

- annotation.table:

  table with annotations with columns (geneID, CHR, start, end,
  annotation). start and end correspond to base-pair coordinates of
  start and end of gene. CHR is chromosome of gene.

- middle.snp:

  character, SNP name in form "CHR-POS" center of window. Often the
  key.snp output of
  [get_ld_in_window](https://danforthcenter.github.io/panvaR/reference/get_ld_in_window.md)

- window:

  integer, kilobases on either side of middle.snp to plot

- gene.color:

  character, color to plot genes

- point.color:

  character, variable in annotation.table that indicates how to color
  points plotted next to gene descriptions. If not supplied, no points
  are plotted. The input "LD" is reserved to give functionality to
  [plot_panvar](https://danforthcenter.github.io/panvaR/reference/plot_panvar.md).
  If used, legend will not be displayed.

- point.fill.scale:

  ggplot2 scale object, a fill scale to customize how point.color is
  displayed.

## Value

GGplot object

## Examples

``` r
#' # organize options
tag.snp <- "Chr_05-6857045"
gwas.df <- read.csv(system.file(
    "extdata",
    "PanvarExample_GLM_GWASresults.csv",
    package = "panvaR"))
annotation.table <- read.csv(system.file(
    "extdata",
    "Setaria_shattering_annotation.csv",
    package = "panvaR"))
plink.path <- bigsnpr::download_plink2()
temp.dir <- file.path(tempdir(), "panvar_ex")
dir.create(temp.dir, showWarnings = FALSE)
geno.bed.filename <- "Setaria_shattering_example_pruned.bed"
geno.bed.directory <- system.file("extdata", package="panvaR")

# make input tables
tables <- make_panvar_tables(
  gwas.res = gwas.df,
  tag.snp = tag.snp,
  annotation.table = annotation.table,
  plink.path = plink.path,
  pvals.in.log = F,
  geno.bed.filename = geno.bed.filename,
  geno.bed.directory = geno.bed.directory,
  window = 25,
  temp.dir = temp.dir,
  compute.scores = FALSE,
  snp.to.gene.buffer = 0)
#> Calculating LD
#> Generating snp to gene correspondence
  
# plot
plotly_gene_annotation(
  panvar.table.list = tables,
  window = 25)
#> Warning: Ignoring unknown aesthetics: text
#> Warning: Ignoring unknown aesthetics: text

{"x":{"data":[{"x":[0.5,0.5,null,0.5,0.5,null,0.5,0.5,null,0.5,0.5,null,0.5,0.5,null,0.5,0.5,null,0.5,0.5],"y":[-6829932,-6832531,null,-6837639,-6838969,null,-6847970,-6850236,null,-6859612,-6862290,null,-6862288,-6866177,null,-6866196,-6868255,null,-6867108,-6873803],"text":["Sevir.5G085300","Sevir.5G085300",null,"Sevir.5G085350","Sevir.5G085350",null,"Sevir.5G085400","Sevir.5G085400",null,"Sevir.5G085500","Sevir.5G085500",null,"Sevir.5G085600","Sevir.5G085600",null,"Sevir.5G085700","Sevir.5G085700",null,"Sevir.5G085800","Sevir.5G085800"],"type":"scatter","mode":"lines","line":{"width":7.559055118110237,"color":"rgba(0,0,255,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.5,0.55000000000000004,null,0.5,0.55000000000000004,null,0.5,0.55000000000000004,null,0.5,0.55000000000000004,null,0.5,0.55000000000000004,null,0.5,0.55000000000000004,null,0.5,0.55000000000000004],"y":[-6831231.5,-6836420,null,-6838304,-6843295,null,-6849103,-6850170,null,-6860951,-6857045,null,-6864232.5,-6863920,null,-6867225.5,-6870795,null,-6870455.5,-6877670],"text":"","type":"scatter","mode":"lines","line":{"width":1.8897637795275593,"color":"rgba(0,0,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.55500000000000005,0.55500000000000005,0.55500000000000005,0.55500000000000005,0.55500000000000005,0.55500000000000005,0.55500000000000005],"y":[-6836420,-6843295,-6850170,-6857045,-6863920,-6870795,-6877670],"text":["(1 of 2) PTHR20961//PTHR20961:SF22 - GLYCOSYLTRANSFERASE // SUBFAMILY NOT NAMED","No gene description.","(1 of 1) PTHR10641//PTHR10641:SF510 - MYB-LIKE DNA-BINDING PROTEIN MYB // SUBFAMILY NOT NAMED","(1 of 2) PTHR33146:SF2 - ENDONUCLEASE 2","(1 of 2) PTHR33146:SF2 - ENDONUCLEASE 2","(1 of 1) PTHR34543//PTHR34543:SF1 - FAMILY NOT NAMED // ABSCISIC ACID (ABA)-DEFICIENT 4 PROTEIN","(1 of 1) KOG4467 - Uncharacterized conserved protein"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(127,127,127,1)","opacity":1,"size":11.338582677165356,"symbol":"square","line":{"width":1.8897637795275593,"color":"rgba(127,127,127,1)"}},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.56000000000000005,0.56000000000000005,0.56000000000000005,0.56000000000000005,0.56000000000000005,0.56000000000000005,0.56000000000000005],"y":[-6836420,-6843295,-6850170,-6857045,-6863920,-6870795,-6877670],"text":["Sevir.5G085300","Sevir.5G085350","Sevir.5G085400","Sevir.5G085500","Sevir.5G085600","Sevir.5G085700","Sevir.5G085800"],"hovertext":["","","","","","",""],"textfont":{"size":14.611872146118722,"color":"rgba(0,0,0,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","textposition":"right","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":10.958904109589042,"l":63.561643835616451},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.495,0.60499999999999998],"tickmode":"array","ticktext":["0.500","0.525","0.550","0.575","0.600"],"tickvals":[0.5,0.52500000000000002,0.55000000000000004,0.57499999999999996,0.59999999999999998],"categoryorder":"array","categoryarray":["0.500","0.525","0.550","0.575","0.600"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.6529680365296811,"tickwidth":0,"showticklabels":false,"tickfont":{"color":null,"family":null,"size":0},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-6887295,-6826795],"tickmode":"array","ticktext":["-25 KB","-18.75 KB","-12.5 KB","-6.25 KB","0 KB","6.25 KB","12.5 KB","18.75 KB","25 KB"],"tickvals":[-6832045,-6838295,-6844545,-6850795,-6857045,-6863295,-6869545,-6875795,-6882045],"categoryorder":"array","categoryarray":["-25 KB","-18.75 KB","-12.5 KB","-6.25 KB","0 KB","6.25 KB","12.5 KB","18.75 KB","25 KB"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"1904731f2a5b":{"x":{},"y":{},"yend":{},"xend":{},"text":{},"type":"scatter"},"190434165778":{"x":{},"y":{},"xend":{},"yend":{}},"19043f55bec3":{"x":{},"y":{},"text":{}},"190428358d62":{"x":{},"y":{},"label":{}}},"cur_data":"1904731f2a5b","visdat":{"1904731f2a5b":["function (y) ","x"],"190434165778":["function (y) ","x"],"19043f55bec3":["function (y) ","x"],"190428358d62":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}  
# clean up 
unlink(temp.dir, recursive = TRUE)
```
