#' ---
#' title: "GSEA"
#' author: "Nicolas Delhomme for Bn Bioinformatics"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'   fig_width: 9
#'   fig_height: 6
#'   toc: TRUE
#'   number_sections: TRUE
#'   toc_depth: 3
#'   toc_float:
#'    collapsed: TRUE
#'    smooth_scroll: TRUE
#'   code_folding: hide
#'   theme: "flatly" 
#' ---
#' 
#' # Synopsis
#' Reproduce the analysis from the [Alhamdoosh _et al._, 2017](https://academic.oup.com/bioinformatics/article/33/3/414/2875813) manuscript.
#' 
#' # Setup
#' 
#' * Libraries and Helpers
#' 
#' Dynamic programming to install the packages if needed
#' ```{r test-install, echo=TRUE, eval=FALSE}
#' lapply(c("limma","edgeR","EGSEA","EGSEAdata"),function(l){
#'   if(!eval(parse(text=paste("require(",l,")")))){
#'     BiocManager::install(l)
#'   }
#' })
#'```
#'

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(EGSEA)
  library(EGSEAdata)
})

#' # Data
#' ## Count data
#' The URL for the data is listed in the [_F1000 Research_ tutorial](https://f1000research.com/articles/6-2010/v1) associated to the manuscript aforementioned.
#' 
#' We can directly load it using the url 
load(url("http://bioinf.wehi.edu.au/EGSEA/mam.rnaseq.rdata"))

#' A quick peek
names(mam.rnaseq.data)
dim(mam.rnaseq.data)

#' ### Normalisation and design
#' * library sequencing depth
x = calcNormFactors(mam.rnaseq.data, method = "TMM")

#' * experimental design
design = model.matrix(~0+x$samples$group+x$samples$lane)
colnames(design) = gsub("x\\$samples\\$group", "", colnames(design))
colnames(design) = gsub("x\\$samples\\$lane", "", colnames(design))
head(design)

#' * contrasts
contr.matrix = makeContrasts(
  BasalvsLP = Basal-LP,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design))
head(contr.matrix)

#' * normalisation (_voom_)
v = voom(x, design, plot=FALSE)
names(v)

#' ## EGSEA annotations
#' ### Retrieving the mouse data 
egsea.data("mouse")
info = egsea.data("mouse", returnInfo = TRUE)

#' * a quick peek
names(info)
info$msigdb$info$collections

#' ### build and index a selection
gs.annots = buildIdx(entrezIDs=v$genes$ENTREZID, species="mouse",
                     msigdb.gsets=c("c2", "c5"), go.part = TRUE)

#' * a quick peek
names(gs.annots)
class(gs.annots$c2)
summary(gs.annots$c2)
show(gs.annots$c2)

#' * Looking into the details
s = getSetByName(gs.annots$c2, "SMID_BREAST_CANCER_LUMINAL_A_DN")

class(s)
names(s)
names(s$SMID_BREAST_CANCER_LUMINAL_A_DN)
slotNames(gs.annots$c2)
colnames

#' Get the gene ID - gene symbol mapping
symbolsMap = v$genes[, c(1, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

#' # EGSEA
#' ## Preparation
egsea.base()
baseMethods = egsea.base()[-2]
baseMethods
egsea.combine()
egsea.sort()

#' * Run it! It takes about five minutes to complete.
gsa = egsea(voom.results=v, contrasts=contr.matrix,
            gs.annots=gs.annots, symbolsMap=symbolsMap,
            baseGSEAs=baseMethods, sort.by="med.rank",
            num.threads = 8, report = FALSE)

#' * a quick peek
show(gsa)
summary(gsa)

#' ## Results
topSets(gsa, gs.label="c2", contrast = "comparison", names.only=TRUE)
t = topSets(gsa, contrast = "comparison",
            names.only=FALSE, number = Inf, verbose = FALSE)
t[grep("LIM_", rownames(t)), c("p.adj", "Rank", "med.rank", "vote.rank")]

topSets(gsa, gs.label="kegg", contrast="BasalvsLP", sort.by="med.rank")
topSets(gsa, gs.label="kegg", contrast="comparison", sort.by="med.rank")

#' ## Visualisation
plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_UP", gs.label="c2",
            contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_UP")

plotHeatmap(gsa, gene.set="LIM_MAMMARY_STEM_CELL_DN", gs.label="c2",
            contrast = "comparison", file.name = "hm_cmp_LIM_MAMMARY_STEM_CELL_DN")

plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
            contrast = "BasalvsLP", gs.label = "kegg",
            file.name = "Vascular_smooth_muscle_contraction")


plotPathway(gsa, gene.set = "Vascular smooth muscle contraction",
            contrast = "comparison", gs.label = "kegg",
            file.name = "Vascular_smooth_muscle_contraction_cmp")

plotMethods(gsa, gs.label = "c2", contrast = "BasalvsLP", file.name = "mds_c2_BasalvsLP")
plotMethods(gsa, gs.label = "c5BP", contrast = "BasalvsLP", file.name = "mds_c5_BasalvsLP")

plotSummary(gsa, gs.label = 3, contrast = 3, file.name = "summary_kegg_LPvsML")
            
plotSummary(gsa, gs.label = 1, contrast = 3, file.name = "summary_c2_LPvsML",
            x.axis = "med.rank")

plotSummary(gsa, gs.label = 1, contrast = 3, file.name = "summary_sig_c2_LPvsML", x.axis = "med.rank", x.cutoff=300)

plotSummary(gsa, gs.label = "kegg", contrast = c(1,2),file.name = "summary_kegg_1vs2")

plotGOGraph(gsa, gs.label="c5BP", contrast = 1, file.name="BasalvsLP-c5BP-top-")

plotGOGraph(gsa, gs.label="c5CC", contrast = 1, file.name="BasalvsLP-c5CC-top-")

plotBars(gsa, gs.label = "c2", contrast = "comparison", file.name="comparison-c2-bars")

plotSummaryHeatmap(gsa, gs.label="c2", hm.vals = "avg.logfc.dir", file.name="summary_heatmaps_c2")

plotSummaryHeatmap(gsa, gs.label="kegg", hm.vals = "avg.logfc.dir", file.name="summary_heatmaps_kegg")

#' ## Candidate genes
t = limmaTopTable(gsa, contrast=1)
head(t)

#' # Report
#' This takes a long time...
#' ```{r report, eval=FALSE}
#' generateReport(gsa, number = 20, report.dir="./mam-rnaseq-egsea-report")
#' ```
#'
#' # Session Info
#' ```{r session-info, echo=FALSE}
#' sessionInfo()
#' ```
#'
