library(TrenaProjectErythropoiesis)
library(igvR)
library(later)
library(trenaSGM)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
library (RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
currentColorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
if(!exists("state"))
   state <- new.env(parent=emptyenv())
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")){
   tp <- TrenaProjectErythropoiesis()
   setTargetGene(tp, "GATA2")
   }

if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   later(function() showGenomicRegion(igv, "GATA2"), 3)
   #later(function(){
   #         track <- DataFrameQuantitativeTrack("gh", tbl.enhancers[, c("chrom", "start", "end", "combinedScore")],
   #                                             autoscale=FALSE, color="blue", min=0, max=50)
   #         displayTrack(igv, track)
   #         }, 4)
   } # igv
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function(chromosome, start.loc, end.loc)
{
   directory <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
display_ATACseq_in_enhancers <- function(tbl.enhancers=data.frame(), trackName="atac")
{
   chrom.loc <- getGenomicRegion(igv)
   if(nrow(tbl.enhancers) > 0)
      gr.region <- GRanges(tbl.enhancers.myeloid[, c("chrom", "start", "end")])
   else
      gr.region <- with(chrom.loc, GRanges(seqnames=chrom, IRanges(start, end)))

   tbl.atac <- with(chrom.loc, getATACseq(chrom=chrom, start=start, end=end))
   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.atac), gr.region, type="any"))
   tbl.atac.trimmed <- tbl.atac[tbl.ov$queryHits,]
   dim(tbl.atac.trimmed)
   samples <- unique(tbl.atac.trimmed$sample)
    # if the area above has enough span, then open chromatin is found in 11 samples:
    # "d04_rep1" "d04_rep2" "d08_rep1" "d10_rep1" "d10_rep2" "d11_rep1" "d11_rep2" "d12_rep1" "d12_rep2" "d16_rep1" "d16_rep2"
   atac.regions.by.sample <- list()
   for(sample.x in samples){
      tbl.sample <- subset(tbl.atac.trimmed, sample==sample.x)[, c("chrom", "start", "end")]
      atac.regions.by.sample[[sample.x]] <- tbl.sample
      dim(tbl.sample)
      currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
      color <- colors[currentColorNumber]
      track <- DataFrameAnnotationTrack(sample.x, tbl.sample, color=color)
      displayTrack(igv, track)
      #write.table(tbl.sample, file=sprintf("tbl.%s.bed", sample.x), quote=FALSE, row.names=FALSE)
      } # for sample.x

   tbl.regions.condensed <- as.data.frame(union(GRanges(tbl.atac.trimmed[, c("chrom", "start", "end")]),
                                                GRanges(tbl.atac.trimmed[, c("chrom", "start", "end")])))[, c("seqnames", "start", "end")]

   colnames(tbl.regions.condensed) <- c("chrom", "start", "end")
   tbl.regions.condensed$chrom <- as.character(tbl.regions.condensed$chrom)
   track <- DataFrameAnnotationTrack("atac combined", tbl.regions.condensed, color="black")
   displayTrack(igv, track)

   state$atac.regions.by.sample <- atac.regions.by.sample

} # display_ATACseq_inEnhancers
#------------------------------------------------------------------------------------------------------------------------
test_display_ATACseq_in_enhancers <- function()
{
   message(sprintf("--- test_display_ATACseq_in_enhancers"))
   showGenomicRegion(igv, "chr3:128,602,838-128,611,664")
   setTargetGene(tp, "GATA2")
   all.tissues <- getEnhancerTissues(tp)
   myeloid.tissues <- grep("myeloid", all.tissues, ignore.case=TRUE, value=TRUE)

   tbl.enhancers <- getEnhancersForTissues(myeloid.tissues)
   dim(tbl.enhancers)
   display_ATACseq_in_enhancers(tbl.enhancers)

   display_ATACseq_in_enhancers()


} # test_display_ATACseq_in_enhancers
#------------------------------------------------------------------------------------------------------------------------
makeAtacBasedModels <- function()
{
  stopifnot("atac.regions.by.sample" %in% names(state))
  models <- list()
  samples <- names(state$atac.regions.by.sample)
  for(sample.x in samples[1:2]){
     printf("--- build model with sample %s", sample.x)
     tbl.regions <- state$atac.regions.by.sample[[sample.x]]
     models[[sample.x]] <- buildModel(tbl.regions)
     } # for sample.x

 state$models <- models

} # makeAtacBasedModels
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(tbl.regions)
{
   genome <- "hg38"
   tss <- getTranscriptsTable(tp)$tss[1]
   mtx.name <- "brandLabDifferentiationTimeCourse-27171x28"
   stopifnot(mtx.name %in% getExpressionMatrixNames(tp))
   mtx <- getExpressionMatrix(tp, mtx.name)
   dim(mtx)

   #db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs", "gata2.gh.fimoBindingSites.sqlite")
   db.name <- "fimoResults-10e-3-chr3-128383794-128647775.sqlite"
   stopifnot(file.exists(db.name))

   recipe <- list(title="fimo.atacseq",
                  type="fimo.database",
                  regions=tbl.regions,
                  geneSymbol=getTargetGene(tp),
                  tss=tss,
                  matrix=mtx,
                  db.host="localhost",
                  db.port=NA_integer_,
                  databases=list(db.name),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="fimoDatabase",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.4,
                  maxModelSize=100,
                  fimoPValueThreshold=1e-5,,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- FimoDatabaseModelBuilder(genome, getTargetGene(tp),  recipe, quiet=TRUE)
   x <- build(builder)

   x

} # buildModel
#------------------------------------------------------------------------------------------------------------------------
displayMotifsForTF <- function(tf="TBX15")
{
  #pfms <- query(MotifDb, c("hsapiens", tf), c("swissregulon", "jaspar2018", "hocomoco")) # 537
  pfms <- query(MotifDb, c("hsapiens", tf), c("jaspar2018")) # 537
  mm <- MotifMatcher("hg38", as.list(pfms), quiet=TRUE)
  loc <- getGenomicRegion(igv)
  tbl.loc <- with(loc, data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))

  tbl.matches <- findMatchesByChromosomalRegion(mm, tbl.loc, pwmMatchMinimumAsPercentage=85)
   # consensusString(pfms[[1]])  [1] "AGGTGTGA"
   #  consensusString(pfms[[2]]) [1] "AGGTGTGA"
   #  consensusString(pfms[[3]]) [1] "GGGGGGGGG?GGGTGGG??"
  browser()
  dim(tbl.matches)
  trackName <- sprintf("%s", tf)
  track <- DataFrameQuantitativeTrack(trackName, tbl.matches[, c("chrom", "motifStart", "motifEnd", "motifRelativeScore")],
                                      color="red", autoscale=FALSE, min=0, max=1, trackHeight=25)
  displayTrack(igv, track)

} # displayMotifsForTF
#------------------------------------------------------------------------------------------------------------------------
# see findLowerFidelityBindingSites function in
#   ~/github/TrenaProjectErythropoiesis/prep/gata-switch/tbx15.R
#  tbl.matches <- findMatchesByChromosomalRegion(mm, state$tbl.regions.condensed, 85)
#
# with these 5 hits in marjorie's atac-seq regions
#   motifName chrom motifStart  motifEnd strand motifScore motifRelativeScore    match chromStart  chromEnd                  seq status shortMotif
#   2  Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128475449 128475456      +   6.710003          0.8856179 GGGTGTGA  128475407 128475590 TCCCTAAGAGTTGCAGA...     wt   MA0803.1
#   12 Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128486955 128486962      +   6.661377          0.8792000 TGGTGTGA  128486844 128487055 GGGTTGGCATAGTAGGG...     wt   MA0803.1
#   21 Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128483451 128483458      +   6.658030          0.8787582 CGGTGTGA  128482966 128483521 AGAAGGGACAGAGGGAC...     wt   MA0803.1
#   13 Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128497527 128497534      +   6.658030          0.8787582 CGGTGTGA  128497484 128498222 GGGTGGACTCCGGGCTG...     wt   MA0803.1
#   11 Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128483002 128483009      +   6.642737          0.8767398 AGGGGTGA  128482966 128483521 AGAAGGGACAGAGGGAC...     wt   MA0803.1
#   1  Hsapiens-jaspar2018-TBX15-MA0803.1  chr3  128461602 128461609      +   6.636363          0.8758986 AGGCGTGA  128461467 128461684 CACCACGCCCGGCTAAT...     wt   MA0803.1
#
#      and this igv plot:
#         [[file:/Users/paul/images/gata2-tbx15-bindingSites.png]]]]
#
#    and this expression correlation plot:
#       [[file:/Users/paul/images/gata2-tbx15-expression-correlation.png]]]]
#
hopeRestored <- function()
{
   tbl.hope.restored <- data.frame(chrom=rep("chr3", 6),
                                   start=c(128475449, 128486955, 128483451, 128497527, 128483002, 128461602),
                                   end=c(128475456, 128486962, 128483458, 128497534, 128483009, 128461609),
                                   stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack("TBX15/GeneHancer/atac", tbl.hope.restored, color="blue", trackHeight=30)
   displayTrack(igv, track)


} # hopeRestored
#------------------------------------------------------------------------------------------------------------------------
# uses targetGene already set for TrenProjectErythropoiesis
getEnhancersForTissues <- function(tissues, display=FALSE, trackName=NA)
{
  if(!tissues[1] == "all")
     stopifnot(all(tissues %in% getEnhancerTissues(tp)))

  tbl.enhancers <- getEnhancers(tp, tissues=tissues)

  if(display){
     currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
     color <- colors[currentColorNumber]
     track.name <- "gh.all"
     if(!is.na(trackName))
        track.name <- trackName
     track <- DataFrameQuantitativeTrack(track.name, tbl.enhancers[, c("chrom", "start", "end", "combinedscore")],
                                         autoscale=TRUE, color=color)
     displayTrack(igv, track)
     } # display

   tbl.enhancers

} # getEnhancersForTissues
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancersForTissues <- function()
{
   message(sprintf("--- test_getEnhancersForTissues"))

   setTargetGene(tp, "GATA2")
   all.tissues <- getEnhancerTissues(tp)
   myeloid.tissues <- grep("myeloid", all.tissues, ignore.case=TRUE, value=TRUE)

   tbl.enhancers.myeloid <- getEnhancersForTissues(tissues=myeloid.tissues)
   checkEquals(dim(tbl.enhancers.myeloid), c(9, 16))

} # test_getEnhancersForTissues
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   showGenomicRegion(igv, "chr3:128,602,838-128,611,664")
   showGenomicRegion(igv, "chr3:128071944-128623958")
   setTargetGene(tp, "GATA2")
   all.tissues <- getEnhancerTissues(tp)
   myeloid.tissues <- grep("myeloid", all.tissues, ignore.case=TRUE, value=TRUE)

   tbl.enhancers.myeloid <- getEnhancersForTissues(myeloid.tissues, display=TRUE, trackName="myeloid enhancers")
   tbl.enhancers.all <- getEnhancersForTissues("all", display=TRUE, trackName="all enhancers")

   display_ATACseq_in_enhancers(tbl.enhancers.myeloid, trackName="myeloid enhancers")

   display_ATACseq_in_enhancers()
   displayMotifsForTF()
   hopeRestored()
   showGenomicRegion(igv, "chr3:128,479,295-128,503,294")

   #display_ATACseq_in_enhancers(tbl.enhancers.all)
   #makeAtacBasedModels()

} # run
#------------------------------------------------------------------------------------------------------------------------
