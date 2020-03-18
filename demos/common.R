library(TrenaProjectHG38.generic)
library(org.Hs.eg.db)
library(GO.db)
if(interactive()) require(igvR)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38
library(TrenaValidator)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProjectErythropoiesis)
library(trenaSGM)
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

if(interactive() && !exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   setBrowserWindowTitle(igv, "GATA2")
  showGenomicRegion(igv, "GATA2")
   #later(function(){
   #         track <- DataFrameQuantitativeTrack("gh", tbl.enhancers[, c("chrom", "start", "end", "combinedScore")],
   #                                             autoscale=FALSE, color="blue", min=0, max=50)
   #         displayTrack(igv, track)
   #         }, 4)
   } # igv
#------------------------------------------------------------------------------------------------------------------------
# if(!exists("tv")) {
#   benchmark.full <- "~/github/trena/misc/saez-benchmark-paper/GarciaAlonso_Supplemental_Tables/database.csv"
#   tbl.bm <-read.table(benchmark.full, sep=",", as.is=TRUE, header=TRUE, nrow=-1)
#   message(sprintf("--- creating instance of TrenaValidator"))
#   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
#   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
#   #mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
#   #setMatrix(tv, mtx)
#   tv <- TrenaValidator(TF="TWIST1", "GATA2", tbl.benchmark);
#   tp.hg38 <- TrenaProjectHG38.generic()
#   }
#------------------------------------------------------------------------------------------------------------------------
conservationTrack <- function()
{
   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)
   tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.blocks)), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

} # conservationTrack
#------------------------------------------------------------------------------------------------------------------------
chipTrack <- function(tf, tissueRestriction=NA_character_, peaksAlso=FALSE)
{
   chrom.loc <- getGenomicRegion(igv)
   tbl.chip <- with(chrom.loc, getChipSeq(tp.hg38, chrom, start, end, tf))

   if(nrow(tbl.chip) == 0){
      printf("--- no ChIP found for %s", tf)
      return()
      }

   if(!is.na(tissueRestriction))
      tbl.chip <- tbl.chip[grep(tissueRestriction, tbl.chip$name),]

   if(nrow(tbl.chip) == 0){
      printf("--- no ChIP found for %s", tf)
      return()
      }

   printf("tf %s,  %d peaks", tf,  nrow(tbl.chip))

   if(nrow(tbl.chip) == 0){
      #printf("no ChIP for TF %s in %d bases", tf, chrom.loc$end - chrom.loc$start)
      return(data.frame())
      }
   tbl.track <- tbl.chip[, c("chrom", "start", "endpos", "name")]
   trackName <- sprintf("Ch-%s", tf)
   track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random", displayMode="squished", trackHeight=25)
   displayTrack(igv, track)

   if(peaksAlso){
      tbl.track <- tbl.chip[, c("chrom", "peakStart", "peakEnd", "name")]
      trackName <- sprintf("peaks-%s", tf)
      track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random", displayMode="squished", trackHeight=25)
      displayTrack(igv, track)
      } # peaksAlso

} # chipTrack
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
      track <- DataFrameAnnotationTrack(sample.x, tbl.sample, color=color, trackHeight=23)
      displayTrack(igv, track)
      #write.table(tbl.sample, file=sprintf("tbl.%s.bed", sample.x), quote=FALSE, row.names=FALSE)
      } # for sample.x

   tbl.regions.condensed <- as.data.frame(union(GRanges(tbl.atac.trimmed[, c("chrom", "start", "end")]),
                                                GRanges(tbl.atac.trimmed[, c("chrom", "start", "end")])))[, c("seqnames", "start", "end")]

   colnames(tbl.regions.condensed) <- c("chrom", "start", "end")
   tbl.regions.condensed$chrom <- as.character(tbl.regions.condensed$chrom)
   track <- DataFrameAnnotationTrack("atac combined", tbl.regions.condensed, color="black", trackHeight=23)
   displayTrack(igv, track)

   state$atac.regions.by.sample <- atac.regions.by.sample
   return(tbl.regions.condensed)

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
