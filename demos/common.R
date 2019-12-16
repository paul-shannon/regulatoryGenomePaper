library(TrenaProjectHG38.generic)
library(org.Hs.eg.db)
library(GO.db)
library(igvR)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38
library(TrenaValidator)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   benchmark.full <- "~/github/trena/misc/saez-benchmark-paper/GarciaAlonso_Supplemental_Tables/database.csv"
   tbl.bm <-read.table(benchmark.full, sep=",", as.is=TRUE, header=TRUE, nrow=-1)
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   #mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   #setMatrix(tv, mtx)
   tv <- TrenaValidator(TF="TWIST1", "GATA2", tbl.benchmark);
   tp.hg38 <- TrenaProjectHG38.generic()
   }
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
