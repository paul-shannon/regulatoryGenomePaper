source("../common.R")
#------------------------------------------------------------------------------------------------------------------------
buildModel <- function(tfs)
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
gata2.brand <- function()
{
  tbl.enhancers.all <- getEnhancersForTissues("all", display=TRUE, trackName="all enhancers")
  roi <- with(tbl.enhancers.all, sprintf("%s:%d-%d", chrom[1], min(start)-5000, max(end)+5000))
  showGenomicRegion(igv, roi)

  motifs.jaspar2018 <- MotifDb::query(MotifDb, c("hsapiens", "jaspar2018"))
  length(motifs.jaspar2018)
  meme.file.jaspar2018 <- "human.jaspar2018.meme"

  export(motifs.jaspar2018, con=meme.file.jaspar2018, format="meme")

  tbl.fimo <- fimoBatch(tbl.enhancers.all, matchThreshold=1e-4, "hg38", meme.file.jaspar2018)
  tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.fimo)), stringsAsFactors=FALSE)
  tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
  tbl.filtered <-  subset(tbl.cons7, default > 0.9)
  tfs <- unique(tbl.filtered$tf)
  "TBX15" %in% tfs

  mtx.name <- "brandLabDifferentiationTimeCourse-27171x28"
  stopifnot(mtx.name %in% getExpressionMatrixNames(tp))
  mtx <- getExpressionMatrix(tp, mtx.name)
  dim(mtx)

  solverNames <- c("lasso", "ridge", "pearson", "spearman", "randomForest", "xgboost")
  solver <- EnsembleSolver(mtx, targetGene="GATA2", candidateRegulators=tfs, solverNames=solverNames)
  tbl.model <- run(solver)
  tbl.model <- tbl.model[order(abs(tbl.model$pearsonCoeff), decreasing=TRUE),]


  conservationTrack()

} # gata2.brand
#------------------------------------------------------------------------------------------------------------------------
gata2.buenrosto <- function()
{
  tbl.enhancers.all <- getEnhancersForTissues("all", display=TRUE, trackName="all enhancers")
  roi <- with(tbl.enhancers.all, sprintf("%s:%d-%d", chrom[1], min(start)-5000, max(end)+5000))
  showGenomicRegion(igv, roi)

  motifs.jaspar2018 <- MotifDb::query(MotifDb, c("hsapiens", "jaspar2018"))
  length(motifs.jaspar2018)
  meme.file.jaspar2018 <- "human.jaspar2018.meme"

  export(motifs.jaspar2018, con=meme.file.jaspar2018, format="meme")

  tbl.fimo <- fimoBatch(tbl.enhancers.all, matchThreshold=1e-4, "hg38", meme.file.jaspar2018)
  tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.fimo)), stringsAsFactors=FALSE)
  tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
  tbl.cons7 <-  subset(tbl.cons7, default > 0.0)
  tfs <- sort(unique(tbl.cons7$tf))
  length(tfs)

  tbl.atac <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74912-atacSeq.RData"))
  colnames(tbl.atac)[1:3] <- c("chrom", "start", "end")
  tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.cons7[, c("seqnames", "start", "end")]),
                                       GRanges(tbl.atac[,c("chrom", "start", "end")])))
  dim(tbl.ov)
  tbl.cons7.atac <- tbl.cons7[unique(tbl.ov$queryHits),]
  tfs <- sort(unique(tbl.cons7.atac$tf))
  length(tfs)

  mtx <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/buenrostoErythropoiesisAsinh.RData"))

  solverNames <- c("lasso", "ridge", "pearson", "spearman", "randomForest", "xgboost")
  solver <- EnsembleSolver(mtx, targetGene="GATA2", candidateRegulators=tfs, solverNames=solverNames)
  tbl.model <- run(solver)
  tbl.model <- tbl.model[order(abs(tbl.model$rfScore), decreasing=TRUE),]
  rownames(tbl.model) <- NULL
  dim(tbl.model)

  # count binding sites
  tfbs.list <- as.list(table(subset(tbl.cons7, tf %in% tbl.model$gene)$tf))
  tbl.model$bindingSites <- as.numeric(tfbs.list[tbl.model$gene])
  tbl.model <- head(tbl.model, n=10)
  roundNumericColumnsInDataframe(tbl.model, 2)
  conservationTrack()

  # MEIS1 regulates early erythroid and megakaryocytic cell fate.
  # https://www.ncbi.nlm.nih.gov/pubmed/25107888
  # MEIS1 increased expression of genes associated with megakaryopoiesis and erythropoiesis such as KLF1, HBD, HBG, SLC40A1, THBS1, GPIb, VWA5A and GATA2

} # gata2.buenrosto
#------------------------------------------------------------------------------------------------------------------------
tal1.buenrosto <- function()
{
  tbl.enhancers.all <- getEnhancersForTissues("all", display=TRUE, trackName="all enhancers")
  roi <- with(tbl.enhancers.all, sprintf("%s:%d-%d", chrom[1], min(start)-5000, max(end)+5000))
  showGenomicRegion(igv, roi)

  motifs.jaspar2018 <- MotifDb::query(MotifDb, c("hsapiens", "jaspar2018"))
  length(motifs.jaspar2018)
  meme.file.jaspar2018 <- "human.jaspar2018.meme"

  export(motifs.jaspar2018, con=meme.file.jaspar2018, format="meme")

  tbl.fimo <- fimoBatch(tbl.enhancers.all, matchThreshold=1e-4, "hg38", meme.file.jaspar2018)
  tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.fimo)), stringsAsFactors=FALSE)
  tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
  tbl.cons7 <-  subset(tbl.cons7, default > 0.5)
  tfs <- sort(unique(tbl.cons7$tf))
  length(tfs)
  "TBX15" %in% tfs

  mtx <- get(load("~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/buenrostoErythropoiesisAsinh.RData"))

  solverNames <- c("lasso", "ridge", "pearson", "spearman", "randomForest", "xgboost")
  solver <- EnsembleSolver(mtx, targetGene="FLI1", candidateRegulators=tfs, solverNames=solverNames)
  tbl.model <- run(solver)
  tbl.model <- tbl.model[order(abs(tbl.model$rfScore), decreasing=TRUE),]
  rownames(tbl.model) <- NULL
  dim(tbl.model)

  # count binding sites
  tfbs.list <- as.list(table(subset(tbl.cons7, tf %in% tbl.model$gene)$tf))
  tbl.model$bindingSites <- as.numeric(tfbs.list[tbl.model$gene])
  tbl.model
  conservationTrack()

  # MEIS1 regulates early erythroid and megakaryocytic cell fate.
  # https://www.ncbi.nlm.nih.gov/pubmed/25107888
  # MEIS1 increased expression of genes associated with megakaryopoiesis and erythropoiesis such as KLF1, HBD, HBG, SLC40A1, THBS1, GPIb, VWA5A and GATA2

} # tal1.buenrosto
#------------------------------------------------------------------------------------------------------------------------
