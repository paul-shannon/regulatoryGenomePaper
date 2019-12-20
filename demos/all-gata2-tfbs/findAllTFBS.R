source("../common.R")
tbl.enhancers.all <- getEnhancersForTissues("all", display=TRUE, trackName="GeneHancer")
conservationTrack()
tbl.atacAll <- display_ATACseq_in_enhancers()

tbl.atacRich <- tbl.atacAll[c(1,3,4),]

motifs.hsapiens <- MotifDb::query(MotifDb, c("hsapiens"), c("jaspar2018", "hocomoco"))
length(motifs.hsapiens)  # 1177
meme.file <- "human.meme"
export(motifs.hsapiens, con=meme.file, format="meme")

fimo <- 1e-4
phast7 <- 0.9

tbl.tfbs <- getTFBS.fimo(tv, tbl.atacRich, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
dim(tbl.tfbs)
tbl.freq <- as.data.frame(sort(table(tbl.tfbs$tf)), stringsAsFactors=FALSE)
colnames(tbl.freq) <- c("tf", "count")
deleters <- setdiff(tbl.freq$tf, rownames(mtx))
length(deleters)
if(length(deleters) > 0){
   delete.indices <- match(deleters, tbl.freq$tf)
   tbl.freq <- tbl.freq[-delete.indices,]
   }
correlations <- unlist(lapply(tbl.freq$tf, function(tf) cor(mtx["GATA2",], mtx[tf,])))
tbl.freq$cor <- correlations
tbl.freqByCor <- tbl.freq[order(abs(tbl.freq$cor), decreasing=TRUE),]
tail(tbl.freq, n=10)
mtx <- getExpressionMatrix(tp, "brandLabDifferentiationTimeCourse-27171x28")
gata2 <- as.numeric(mtx["GATA2",])
gene2.name <- "HOXC13"
gene2 <- as.numeric(mtx[gene2.name,])
cor(gata2, gene2)


plot(gata2, pch=16, type="b", xaxt="n", xlab="day", ylab="mRNA", main=sprintf("%s vs. GATA2", gene2.name), ylim=c(0,12))
x.axis.labels <- c("0.1", "0.2", "2.1", "2.2", "4.1", "4.2", "6.1", "6.2", "7.1", "7.2", "8.1", "8.2", "8.1", "8.2",
                   "10.1", "10.2", "10.1", "10.2", "11.1", "11.2", "11.1", "11.2", "12.1", "12.2", "14.1", "14.2",
                   "16.1", "16.2")

x.locs  <- c(0, 2, 4, 6, 8, 10, 14, 18, 22, 24, 26)
days    <- c(0, 2, 4, 6, 7,  8, 10, 11, 12, 14, 16)
axis(side=1, at=x.locs, labels=days)
lines(gene2, type="b", pch=16, col="blue")

legend(23, 12, c("GATA2", gene2.name), c("black", "blue"))
text(1.5, 7, "MPP")
text(5, 7, "MEP")
text(12, 7, "CPU-E")
text(20, 7, "ProEB")
text(25, 7, "BasoEB")
