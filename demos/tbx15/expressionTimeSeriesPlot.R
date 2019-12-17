library(TrenaProjectErythropoiesis)
tpe <- TrenaProjectErythropoiesis()
getExpressionMatrixNames(tpe)
mtx <- getExpressionMatrix(tpe, "brandLabDifferentiationTimeCourse-27171x28")
dim(mtx)

gata2 <- as.numeric(mtx["GATA2",])
tbx15 <- as.numeric(mtx["TBX15",])
cor(tbx15, gata2)

plot(gata2, pch=16, type="b", xaxt="n", xlab="day", ylab="mRNA", main="TBX15 vs. GATA2", ylim=c(0,12))
x.axis.labels <- c("0.1", "0.2", "2.1", "2.2", "4.1", "4.2", "6.1", "6.2", "7.1", "7.2", "8.1", "8.2", "8.1", "8.2",
                   "10.1", "10.2", "10.1", "10.2", "11.1", "11.2", "11.1", "11.2", "12.1", "12.2", "14.1", "14.2",
                   "16.1", "16.2")

x.locs  <- c(0, 2, 4, 6, 8, 10, 14, 18, 22, 24, 26)
days    <- c(0, 2, 4, 6, 7,  8, 10, 11, 12, 14, 16)
axis(side=1, at=x.locs, labels=days)
lines(tbx15, type="b", pch=16, col="blue")
legend(23, 12, c("GATA2", "TBX15"), c("black", "blue"))
text(1.5, 7, "MPP")
text(5, 7, "MEP")
text(12, 7, "CPU-E")
text(20, 7, "ProEB")
text(25, 7, "BasoEB")

