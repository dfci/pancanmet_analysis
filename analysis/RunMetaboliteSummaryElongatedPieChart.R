library(RColorBrewer)
library(treemap)

x <- read.table(file.path("scratch", "scratch.txt"))
tmp <- as.matrix(x$V1)

#image(scale(tmp, center=FALSE), col=brewer.pal(10, "Spectral"), axes=FALSE)
#image(as.matrix(1:length(tmp)), col=brewer.pal(10, "Spectral"), axes=FALSE)
#x1 <- cut(x$V1, seq(min(x$V1), max(x$V1), length.out=10), right=FALSE, labels=c(1:9))

tx <- as.data.frame(tmp)

colnames(tx) <- "freq"

pdf(file="tree.pdf", width=21, height=1)
treemap(tx, index=c("freq"), vSize = "freq", title="")
dev.off()
