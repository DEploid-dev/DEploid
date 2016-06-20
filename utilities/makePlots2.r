rm(list= ls())

plot.prob <-function (tmpProp, title){
    rainbowColorBin = 16
    barplot(t(tmpProp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="SNP index", ylab="Posterior probabilities", main=title)
}

plot.postProb.ofCase <- function ( case ){
    png(paste(prefix, ".", case, ".png", sep = ""), width = 3500, height = 2000)
    obj = read.table( paste(prefix, ".", case, sep = ""), header=T)
    chromName = levels(obj$CHROM)
    par(mfrow = c(length(chromName)/2,2))
    for ( chromI in chromName ){
        plot.prob ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])], "")
    }
    dev.off()
}

args=(commandArgs(TRUE))

prefix = args[1]

plot.postProb.ofCase("single1")
plot.postProb.ofCase("single2")
plot.postProb.ofCase("pair")

