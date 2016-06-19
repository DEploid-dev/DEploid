rm(list= ls()

plot.prob <-function (tmpProp, title){
    rainbowColorBin = 16
#    rainbowColorBin = 5
    barplot(t(tmpProp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="SNP index", ylab="Posterior probabilities", main=title)
}

prefix = args[1]

png(paste(prefix, ".single1.png", sep = ""), width = 3500, height = 2000)
single1 = read.table( paste(prefix, ".single1", sep = ""), header=T)
chromName = levels(single1$CHROM)
par(mfrow = c(length(chromName)/2,2))
for ( chromI in chromName ){
    plot.prob ( single1[which( chromI == single1$CHROM),c(3:dim(single1)[2])], "")
}
dev.off()


png(paste(prefix, ".single2.png", sep = ""), width = 3500, height = 2000)
single2 = read.table( paste(prefix, ".single2", sep = ""), header=T)
chromName = levels(single2$CHROM)
par(mfrow = c(length(chromName)/2,2))
for ( chromI in chromName ){
    plot.prob ( single2[which( chromI == single2$CHROM),c(3:dim(single2)[2])], "")
}


png(paste(prefix, ".pair.png", sep = ""), width = 3500, height = 2000)
pair = read.table( paste(prefix, ".pair", sep = ""), header=T)
chromName = levels(pair$CHROM)
par(mfrow = c(length(chromName)/2,2))
for ( chromI in chromName ){
    plot.prob ( pair[which( chromI == pair$CHROM),c(3:dim(pair)[2])], "")
}
dev.off()
