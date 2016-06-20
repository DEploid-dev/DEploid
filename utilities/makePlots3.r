rm(list= ls())

plot.switchMissCopy.ofCase <- function ( case ) {
    obj = read.table(paste(prefix, ".",case, sep=""), header=F)
    chromName = levels(obj$V1)
    png(paste(prefix, ".", case, ".png", sep = ""), width = 3500, height = 2000)

    par(mfrow = c(length(chromName)/2,2))
    for ( chromI in chromName ){
        tmpTable = table(obj$V2[which( chromI == obj$V1)])/4800
        index = as.numeric(names(tmpTable))
        myfreq = rep(0, length(which(CHROM == chromI)))
        tmpPOS = POS[which(CHROM == chromI)]
        myfreq[which(tmpPOS %in% index)] = tmpTable
        barplot(myfreq, ylim=c(0,1))
    }
    dev.off()
}

args=(commandArgs(TRUE))

prefix = args[1]
#prefix="PG0391"

# load Index
single1 = read.table("tests/testData/PG0410.C_ref.txt", header=T)
#single1 = read.table( paste(prefix, ".single1", sep = ""), header=T)
CHROM = single1$CHROM
POS = single1$POS

plot.switchMissCopy.ofCase ( "oneSwitchOne" )
plot.switchMissCopy.ofCase ( "oneMissCopyOne" )
plot.switchMissCopy.ofCase ( "twoSwitchOne" )
plot.switchMissCopy.ofCase ( "twoMissCopyOne" )
plot.switchMissCopy.ofCase ( "twoSwitchTwo" )
plot.switchMissCopy.ofCase ( "twoMissCopyTwo" )


