tmp = read.table("tmp.prop", header=F)
rainbowColorBin = 12
barplot(t(tmp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="Iteration", ylab="Component Freq", main="Components")
