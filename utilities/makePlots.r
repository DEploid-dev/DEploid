tmpProp = read.table("tmp.prop", header=F)
rainbowColorBin = 12
barplot(t(tmpProp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="Iteration", ylab="Component Freq", main="Components")


ref = read.table("tests/PG0393_ref.txt", header=T)$V3
alt = read.table("tests/PG0393_alt.txt", header=T)$V3

prop = as.numeric(tmpProp[dim(tmpProp)[1],])
hap = as.matrix(read.table("tmp.hap", header=F))

obsWSAF = alt/(alt+ref+0.00000001)
expWSAF = hap %*%prop


currentWSAFcov  = cov(obsWSAF, expWSAF)
currentWSAFcorr = cor(obsWSAF, expWSAF)

plot(obsWSAF, expWSAF, pch=19, col="blue", xlab="Observed WSAF (ALT/(ALT+REF))", ylab="Expected WSAF (h%*%p)",
     main=paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"),
     xlim = c(-0.05, 1.05), cex = 0.5, ylim = c(-0.05, 1.05));
abline(0,1,lty="dotted");


llkTable = read.table("tmp.llk", header=F)
llk = llkTable$V2
llkEvent = llkTable$V1
llk_sd = sd(llk)
llk_range = range(llk)
plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="LLK", main=paste("LLK sd:", round(llk_sd, digits = 4),", range:",round(llk_range[1]),",",round(llk_range[2])));
updateSingleAt = which(llkEvent == 1)
updateBothAt = which(llkEvent == 2)
updatePropAt = which(llkEvent == 0)
index = c(1:length(llk))
points(index[updateSingleAt], llk[updateSingleAt], cex = 0.6, col="red")
points(index[updateBothAt], llk[updateBothAt], cex = 0.6, col="blue")
points(index[updatePropAt], llk[updatePropAt], cex = 0.6, col="green")
