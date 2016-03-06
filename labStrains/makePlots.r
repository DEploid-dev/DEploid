for (i in 1:30){
     dataDir = "../"
     currentDir = system("echo ${PWD##*/}", intern = T)

     prefix = paste(currentDir,"_seed", i, sep="")
     #prefix = "PG0394_ind_seed1"

     unlink ( paste( prefix, "*.png", sep= ""))
     png(paste( prefix, ".png", sep= ""), width = 1000, height = 1000)

     par (mfrow = c(2,2))

     tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
     rainbowColorBin = 12
     barplot(t(tmpProp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="Iteration", ylab="Component Freq", main="Components")

     ref = read.table(paste(dataDir, currentDir,"_ref.txt",sep = ""), header = T)[,3]
     alt = read.table(paste(dataDir, currentDir,"_alt.txt",sep = ""), header = T)[,3]
     prop = as.numeric(tmpProp[dim(tmpProp)[1],])
     hap = as.matrix(read.table(paste(prefix,".hap",sep=""), header=F))

     obsWSAF = alt/(alt+ref+0.00000001)
     expWSAF = hap %*%prop

     currentWSAFcov  = cov(obsWSAF, expWSAF)
     currentWSAFcorr = cor(obsWSAF, expWSAF)

     plot(obsWSAF, expWSAF, pch=19, col="blue", xlab="Observed WSAF (ALT/(ALT+REF))", ylab="Expected WSAF (h%*%p)",
          main=paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"),
          xlim = c(-0.05, 1.05), cex = 0.5, ylim = c(-0.05, 1.05));
     abline(0,1,lty="dotted");

     plaf = read.table(paste(dataDir, "labStrains_samples_PLAF.txt",sep = ""), header = T)[,3]
     plot ( plaf, obsWSAF, col = "red")
     points ( plaf, expWSAF, col = "blue")
#     dev.off()

#     png(paste( prefix, ".llk.png", sep= ""))
     llkTable = read.table( paste( prefix, ".llk", sep=""), header=F)
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
     dev.off()
}
