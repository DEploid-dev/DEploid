rm(list= ls())

fun.llk<-function(cov.ref, cov.alt, f.samp, err=0.01, fac=100) {
    f.samp<-f.samp+err*(1-2*f.samp);
    llk<-lbeta(cov.alt+f.samp*fac, cov.ref+(1-f.samp)*fac)-lbeta(f.samp*fac,(1-f.samp)*fac);
    if (sum(is.nan(llk))>1){
      print("f.samp = ")
    }
#    llk<-lgamma(fac*f.samp+cov.alt)+lgamma(fac*(1-f.samp)+cov.ref)-lgamma(fac*f.samp)-lgamma(fac*(1-f.samp));
    return(llk);
}

fun.dic.by.llk.var = function ( tmpllk ){
     return (  mean(-2*tmpllk) + var(-2*tmpllk)/2 )# D_bar + 1/2 var (D_theta), where D_theta = -2*tmpllk, and D_bar = mean(D_theta)
}

fun.dic.by.theta = function ( tmpllk, thetallk ){
    DIC.WSAF.bar = -2 * sum(thetallk)
    return (  mean(-2*tmpllk) + (mean(-2*tmpllk) - DIC.WSAF.bar) ) # D_bar + pD, where pD = D_bar - D_theta, and D_bar = mean(D_theta)
}

for (i in 1:15){
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
     dic.by.var = fun.dic.by.llk.var (llk)
     dic.by.theta = fun.dic.by.theta ( llk, fun.llk(ref, alt, expWSAF))
     plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="LLK", main=paste("LLK sd:", round(llk_sd, digits = 4),
                                                                                      ", dic.by.var: ",round(dic.by.var),
                                                                                      ", dic.by.theta: ",round(dic.by.theta)));
     updateSingleAt = which(llkEvent == 1)
     updateBothAt = which(llkEvent == 2)
     updatePropAt = which(llkEvent == 0)
     index = c(1:length(llk))
     points(index[updateSingleAt], llk[updateSingleAt], cex = 0.6, col="red")
     points(index[updateBothAt], llk[updateBothAt], cex = 0.6, col="blue")
     points(index[updatePropAt], llk[updatePropAt], cex = 0.6, col="green")
     dev.off()
}
