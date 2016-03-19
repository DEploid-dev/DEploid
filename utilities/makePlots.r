rm(list= ls())
#R --slave "--args prefix refPath altPath plafPath excludePath" < makePlots.r > plotLLK.suffix.rout
# R --slave "--args PG0389.C.excluded.noPanel.seed1 labStrains/PG0389.C_ref.txt labStrains/PG0389.C_alt.txt labStrains/labStrains_samples_PLAF.txt labStrains/PG0389.C.exclude.csv " < utilities/makePlots.r
# R --slave "--args PG0389.C.noPanel.seed1 labStrains/PG0389.C_ref.txt labStrains/PG0389.C_alt.txt labStrains/labStrains_samples_PLAF.txt " < utilities/makePlots.r
# R --slave "--args PG0389.C.excluded.panel.seed1 labStrains/PG0389.C_ref.txt labStrains/PG0389.C_alt.txt labStrains/labStrains_samples_PLAF.txt labStrains/PG0389.C.exclude.csv "< makePlots.r

fun.llk <- function(cov.ref, cov.alt, f.samp, err=0.01, fac=100) {
    f.samp<-f.samp+err*(1-2*f.samp);
    llk<-lbeta(cov.alt+f.samp*fac, cov.ref+(1-f.samp)*fac)-lbeta(f.samp*fac,(1-f.samp)*fac);
    #  llk<-lgamma(fac*f.samp+cov.alt)+lgamma(fac*(1-f.samp)+cov.ref)-lgamma(fac*f.samp)-lgamma(fac*(1-f.samp));
    if (sum(is.nan(llk))>1){
      print("f.samp = ")
    }
    return(llk);
}

fun.dic.by.llk.var <- function ( tmpllk ){
     return (  mean(-2*tmpllk) + var(-2*tmpllk)/2 )# D_bar + 1/2 var (D_theta), where D_theta = -2*tmpllk, and D_bar = mean(D_theta)
}

fun.dic.by.theta <- function ( tmpllk, thetallk ){
    DIC.WSAF.bar = -2 * sum(thetallk)
    return (  mean(-2*tmpllk) + (mean(-2*tmpllk) - DIC.WSAF.bar) ) # D_bar + pD, where pD = D_bar - D_theta, and D_bar = mean(D_theta)
}

plot.prop <-function (propMat, title){
    rainbowColorBin = 5
    barplot(t(tmpProp), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="Iteration", ylab="Component Freq", main=title)
}

plot.altVsRef <- function ( ref, alt, title, exclude.ref = c(), exclude.alt = c() ){
    tmp.range = 1.1*mean(max(alt), max(ref))
    plot ( ref, alt, xlim=c(0,tmp.range), ylim=c(0,tmp.range), cex = 0.5, xlab = "REF", ylab = "ALT", main = title)
    points (exclude.ref, exclude.alt, col = "red")
    abline(v =50, untf = FALSE, lty = 2)
    abline(h =50, untf = FALSE, lty = 2)

    abline(h =150, untf = FALSE, lty = 2)
    abline(v =150, untf = FALSE, lty = 2)
}

plot.wsaf.hist <- function ( obsWSAF, title ="Histogram 0<WSAF<1" ){
    hist(obsWSAF, main=title, breaks = seq(0, 1, by =0.1), xlab = "WSAF")
}

plot.llk <- function (llkTable, ref, alt, expWSAF, logFileName ){
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

    cat ( "dic.by.var: ", dic.by.var, "\n", file = logFileName, append = T)
    cat ( "dic.by.theta: ", dic.by.theta, "\n", file = logFileName, append = T)

}

plot.plaf.vs.wsaf <- function ( plaf, obsWSAF, expWSAF, title = "PLAF vs WSAF" ){
    plot ( plaf, obsWSAF, col = "red", main = title )
    points ( plaf, expWSAF, col = "blue")
}


plot.wsaf <- function (obsWSAF, expWSAF, logFileName){
    currentWSAFcov  = cov(obsWSAF, expWSAF)
    currentWSAFcorr = cor(obsWSAF, expWSAF)

    plot(obsWSAF, expWSAF, pch=19, col="blue", xlab="Observed WSAF (ALT/(ALT+REF))", ylab="Expected WSAF (h%*%p)",
         main=paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"),
         xlim = c(-0.05, 1.05), cex = 0.5, ylim = c(-0.05, 1.05));
    abline(0,1,lty="dotted");

    cat ( "corr: ",  currentWSAFcorr, "\n", file = logFileName)
}

args=(commandArgs(TRUE))
print(args)

prefix = args[1]
refFile = args[2]
altFile = args[3]
ref.full = read.table(refFile, header = T)
alt.full = read.table(altFile, header = T)

plafFile = args[4]

logFileName  = paste(prefix, "dic.log", sep = "")
propFileName = paste(prefix, ".prop",   sep = "")
hapFileName  = paste(prefix, ".hap",    sep = "")
llkFileName  = paste(prefix, ".llk",    sep = "")


# Extract and compute stuff
if ( length(args) == 5 ){
    excludeFile = args[5]
    exclude = read.csv(excludeFile, header = T)
    excludeLogic = ( paste(ref.full$CHROM, ref.full$POS) %in% paste(exclude$CHROM, exclude$POS) )
    excludeindex = which(excludeLogic)
    includeindex = which(!excludeLogic)

    ref = ref.full[includeindex,3]
    alt = alt.full[includeindex,3]
    exclude.ref = ref.full[excludeindex,3]
    exclude.alt = alt.full[excludeindex,3]

    plaf = read.table(plafFile, header = T)[includeindex,3]
} else {
    ref = ref.full[,3]
    alt = alt.full[,3]
    plaf = read.table(plafFile, header = T)[,3]
}

tmpProp = read.table(propFileName, header=F)
prop = as.numeric(tmpProp[dim(tmpProp)[1],])

hap = as.matrix(read.table(hapFileName, header=T)[,-c(1,2)] )
obsWSAF = alt/(alt+ref+0.00000001)
expWSAF = hap %*%prop
llkTable = read.table( llkFileName, header=F)

# plotting
unlink ( paste( prefix, "*.png", sep= ""))
png(paste( prefix, ".png", sep= ""), width = 1500, height = 1000)
par (mfrow = c(2,3))

if ( length(args) == 5 ){
    plot.altVsRef ( ref, alt, "Alt vs Ref", exclude.ref, exclude.alt )
} else {
    plot.altVsRef ( ref, alt, "Alt vs Ref" )
}
plot.wsaf.hist ( obsWSAF )
plot.prop( tmpProp, "Components" )
plot.wsaf ( obsWSAF, expWSAF, logFileName )
plot.plaf.vs.wsaf (plaf, obsWSAF, expWSAF)
plot.llk( llkTable, ref, alt, expWSAF, logFileName)

dev.off()

