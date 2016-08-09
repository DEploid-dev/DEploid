fun.parse <- function( args ){
    fun.local.checkAndIncreaseArgI <- function ( ){
        arg_i = arg_i+1
    }

    outPrefix = "dataExplore"
    vcfFileName = ""
    plafFileName = ""
    excludeFileName = ""
    dEploidPrefix = ""
    excludeBool = FALSE

    arg_i = 1
    while ( arg_i < length(args) ){
        argv = args[arg_i]
        if ( argv == "-vcf" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            vcfFileName = args[arg_i]
        } else if ( argv == "-plaf" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            plafFileName = args[arg_i]
        } else if ( argv == "-exclude" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            excludeFileName = args[arg_i]
            excludeBool = TRUE
        } else if ( argv == "-o" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            outPrefix = args[arg_i]
        } else if ( argv == "-dEprefix" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            dEploidPrefix = args[arg_i]
        } else {
            cat ("Unknow flag: ", argv, "\n")
        }

        arg_i = arg_i + 1
    }

    if ( vcfFileName == "" ){
        stop ("Vcf File name not specified!")
    }

    if ( plafFileName == "" ){
        stop ("Plaf File name not specified!")
    }

#    cat ("vcfFileName: ", vcfFileName, "\n")
#    cat ("plafFileName: ", plafFileName, "\n")

    return ( list ( vcfFileName = vcfFileName,
                    plafFileName = plafFileName,
                    outPrefix = outPrefix,
                    dEploidPrefix = dEploidPrefix,
                    excludeFileName = excludeFileName) )
}


fun.dEploidPrefix <- function ( prefix ){
    if ( prefix == "" ){
        stop ("dEprefix ungiven!!!")
    }

    return ( list ( propFileName = paste(prefix, ".prop",   sep = ""),
                    hapFileName  = paste(prefix, ".hap",    sep = ""),
                    llkFileName  = paste(prefix, ".llk",    sep = "") ) )
}


fun.extract.vcf <- function ( vcfName, ADFieldIndex = 2 ){
    # Assume that AD is the second field
    skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), T))
    vcf  = read.table( vcfName, skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

    sampleName = names(vcf)[10]

    tmp = vcf[[sampleName]]
    field = strsplit(as.character(tmp),":")

    tmpCovStr = unlist(lapply(field, `[[`, ADFieldIndex))
    tmpCov = strsplit(as.character(tmpCovStr),",")

    refCount = as.numeric(unlist(lapply(tmpCov, `[[`, 1)))
    altCount = as.numeric(unlist(lapply(tmpCov, `[[`, 2)))

    return ( data.frame( CHROM = vcf[,1],
                         POS = vcf[,2],
                         refCount = refCount,
                         altCount = altCount )
           )
}


fun.extract.plaf <- function ( plafName ){
    return ( read.table(plafName, header=T) )
}


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


plot.altVsRef <- function ( ref, alt, title = "Alt vs Ref", exclude.ref = c(), exclude.alt = c() ){
    tmp.range = 1.1*mean(max(alt), max(ref))
    plot ( ref, alt, xlim=c(0,tmp.range), ylim=c(0,tmp.range), cex = 0.5, xlab = "REF", ylab = "ALT", main = title)
    points (exclude.ref, exclude.alt, col = "red")
    abline(v =50, untf = FALSE, lty = 2)
    abline(h =50, untf = FALSE, lty = 2)

    abline(h =150, untf = FALSE, lty = 2)
    abline(v =150, untf = FALSE, lty = 2)
}


plot.wsaf.hist <- function ( obsWSAF, indexed = TRUE, title ="Histogram 0<WSAF<1" ){
    tmpWSAF_index = 1:length(obsWSAF)
    if ( indexed ){
        tmpWSAF_index = which(((obsWSAF<1) * (obsWSAF>0) ) == 1)
    }
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

#    cat ( "dic.by.var: ", dic.by.var, "\n", file = logFileName, append = T)
#    cat ( "dic.by.theta: ", dic.by.theta, "\n", file = logFileName, append = T)

}


plot.plaf.vs.wsaf <- function ( plaf, obsWSAF, expWSAF = c(), title = "PLAF vs WSAF" ){
    plot ( plaf, obsWSAF, cex = 0.5, xlim = c(0,1), ylim = c(0,1), col = "red", main = title, xlab = "PLAF", ylab = "WSAF" )
    if ( length(expWSAF) > 0 ){
        points ( plaf, expWSAF, cex = 0.5, col = "blue")
    }
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


fun.calc.obsWSAF <- function (alt, ref) {
    return ( alt / (ref + alt + 0.00000001) )
}


fun.dataExplore <- function (vcfInfo, plafInfo, prefix = "") {
    PLAF = plafInfo$PLAF
    ref = vcfInfo$refCount
    alt = vcfInfo$altCount

    png ( paste ( prefix, "altVsRefAndWSAFvsPLAF.png", sep = "" ), width = 1800, height = 600)
    par( mfrow = c(1,3) )

    plot.altVsRef ( ref, alt )

    WSAF = fun.calc.obsWSAF ( alt, ref )

    plot.wsaf.hist ( WSAF )

    plot.plaf.vs.wsaf ( PLAF, WSAF )

    dev.off()
}


fun.interpretDEploid <- function (vcfInfo, plafInfo, prefix = "") {
    PLAF = plafInfo$PLAF
    ref = vcfInfo$refCount
    alt = vcfInfo$altCount

    png ( paste ( prefix, ".interpretDEploidFigure.1.png", sep = "" ),  width = 1500, height = 1000)
    par( mfrow = c(2,3) )

    plot.altVsRef ( ref, alt )

    obsWSAF = fun.calc.obsWSAF ( alt, ref )
    plot.wsaf.hist ( obsWSAF )
    plot.plaf.vs.wsaf ( PLAF, obsWSAF )

plot.prop( tmpProp, "Components" )
plot.wsaf ( obsWSAF, expWSAF, logFileName )
plot.plaf.vs.wsaf (plaf, obsWSAF, expWSAF)
plot.llk( llkTable, ref, alt, expWSAF, logFileName)



    dev.off()
}