fun.parse <- function( args ){
    fun.local.checkAndIncreaseArgI <- function ( ){
        arg_i = arg_i+1
    }

    outPrefix = "dataExplore"
    vcfFileName = ""
    refFileName = ""
    altFileName = ""
    plafFileName = ""
    excludeFileName = ""
    dEploidPrefix = ""
    excludeBool = FALSE
    ADFieldIndex = 2
    arg_i = 1
    while ( arg_i < length(args) ){
        argv = args[arg_i]
        if ( argv == "-vcf" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            vcfFileName = args[arg_i]
        } else if ( argv == "-plaf" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            plafFileName = args[arg_i]
        } else if ( argv == "-ref" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            refFileName = args[arg_i]
        } else if ( argv == "-alt" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            altFileName = args[arg_i]
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
        } else if ( argv == "-ADFieldIndex" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            ADFieldIndex = as.numeric(args[arg_i])
        } else {
            cat ("Unknow flag: ", argv, "\n")
        }

        arg_i = arg_i + 1
    }

#    if ( vcfFileName == "" || ( refFileName == "" && altFileName == "") ){
#        stop ("Vcf File name not specified!")
#    }

    if ( plafFileName == "" ){
        stop ("Plaf File name not specified!")
    }

#    cat ("vcfFileName: ", vcfFileName, "\n")
#    cat ("plafFileName: ", plafFileName, "\n")

    return ( list ( vcfFileName = vcfFileName,
                    refFileName = refFileName,
                    altFileName = altFileName,
                    plafFileName = plafFileName,
                    outPrefix = outPrefix,
                    dEploidPrefix = dEploidPrefix,
                    excludeFileName = excludeFileName,
                    excludeBool = excludeBool,
                    ADFieldIndex = ADFieldIndex) )
}


fun.dEploidPrefix <- function ( prefix ){
    if ( prefix == "" ){
        stop ("dEprefix ungiven!!!")
    }

    return ( list ( propFileName = paste(prefix, ".prop",   sep = ""),
                    hapFileName  = paste(prefix, ".hap",    sep = ""),
                    llkFileName  = paste(prefix, ".llk",    sep = ""),
                    dicLogFileName  = paste(prefix, "dic.log", sep = "") ) )
}


fun.extract.coverage <- function ( inputs ){
    if ( inputs$vcfFileName != "" ){
        return (extractCoverageFromVcf (inputs$vcfFileName, inputs$ADFieldIndex ))
    } else {
        return (extractCoverageFromTxt (inputs$refFileName, inputs$altFileName))
    }

}


fun.extract.exclude <- function (excludeFileName, excludeBool){
    if ( excludeBool ) {
        return ( list ( excludeBool = excludeBool,
                        excludeTable = read.table(excludeFileName, header = TRUE, comment.char = "")))
    } else {
        return ( list ( excludeBool = excludeBool ))
    }
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


plot.llk <- function (llkTable, ref, alt, expWSAF, title = "", cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    llk = llkTable$V2
    llkEvent = llkTable$V1
#    llk_sd = sd(llk)
#    llk_range = range(llk)

#    dic.by.var = fun.dic.by.llk.var (llk)
#    dic.by.theta = fun.dic.by.theta ( llk, fun.llk(ref, alt, expWSAF))

    plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="LLK", main=title,
        cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    updateSingleAt = which(llkEvent == 1)
    updateBothAt = which(llkEvent == 2)
    updatePropAt = which(llkEvent == 0)
    index = c(1:length(llk))
    points(index[updateSingleAt], llk[updateSingleAt], cex = 0.6, col="red")
    points(index[updateBothAt], llk[updateBothAt], cex = 0.6, col="blue")
    points(index[updatePropAt], llk[updatePropAt], cex = 0.6, col="green")

}


fun.getllk.dic <- function ( llkTable, ref, alt, expWSAF, logFileName ){
    llk = llkTable$V2
    llkEvent = llkTable$V1
    llk_sd = sd(llk)
    llk_range = range(llk)
    dic.by.var = fun.dic.by.llk.var (llk)
    dic.by.theta = fun.dic.by.theta ( llk, fun.llk(ref, alt, expWSAF))

    cat ( "dic.by.var: ", dic.by.var, "\n", file = logFileName, append = T)
    cat ( "dic.by.theta: ", dic.by.theta, "\n", file = logFileName, append = T)
    return (paste("LLK sd:", round(llk_sd, digits = 4),
                  ",\ndic.by.var: ",round(dic.by.var),
                  ", dic.by.theta: ",round(dic.by.theta)))
}


fun.getWSAF.corr <- function( obsWSAF, expWSAF, dicLogFileName ){
    currentWSAFcov  = cov(obsWSAF, expWSAF)
    currentWSAFcorr = cor(obsWSAF, expWSAF)
    cat ( "corr: ",  currentWSAFcorr, "\n", file = dicLogFileName)

    return (paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"))
}


fun.dataExplore <- function (coverage, PLAF, prefix = "") {
#    PLAF = plafInfo$PLAF
    ref = coverage$refCount
    alt = coverage$altCount

    png ( paste ( prefix, "altVsRefAndWSAFvsPLAF.png", sep = "" ), width = 1800, height = 600)
    par(mar = c(5,7,7,4))
    par( mfrow = c(1,3) )

    plotAltVsRef ( ref, alt, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    obsWSAF = computeObsWSAF ( alt, ref )

    histWSAF ( obsWSAF, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    plotWSAFvsPLAF ( PLAF, obsWSAF, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    dev.off()
}


fun.interpretDEploid.1 <- function (coverage, PLAF, dEploidPrefix, prefix = "", exclude ) {

#    PLAF = plafInfo$PLAF
    ref = coverage$refCount
    alt = coverage$altCount

    dEploidOutput = fun.dEploidPrefix ( dEploidPrefix )
    tmpProp = read.table(dEploidOutput$propFileName, header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
    hap = as.matrix(read.table(dEploidOutput$hapFileName, header=T)[,-c(1,2)] )
    expWSAF = hap %*%prop
    llkTable = read.table( dEploidOutput$llkFileName, header=F)

#    png ( paste ( prefix, ".interpretDEploidFigure.1.png", sep = "" ),  width = 1500, height = 500, bg="transparent")
    png ( paste ( prefix, ".interpretDEploidFigure.1.png", sep = "" ),  width = 1500, height = 1000)
    par(mar = c(5,7,7,4))
    par( mfrow = c(2,3) )
    plotAltVsRef ( ref, alt, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    obsWSAF = computeObsWSAF ( alt, ref )
    histWSAF ( obsWSAF, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    if (exclude$excludeBool){
        excludeLogic = ( paste(coverage$CHROM, coverage$POS) %in% paste(exclude$excludeTable$CHROM, exclude$excludeTable$POS) )
        excludeindex = which(excludeLogic)
        includeindex = which(!excludeLogic)
        obsWSAF = obsWSAF[includeindex]
        PLAF = PLAF[includeindex]
        ref = ref[includeindex]
        alt = alt[includeindex]
    }
    plotWSAFvsPLAF ( PLAF, obsWSAF, expWSAF, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    plotProportions( tmpProp, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    tmpTitle = fun.getWSAF.corr (obsWSAF, expWSAF, dEploidOutput$dicLogFileName)
    plotObsExpWSAF ( obsWSAF, expWSAF, tmpTitle, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    tmpTitle = fun.getllk.dic (llkTable, ref, alt, expWSAF, dEploidOutput$dicLogFileName )
    plot.llk( llkTable, ref, alt, expWSAF, tmpTitle, cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5 )

    dev.off()
}


plot.wsaf.vs.index <- function ( coverage, expWSAF = c(), expWSAFChrom = c(), exclude, titlePrefix = "" ){
    chromList = levels(coverage$CHROM)
    ref = coverage$refCount
    alt = coverage$altCount
    obsWSAF = computeObsWSAF ( alt, ref )
    nFigures = length(chromList)
    for ( chromI in chromList ){
        plot( obsWSAF[coverage$CHROM==chromI], col="red", ylim=c(0,1), main = paste(titlePrefix, chromI, "WSAF"), ylab = "WSAF", cex.axis = 2*nFigures/8, cex.lab = 2*nFigures/8,
        cex.main = 2*nFigures/6)

        if ( length(expWSAF) > 0 ){
            plotIndex = c()
            if (exclude$excludeBool){
                tmpCoveragePos = coverage$POS[coverage$CHROM==chromI]
                tmpExcludePos = exclude$excludeTable$POS[exclude$excludeTable$CHROM==chromI]
                excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
                excludeindex = which(excludeLogic)
                plotIndex = which(!excludeLogic)
            } else {
                plotIndex = c(1:length(obsWSAF[coverage$CHROM==chromI]))
            }
            points(plotIndex, expWSAF[expWSAFChrom == chromI], col="blue")
        }
    }
}


fun.interpretDEploid.2 <- function ( coverage, dEploidPrefix, prefix = "", exclude ){
    dEploidOutput = fun.dEploidPrefix ( dEploidPrefix )
    tmpProp = read.table(dEploidOutput$propFileName, header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
    hapInfo = read.table(dEploidOutput$hapFileName, header=T)
    hapChrom = hapInfo[,1]
    hap = as.matrix(hapInfo[,-c(1,2)])
    expWSAF = hap %*%prop

    png(paste( prefix, ".interpretDEploidFigure.2.png", sep= ""), width = 3500, height = 2000)
    chromName = levels(coverage$CHROM)
    ncol = ceiling(length(chromName)/2)
    par(mar = c(5,7,7,4))
    par(mfrow = c(ncol,length(chromName)/ncol))
    plot.wsaf.vs.index ( coverage, expWSAF, hapChrom, exclude )
    dev.off()

}


plot.postProb.ofCase <- function ( inPrefix, outPrefix, case, strainNumber ){
    png(paste(outPrefix, ".", case, ".png", sep = ""), width = 3500, height = 2000)
    obj = read.table( paste(inPrefix, ".", case, sep = ""), header=T)
    chromName = levels(obj$CHROM)
    ncol = ceiling(length(chromName)/2)
    par(mfrow = c(ncol,length(chromName)/ncol))
    par(mar = c(5,7,7,4))
    nFigures = length(chromName)
    for ( chromI in chromName ){
        haplotypePainter ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])], paste("Strain", strainNumber+1, chromI, "posterior probabilities"), 2*nFigures)
    }
    dev.off()
}


fun.interpretDEploid.3 <- function ( inPrefix, outPrefix = "" ){
    strainI = 0
    while ( file.exists(paste(inPrefix, ".single", strainI, sep="")) ){
        plot.postProb.ofCase( inPrefix, outPrefix, paste("single", strainI, sep=""), strainI)
        strainI = strainI+1
    }
}
