fun.print.help <- function(){
    cat("    DEploid R utilities help\n")

    cat("
          Arguments:
               -help  --  Help. List the following content.
            -vcf STR  --  VCF file path.
            -ref STR  --  File path of reference allele count.
            -alt STR  --  File path of alternative allele count.
           -plaf STR  --  File path of population level allele frequencies.
        -exclude STR  --  File path of sites to be excluded.
              -o STR  --  Specify the file name prefix of the output.")
}

fun.print.help.interpret <- function(){
    fun.print.help()
    cat("
       -dEprefix STR  --  Specify DEploid output file prefix.\n")
    cat("
          Example:
          ./dEploid -vcf data/testData/PG0390-C.test.vcf \\
           -plaf data/testData/labStrains.test.PLAF.txt \\
           -o PG0390-CNopanel \\
           -noPanel
          ./utilities/dataExplore.r \\
          -vcf data/testData/PG0390-C.test.vcf.gz \\
          -plaf data/testData/labStrains.test.PLAF.txt \\
          -o PG0390-CNopanel \\
          -dEprefix PG0390-CNopanel\n\n")


    q(save="no")
}

fun.print.help.explore <- function(){
    fun.print.help()
    cat("\n
          Example:
            ./utilities/dataExplore.r \\
             -vcf data/testData/PG0390-C.test.vcf.gz \\
             -plaf data/testData/labStrains.test.PLAF.txt \\
             -o PG0390-C\n\n")
    q(save="no")
}

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
    inbreedingBool = FALSE
    ADFieldIndex = 2
    pdfBool = FALSE
    skip1Bool = FALSE
    ibdBool = FALSE
    helpBool = FALSE
    ringBool = FALSE
    ringDecreasingOrder = TRUE
    transformP = FALSE
    dEploid_v = "classic"
    arg_i = 1
    filter.window = 10
    filter.threshold = 0.995
    trackHeight = 0.6
    while ( arg_i <= length(args) ){
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
        } else if ( argv == "-inbreeding" ){
            inbreedingBool = TRUE
        } else if ( argv == "-ADFieldIndex" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            ADFieldIndex = as.numeric(args[arg_i])
        } else if ( argv == "-filter.threshold" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            filter.threshold = as.numeric(args[arg_i])
            # check range
            if ( filter.threshold < 0 | filter.threshold > 1){
                stop(paste("filter.threshold out of range [0, 1]:", filter.threshold))
            }
        } else if ( argv == "-filter.window" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            filter.window = as.numeric(args[arg_i])
            # check range
            if ( filter.window < 0 ){
                stop(paste("filter.window cannot be negative:", filter.window))
            }
        } else if ( argv == "-pdf" ){
            pdfBool = TRUE
        } else if ( argv == "-skip1" ){
            skip1Bool = TRUE
        } else if ( argv == "-ibd" ){
            ibdBool = TRUE
        } else if ( argv == "-help" ){
            helpBool = TRUE
        } else if ( argv == "-ring" ){
            ringBool = TRUE
        } else if ( argv == "-reverseRing" ){
            ringBool = TRUE
            ringDecreasingOrder = FALSE
        } else if ( argv == "-best" ){
            dEploid_v = "best"
        } else if ( argv == "-trackHeight" ){
            arg_i = fun.local.checkAndIncreaseArgI ( )
            trackHeight = as.numeric(args[arg_i])
            # check range
            if ( trackHeight < 0 | trackHeight > 1){
                stop(paste("trackHeight out of range [0, 1]:", trackHeight))
            }
        } else if ( argv == "-transformP" ){
            transformP = TRUE
        } else {
            stop(paste("Unknown flag:", argv))
        }

        arg_i = arg_i + 1
    }

    if (length(args) == 0){
        helpBool = TRUE
    }

#    if ( vcfFileName == "" || ( refFileName == "" && altFileName == "") ){
#        stop ("Vcf File name not specified!")
#    }
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
                    ADFieldIndex = ADFieldIndex,
                    pdfBool = pdfBool,
                    inbreedingBool = inbreedingBool,
                    skip1Bool = skip1Bool,
                    ibdBool = ibdBool,
                    helpBool = helpBool,
                    filter.threshold = filter.threshold,
                    filter.window = filter.window,
                    ringBool = ringBool,
                    ringDecreasingOrder = ringDecreasingOrder,
                    trackHeight = trackHeight,
                    transformP = transformP,
                    dEploid_v = dEploid_v) )
}


fun.dEploidPrefix <- function (prefix, dEploid_v = "classic" ){
    if ( prefix == "" ){
        stop ("dEprefix ungiven!!!")
    }

    if (dEploid_v == "classic") {
        version_suffix = "classic"
        return ( list ( propFileName = paste(prefix, ".", version_suffix, ".prop",   sep = ""),
                        hapFileName  = paste(prefix, ".", version_suffix, ".hap",    sep = ""),
                        llkFileName  = paste(prefix, ".", version_suffix, ".llk",    sep = "")
    #                    dicLogFileName  = paste(prefix, "dic.log", sep = "")
                     ) )
    } else if (dEploid_v == "best") {
    return ( list ( propFileName.chooseK = paste(prefix, ".chooseK.prop",   sep = ""),
                    hapFileName.chooseK  = paste(prefix, ".chooseK.hap",    sep = ""),
                    llkFileName.chooseK  = paste(prefix, ".chooseK.llk",    sep = ""),
                    propFileName.ibd = paste(prefix, ".ibd.prop",   sep = ""),
                    hapFileName.ibd = paste(prefix, ".ibd.hap",    sep = ""),
                    llkFileName.ibd = paste(prefix, ".ibd.llk",    sep = ""),
                    hapFileName.final = paste(prefix, ".final.hap",    sep = "")
#                    dicLogFileName  = paste(prefix, "dic.log", sep = "")
                    ) )
    } else {
        return(NULL)
    }
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

#    cat ( "dic.by.var: ", dic.by.var, "\n", file = logFileName, append = T)
#    cat ( "dic.by.theta: ", dic.by.theta, "\n", file = logFileName, append = T)
    return (paste("LLK sd:", round(llk_sd, digits = 4),
                  ",\ndic.by.var: ",round(dic.by.var),
                  ", dic.by.theta: ",round(dic.by.theta)))
}


fun.getWSAF.corr <- function( obsWSAF, expWSAF, dicLogFileName ){
    currentWSAFcov  = cov(obsWSAF, expWSAF)
    currentWSAFcorr = cor(obsWSAF, expWSAF)
#    cat ( "corr: ",  currentWSAFcorr, "\n", file = dicLogFileName)

    return (paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"))
}


fun.find.more <- function (outliers.idx, window.size){
    idx.out = c()
    for ( i in 1:length(outliers.idx)){
        near.outliers.idx = which(((outliers.idx[i] - window.size) < outliers.idx) & (outliers.idx < (outliers.idx[i] + window.size)))
        idx.len = length(near.outliers.idx)
        if ( length(near.outliers.idx)>1 ){
            idx.out = c(idx.out, outliers.idx[near.outliers.idx[1]]:outliers.idx[near.outliers.idx[idx.len]])
        } else{
            idx.out = c(idx.out, outliers.idx[near.outliers.idx[1]])
        }
    }
    return(unique(idx.out))
}


plot.total.coverage <- function(ref, alt, chroms, cex.lab = 1, cex.main = 1, cex.axis = 1,  threshold, window.size){
    totalDepth = ref + alt
    x = 1:length(totalDepth)
    tmpQ = quantile(totalDepth, threshold)
    tmpIdx = which((totalDepth > tmpQ ))
    potentialOutliers = fun.find.more(tmpIdx, window.size)

    chromCol = (as.numeric(chroms) %% 2 )
    chromCol[chromCol==1] = NA
    chromCol[chromCol==0] = 8
    plot(x, totalDepth, type="n", cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main, ylab="Coverage depth", xlab="SNP index", main = "Coverage across the sequence")
    rect(x[-1],
         0,
         x[-length(x)],
         max(totalDepth)*1.5, col = chromCol, border = "transparent")
    points(x, totalDepth, pch = 16)
    abline(h = tmpQ, col = "red")
    points(x[potentialOutliers], totalDepth[potentialOutliers], col = "red", pch = "x", cex = 2)
    return (potentialOutliers)
}


fun.dataExplore <- function (coverage, PLAF, prefix = "", pdfBool, threshold = 0.995, window.size = 10) {
#    PLAF = plafInfo$PLAF
    ref = coverage$refCount
    alt = coverage$altCount

    if ( pdfBool == TRUE ){
        cexSize = 3
        pdf ( paste ( prefix, "altVsRefAndWSAFvsPLAF.pdf", sep = "" ), width = 30, height = 20)
    } else {
        cexSize = 2.5
        png ( paste ( prefix, "altVsRefAndWSAFvsPLAF.png", sep = "" ), width = 1800, height = 1500)
    }

    layout(matrix(c(1,1,1,
                    1,1,1,
                    2,3,4,
                    2,3,4,
                    5,5,5), 5, 3, byrow = TRUE))
    par(mar = c(5,7,7,4))

    badGuys = plot.total.coverage(ref, alt, coverage$CHROM, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize, threshold, window.size)

    if ( length(badGuys) > 0 ){
        CHROM = coverage$CHROM[badGuys]
        POS = coverage$POS[badGuys]
        write.table(data.frame(CHROM, POS), file = paste(prefix, "PotentialOutliers.txt", sep=""), sep = "\t", quote = F, row.names = F)
    }

    plotAltVsRef ( ref, alt, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize, potentialOutliers = badGuys )

    obsWSAF = computeObsWSAF ( alt, ref )

    histWSAF ( obsWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    plotWSAFvsPLAF ( PLAF, obsWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize, potentialOutliers = badGuys  )

    plot(obsWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize, ylab = "WSAF", type = "n")
    chromSize = table(coverage$CHROM)
    pf_x = c(0, cumsum(chromSize))
    pf_chromCol = (as.numeric(as.factor(names(chromSize))) %% 2 )
    pf_chromCol[pf_chromCol==1] = NA
    pf_chromCol[pf_chromCol==0] = 8
    rect(pf_x[-1],
         0,
         pf_x[-length(pf_x)],
         1, col = pf_chromCol)
    points(obsWSAF, cex = 0.5)
    dev.off()
}


fun.interpretDEploid.best <- function (coverage, PLAF, dEploidPrefix, prefix = "", exclude, pdfBool ) {
    ref = coverage$refCount
    alt = coverage$altCount

    dEploidOutput = fun.dEploidPrefix(dEploidPrefix, dEploid_v = "best")
    cexSize = 2.5

    png ( paste ( prefix, ".interpretDEploidFigure.1.png", sep = "" ),  width = 1500, height = 3750)
    par(mar = c(5,7,7,4))
    par( mfrow = c(5,2) )
    ##################################################################################################
    hap.chooseK = read.table(dEploidOutput$hapFileName.chooseK, header=T)
    includeLogic.chooseK = which( paste(coverage$CHROM, coverage$POS) %in% paste(hap.chooseK$CHROM, hap.chooseK$POS) )
    plotAltVsRef ( ref[includeLogic.chooseK], alt[includeLogic.chooseK], cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )
    plotAltVsRef ( ref, alt, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    ##################################################################################################
    obsWSAF = computeObsWSAF ( alt, ref )
    histWSAF (obsWSAF[includeLogic.chooseK], cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )
    histWSAF (obsWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    if (exclude$excludeBool){
        excludeLogic = ( paste(coverage$CHROM, coverage$POS) %in% paste(exclude$excludeTable$CHROM, exclude$excludeTable$POS) )
        excludeindex = which(excludeLogic)
        includeindex = which(!excludeLogic)
        obsWSAF = obsWSAF[includeindex]
        PLAF = PLAF[includeindex]
        includeLogic.chooseK = which( paste(coverage$CHROM[includeindex], coverage$POS[includeindex]) %in% paste(hap.chooseK$CHROM, hap.chooseK$POS) )
    }

    ##################################################################################################
    pp = read.table(dEploidOutput$propFileName.ibd, header=F)
    p = read.table(dEploidOutput$propFileName.chooseK, sep ="\t", header=F)
    barplot(t(p), border=F, space=0)
    barplot(t(pp), border=F, space=0)

    ##################################################################################################
    prop = as.numeric(pp[dim(pp)[1],])
    hap = as.matrix(read.table(dEploidOutput$hapFileName.final, header=T)[,-c(1,2)] )
    expWSAF = hap %*% prop
    plotWSAFvsPLAF(PLAF[includeLogic.chooseK], obsWSAF[includeLogic.chooseK], expWSAF[includeLogic.chooseK], cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )
    plotWSAFvsPLAF(PLAF, obsWSAF, expWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    ##################################################################################################
    tmpTitle = fun.getWSAF.corr (obsWSAF[includeLogic.chooseK], expWSAF[includeLogic.chooseK], "")
    plotObsExpWSAF ( obsWSAF[includeLogic.chooseK], expWSAF[includeLogic.chooseK], tmpTitle, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )
    tmpTitle = fun.getWSAF.corr (obsWSAF, expWSAF, "")
    plotObsExpWSAF ( obsWSAF, expWSAF, tmpTitle, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )
    dev.off()
}


fun.interpretDEploid.1 <- function (coverage, PLAF, dEploidPrefix, prefix = "", exclude, pdfBool ) {

#    PLAF = plafInfo$PLAF
    ref = coverage$refCount
    alt = coverage$altCount

    dEploidOutput = fun.dEploidPrefix ( dEploidPrefix )
    tmpProp = read.table(dEploidOutput$propFileName, header=F)
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
    hap = as.matrix(read.table(dEploidOutput$hapFileName, header=T)[,-c(1,2)] )
    expWSAF = hap %*%prop
    llkTable = read.table( dEploidOutput$llkFileName, header=F)


    if ( pdfBool == TRUE ){
        cexSize = 3.5
        pdf ( paste ( prefix, ".interpretDEploidFigure.1.pdf", sep = "" ), width = 30, height = 20)
    } else {
        cexSize = 2.5
        png ( paste ( prefix, ".interpretDEploidFigure.1.png", sep = "" ),  width = 1500, height = 1000)
    }

    par(mar = c(5,7,7,4))
    par( mfrow = c(2,3) )
    plotAltVsRef ( ref, alt, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    obsWSAF = computeObsWSAF ( alt, ref )
    histWSAF ( obsWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    if (exclude$excludeBool){
        excludeLogic = ( paste(coverage$CHROM, coverage$POS) %in% paste(exclude$excludeTable$CHROM, exclude$excludeTable$POS) )
        excludeindex = which(excludeLogic)
        includeindex = which(!excludeLogic)
        obsWSAF = obsWSAF[includeindex]
        PLAF = PLAF[includeindex]
        ref = ref[includeindex]
        alt = alt[includeindex]
    }
    plotWSAFvsPLAF ( PLAF, obsWSAF, expWSAF, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    plotProportions( tmpProp, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    tmpTitle = fun.getWSAF.corr (obsWSAF, expWSAF, dEploidOutput$dicLogFileName)
    plotObsExpWSAF ( obsWSAF, expWSAF, tmpTitle, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    tmpTitle = fun.getllk.dic (llkTable, ref, alt, expWSAF, dEploidOutput$dicLogFileName )
    plot.llk( llkTable, ref, alt, expWSAF, tmpTitle, cex.lab = cexSize, cex.main = cexSize, cex.axis = cexSize )

    dev.off()
}


plot.wsaf.vs.index.ring <- function ( coverage, expWSAF = c(), expWSAFChrom = c(), exclude, titlePrefix = "" ){
    chromCol = (as.numeric(1:length(levels(coverage$CHROM)) %% 2 ))
    chromCol[chromCol==1] = NA
    chromCol[chromCol==0] = 8

    circlize::circos.trackPlotRegion(factor = expWSAFChrom, ylim=c(0,1), track.height = 0.18, bg.col = chromCol, panel.fun=function(x,y){
        name = circlize::get.cell.meta.data("sector.index")
        xlim = circlize::get.cell.meta.data("xlim")
        ylim = circlize::get.cell.meta.data("ylim")
        chromRegion = coverage[coverage$CHROM==name,]

        ref = chromRegion$refCount
        alt = chromRegion$altCount
        obsWSAF = computeObsWSAF ( alt, ref )

        nSnp = dim(chromRegion)[1]

        if (exclude$excludeBool){
            tmpCoveragePos = coverage$POS[coverage$CHROM==name]
            tmpExcludePos = exclude$excludeTable$POS[exclude$excludeTable$CHROM==name]
            excludeLogic = ( tmpCoveragePos %in% tmpExcludePos )
            plotIndex = which(!excludeLogic)
        } else {
            plotIndex = c(1:nSnp)
        }
        circlize::circos.points(plotIndex, obsWSAF[plotIndex], col="red", pch = 16)
        circlize::circos.points(plotIndex, expWSAF[expWSAFChrom == name], col="blue", pch = 16)
    })
}


plot.wsaf.vs.index <- function ( coverage, expWSAF = c(), expWSAFChrom = c(), exclude, titlePrefix = "" ){
    chromList = unique(coverage$CHROM)
    ref = coverage$refCount
    alt = coverage$altCount
    obsWSAF = computeObsWSAF ( alt, ref )
    nFigures = length(chromList)
    totalCoverage = alt + ref
    for ( chromI in chromList ){
        print(chromI)
        tmpWSAF = obsWSAF[coverage$CHROM==chromI]
        colorFrac = (totalCoverage[coverage$CHROM==chromI])/50
        colorTrans = unlist(lapply(colorFrac, function(x){adjustcolor("red", alpha.f=x)}))
        plot(tmpWSAF , col=colorTrans, ylim=c(0,1), main = paste(titlePrefix, chromI, "WSAF"), ylab = "WSAF",
            cex.axis = 3.5, cex.lab = 3.5, cex.main = 4, xaxt = "n", yaxt = "n", pch = 16)
        newXaxt = round(seq(1, length(tmpWSAF), length.out = 6))
        axis(1, at = newXaxt, labels = as.character(newXaxt),
            cex.axis= 3.5)
        newYaxt = seq(0, 1, length.out = 3)
        axis(2, at = newYaxt, labels = as.character(newYaxt),
            cex.axis= 3.5)


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
#            print(length(plotIndex))
#            print(length(expWSAF))
#            print(length(expWSAF[expWSAFChrom == chromI]))
#            print("##########")
            colorTrans = unlist(lapply(colorFrac, function(x){adjustcolor("blue", alpha.f=x)}))[plotIndex]
            points(plotIndex, expWSAF[expWSAFChrom == chromI], col=colorTrans, pch = 16)
        }
    }
}


fun.interpretDEploid.2 <- function ( coverage, dEploidPrefix, prefix = "", exclude, pdfBool, ringBool = FALSE, dEploid_v = "classic"){
    dEploidOutput = fun.dEploidPrefix(dEploidPrefix, dEploid_v)
    if (dEploid_v == "classic") {
        tmpProp = read.table(dEploidOutput$propFileName, header=F)
        hapInfo = read.table(dEploidOutput$hapFileName, header=T)
    } else if (dEploid_v == "best") {
        tmpProp = read.table(dEploidOutput$propFileName.ibd, header=F)
        hapInfo = read.table(dEploidOutput$hapFileName.final, header=T)
    }
    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
    hapChrom = hapInfo[,1]
    hap = as.matrix(hapInfo[,-c(1,2)])
    expWSAF = hap %*%prop

    if ( ringBool ){
        if ( pdfBool == TRUE ){
            cexSize = 3
            pdf ( paste ( prefix, ".interpretDEploidFigure.2.ring.pdf", sep = "" ), width = 45, height = 45)
        } else {
            cexSize = 2.5
            png ( paste ( prefix, ".interpretDEploidFigure.2.ring.png", sep= ""), width = 3500, height = 3500)
        }
        fun.ring.plot.initialize (hapInfo$CHROM)
        plot.wsaf.vs.index.ring ( coverage, expWSAF, hapChrom, exclude )

        circlize::circos.clear();
        dev.off()

    } else {
        if ( pdfBool == TRUE ){
            cexSize = 3
            pdf ( paste ( prefix, ".interpretDEploidFigure.2.pdf", sep = "" ), width = 45, height = 30)
        } else {
            cexSize = 2.5
            png ( paste ( prefix, ".interpretDEploidFigure.2.png", sep= ""), width = 3500, height = 2000)
        }
        chromName = levels(coverage$CHROM)
        ncol = ceiling(length(chromName)/2)
        par(mar = c(5,7,7,4))
        par(mfrow = c(ncol,length(chromName)/ncol))
        plot.wsaf.vs.index ( coverage, expWSAF, hapChrom, exclude )
        dev.off()
    }
}


plot.postProb.ofCase <- function ( inPrefix, outPrefix, case, strainNumber, pdfBool, inbreeding = FALSE  ){
    if ( pdfBool == TRUE ){
        cexSize = 3
        if ( inbreeding ){
            pdf(paste(outPrefix, ".", case, ".inbreeding.pdf", sep = ""), width = 45, height = 30)
        } else {
            pdf(paste(outPrefix, ".", case, ".pdf", sep = ""), width = 45, height = 30)
        }
    } else {
        cexSize = 2.5
        if (inbreeding){
            png(paste(outPrefix, ".", case, ".inbreeding.png", sep = ""), width = 3500, height = 2000)
        } else{
            png(paste(outPrefix, ".", case, ".png", sep = ""), width = 3500, height = 2000)
        }
    }
    obj = read.table( paste(inPrefix, ".", case, sep = ""), header=T)
    chromName = levels(obj$CHROM)
    ncol = ceiling(length(chromName)/2)
    par(mfrow = c(ncol,length(chromName)/ncol))
    par(mar = c(5,7,7,4))
    nFigures = length(chromName)
    for ( chromI in chromName ){
        if ( inbreeding ){
            haplotypePainter ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])],
                title = paste("Strain", strainNumber+1, chromI, "inbreeding probabilities"),
                labelScaling = 2*nFigures,
                numberOfInbreeding = sum(grepl("I", names(obj))))
        } else{
            haplotypePainter ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])],
                title = paste("Strain", strainNumber+1, chromI, "posterior probabilities"),
                labelScaling = 2*nFigures)
        }

    }
    dev.off()
}


fun.interpretDEploid.3 <- function ( inPrefix, outPrefix = "", pdfBool, inbreeding = FALSE  ){
    strainI = 0
    while ( file.exists(paste(inPrefix, ".single", strainI, sep="")) ){
        plot.postProb.ofCase( inPrefix, outPrefix, paste("single", strainI, sep=""), strainI, pdfBool)
        if ( inbreeding ){
            plot.postProb.ofCase( inPrefix, outPrefix, paste("single", strainI, sep=""), strainI, pdfBool, inbreeding = TRUE)
        }
        strainI = strainI+1
    }
}


fun.interpretDEploid.4 <- function ( inPrefix, outPrefix = "", pdfBool ){
    inFile = paste(inPrefix, ".ibd.probs", sep = "")
    if (!file.exists(inFile)){
        print("In file not exist")
	return()
    }

    if ( pdfBool == TRUE ){
        cexSize = 3
        pdf(paste(outPrefix, ".ibd.probs.pdf", sep = ""), width = 45, height = 30)
    } else {
        cexSize = 2.5
        png(paste(outPrefix, ".ibd.probs.png", sep = ""), width = 3500, height = 2000)
    }
    obj = read.table(inFile , header=T)
    chromName = levels(obj$CHROM)
    ncol = ceiling(length(chromName)/2)
    par(mfrow = c(ncol,length(chromName)/ncol))
    par(mar = c(5,7,7,4))
    nFigures = length(chromName)
    for ( chromI in chromName ){
        haplotypePainter ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])],
            title = paste("IBD probabilities of", chromI),
            labelScaling = 2*nFigures)
    }
    dev.off()
}


fun.interpretDEploid.3.ring <- function (inPrefix, outPrefix = "", pdfBool, inbreeding = FALSE, coverage, exclude, ringDecreasingOrder, trackHeight = 0.8, transformP = FALSE ){
    if ( pdfBool == TRUE ){
        cexSize = 3
        if ( inbreeding ){
            pdf(paste(outPrefix, ".inbreeding.ring.pdf", sep = ""), width = 45, height = 45)
        } else {
            pdf(paste(outPrefix, ".ring.pdf", sep = ""), width = 45, height = 45)
        }
    } else {
        cexSize = 2.5
        if (inbreeding){
            png(paste(outPrefix, ".inbreeding.ring.png", sep = ""), width = 3500, height = 3500)
        } else{
            png(paste(outPrefix, ".ring.png", sep = ""), width = 3500, height = 3500)
        }
    }

    dEploidOutput = fun.dEploidPrefix ( inPrefix )
    tmpProp = read.table(dEploidOutput$propFileName, header=F)

    lastProp = as.numeric(tmpProp[dim(tmpProp)[1],])

    orderedProp = sort.int(lastProp, index.return=T, decreasing=ringDecreasingOrder)
    orderedProp.p = orderedProp$x
    myOrder = orderedProp$ix-1

    myOrder = myOrder[orderedProp.p>0.01]
    orderedProp.p = orderedProp.p[orderedProp.p>0.01]

    if (transformP){
        orderedProp.p = orderedProp.p %*% (diag(rep(1,length(orderedProp.p)))+1)
        orderedProp.p = orderedProp.p/sum(orderedProp.p)
        print(orderedProp.p)
    }

    first = myOrder[1]
    idx = 1
    for ( strain in myOrder ){
        readFrom = paste(inPrefix, ".single", strain, sep = "")
        if ( !file.exists(readFrom) ){
            next
        }

        cat("Loading from ", readFrom, " ")
        probs = read.table(readFrom, header=T)

        if ( strain == first ){
            fun.ring.plot.initialize(probs$CHROM)
            if ( inbreeding ){
                hapInfo = read.table(dEploidOutput$hapFileName, header=T)
                hapChrom = hapInfo[,1]
                hap = as.matrix(hapInfo[,-c(1,2)])
                expWSAF = hap %*% lastProp

                plot.wsaf.vs.index.ring ( coverage, expWSAF, hapChrom, exclude )

            }
        }

        circlize::circos.trackPlotRegion(factor = probs$CHROM, ylim=c(0,1), track.height = trackHeight*orderedProp.p[idx],
            panel.fun=function(x,y){
                name = circlize::get.cell.meta.data("sector.index")
                xlim = circlize::get.cell.meta.data("xlim")
                ylim = circlize::get.cell.meta.data("ylim")
                cat(".")
                chromRegion = probs[probs$CHROM==name,]
                if ( inbreeding == T ){
                    numberOfInbreeding = sum(grepl("I", names(probs)))
                    panelSize <- dim(probs)[2]-2-numberOfInbreeding
                    rainbowColors <- c(rep("#46a8e1", panelSize),
                                       rep("#f34747", numberOfInbreeding))
                } else {
                    rainbowColorBin <- 16
                    rainbowColors = rainbow(rainbowColorBin)

                }
                nSnp = dim(chromRegion)[1]
                print(nSnp)
                nhap = dim(chromRegion)[2] - 2
                cumProb = rep(0, nSnp)
                for ( i in 1:nhap ){
                    circlize::circos.rect(0:(nSnp-1), cumProb, 1:nSnp, cumProb + chromRegion[,i+2], col=rainbowColors[i], border=NA)
                    cumProb = cumProb + chromRegion[,i+2]
                }
            }
        )
        cat("\n")
        idx = idx+1
    }
    circlize::circos.clear();
    dev.off()
}


plot.ibd.change <- function(changeAt, titlePrefix,nFigures){
    plot(c(0, dim(changeAt)[1]), c(0, 1), type="n", ylim=c(0,.5), main = paste(titlePrefix, "IBD changes, and LS switches at"), ylab = "Frequency",
        cex.axis = 3.5, cex.lab = 3.5, cex.main = 4, xaxt = "n", yaxt = "n", xlab = "")
    lines(changeAt$IBDpathChangeAt, col = "red", lty=2, lwd=.5)
    lines(changeAt$finalIBDpathChangeAt, col = "red", lty=1, lwd=2)

#    lines(changeAt$siteOfTwoSwitchOne, col = "grey", lty=2)
#    lines(changeAt$finalSiteOfTwoSwitchOne, col = "grey", lty=1)

#    lines(changeAt$siteOfTwoMissCopyOne, col = "cyan", lty=2)
#    lines(changeAt$finalSiteOfTwoMissCopyOne, col = "cyan", lty=1)

    lines(changeAt$siteOfTwoSwitchTwo, col = "blue", lty=2, lwd=.5)
    lines(changeAt$finalSiteOfTwoSwitchTwo, col = "blue", lty=1, lwd=2)

#    lines(changeAt$siteOfTwoMissCopyTwo, col = "magenta", lty=2)
#    lines(changeAt$finalSiteOfTwoMissCopyTwo, col = "magenta", lty=1)

    lines(changeAt$siteOfOneSwitchOne, col = "green", lty=2, lwd=.5)
    lines(changeAt$finalSiteOfOneSwitchOne, col = "green", lty=1, lwd=2)

#    lines(changeAt$siteOfOneMissCopyOne, col = "yellow", lty=2)
#    lines(changeAt$finalSiteOfOneMissCopyOne, col = "yellow", lty=1)


    newXaxt = round(seq(1, dim(changeAt)[1], length.out = 6))
    axis(1, at = newXaxt, labels = as.character(newXaxt),
        cex.axis= 3.5)
    newYaxt = seq(0, 0.5, length.out = 3)
    axis(2, at = newYaxt, labels = as.character(newYaxt),
        cex.axis= 3.5)
}


#fun.interpretDEploid.4 <- function(inPrefix, outPrefix = "", pdfBool){
#    fileName = paste(inPrefix, ".extra", sep = "")
#    if ( file.exists(fileName) ){
#        obj = read.table(fileName , header=T)
#        chromName = levels(obj$CHROM)
#        ncol = ceiling(length(chromName)/2)
#        if ( pdfBool == TRUE ){
#            cexSize = 3
#            pdf ( paste ( outPrefix, ".interpretDEploidFigure.extra.pdf", sep = "" ), width = 45, height = 30)
#        } else {
#            cexSize = 2.5
#            png ( paste ( outPrefix, ".interpretDEploidFigure.extra.png", sep= ""), width = 3500, height = 2000)
#        }
#        par(mfrow = c(ncol,length(chromName)/ncol))
#        par(mar = c(5,7,7,4))
#        nFigures = length(chromName)
#        for ( chromI in chromName ){
#            plot.ibd.change ( obj[which( chromI == obj$CHROM),], chromI, nFigures )
#        }
#        dev.off()
#    }
#}


fun.ring.plot.initialize <- function(chrom, name.suffix = ""){
    circlize::circos.initialize(factor=chrom, xlim = cbind(1, table(chrom)))
    circlize::circos.trackPlotRegion(factor = chrom, ylim=c(0,1), track.height = 0.1, bg.border = NA,
        panel.fun=function(x,y){
            name = circlize::get.cell.meta.data("sector.index")
            xlim = circlize::get.cell.meta.data("xlim")
            ylim = circlize::get.cell.meta.data("ylim")
            circlize::circos.text(mean(xlim), 0.9, paste(name, name.suffix), cex = 4, facing = "inside")
            circlize::circos.axis(h = "bottom", labels.cex=3, direction = "outside", major.at=xlim, minor.ticks=1, labels.away.percentage = 0.15)
        }
    )

}
