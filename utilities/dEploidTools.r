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
                    excludeBool = excludeBool) )
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
        return (extractCoverageFromVcf (inputs$vcfFileName))
    } else {
        return (extractCoverageFromTxt (inputs$refFileName, inputs$altFileName))
    }

}

#' @title Extract read counts from plain text file
#'
#' @description Extract read counts from tab-delimited text files of a single sample.
#'
#' @note The allele count files must be tab-delimited. The allele count files contain three columns: chromosomes, positions and allele count.
#'
#' @param refFileName Path of the reference allele count file.
#'
#' @param altFileName Path of the alternative allele count file.
#'
#' @return A data.frame contains four columns: chromosomes, positions, reference allele count, alternative allele count.
#'
#' @export
#'
#' @examples
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390 = extractCoverageFromTxt(refFile, altFile)
#'
extractCoverageFromTxt <- function ( refFileName, altFileName ){
    ref = read.table(refFileName, header = TRUE, comment.char = "")
    alt = read.table(altFileName, header = TRUE, comment.char = "")
    return ( data.frame( CHROM = ref[,1],
                         POS = ref[,2],
                         refCount = ref[,3],
                         altCount = alt[,3] )
           )
}


fun.extract.exclude <- function (excludeFileName, excludeBool){
    if ( excludeBool ) {
        return ( list ( excludeBool = excludeBool,
                        excludeTable = read.table(excludeFileName, header = TRUE, comment.char = "")))
    } else {
        return ( list ( excludeBool = excludeBool ))
    }
}

#' @title Extract read counts from VCF
#'
#' @description Extract read counts from VCF file of a single sample.
#'
#' @note The VCF file should only contain one sample. If more samples present in the VCF, it only returns coverage for of the first sample.
#'
#' @param vcfName Path of the VCF file.
#'
#' @param ADFieldIndex Index of the AD field of the sample field. For example, if the format is "GT:AD:DP:GQ:PL", the AD index is 2 (by default).
#'
#' @return A data.frame contains four columns: chromosomes, positions, reference allele count, alternative allele count.
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#'
extractCoverageFromVcf <- function ( vcfName, ADFieldIndex = 2 ){
    # Assume that AD is the second field
    catCmd = "cat"
    if ( grepl("gzip", system(paste("file --mime-type", vcfName), T) ) == TRUE ){
        catCmd = "zcat"
    }

    skipNum = as.numeric(system(paste(catCmd, vcfName, " | head -5000 | grep \"##\" | wc -l"), T))
    vcf  = read.table( gzfile(vcfName), skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

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


#' @title Extract PLAF
#'
#' @description Extract population level allele frequency (PLAF) from text file.
#'
#' @note The text file must have header, and population level allele frequency recorded in the "PLAF" field.
#'
#' @param plafName Path of the PLAF text file.
#'
#' @return A numeric array of PLAF
#'
#' @export
#'
#' @examples
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' plaf = extractPLAF(plafFile)
#'
extractPLAF<- function ( plafName ){
    return ( read.table(plafName, header=T)$PLAF )
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



#' @title Plot proportions
#'
#' @description Plot the MCMC samples of the proportion, indexed by the MCMC chain.
#'
#' @param proportions Matrix of the MCMC proportion samples. The matrix size is number of the MCMC samples by the number of strains.
#'
#' @param title Figure title.
#'
#' @return
#'
#' @export
#'
#' @examples
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' panelFile = system.file("extdata", "labStrains.test.panel.txt", package = "DEploid")
#' PG0390 = dEploid(paste("-vcf", vcfFile, "-plaf", plafFile, "-noPanel"))
#' plotProportions( PG0390$Proportions, "PG0390-C proportions" )
#'
plotProportions <-function (proportions, title = "Components"){
    rainbowColorBin = 16
    barplot(t(proportions), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="Iteration", ylab="Component proportion", main=title)
}


#' @title Plot coverage
#'
#' @description Plot alternative allele count vs reference allele count at each site.
#'
#' @param ref Numeric array of reference allele count.
#'
#' @param alt Numeric array of alternative allele count.
#'
#' @param title Figure title, "Alt vs Ref" by default
#'
#' @param exclude.ref Numeric array of reference allele count at sites that are not deconvoluted.
#'
#' @param exclude.alt Numeric array of alternative allele count at sites that are not deconvoluted
#'
#' @return
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#' plotAltVsRef( PG0390$refCount, PG0390$altCount )
#'
plotAltVsRef <- function ( ref, alt, title = "Alt vs Ref", exclude.ref = c(), exclude.alt = c() ){
    tmp.range = 1.1*mean(max(alt), max(ref))
    plot ( ref, alt, xlim=c(0, tmp.range), ylim=c(0,tmp.range), cex = 0.5, xlab = "REF", ylab = "ALT", main = title)
    points (exclude.ref, exclude.alt, col = "red")
    abline(v =50, untf = FALSE, lty = 2)
    abline(h =50, untf = FALSE, lty = 2)

    abline(h =150, untf = FALSE, lty = 2)
    abline(v =150, untf = FALSE, lty = 2)
}


#' @title WSAF histogram
#'
#' @description Produce histogram of the allele frequency within sample.
#'
#' @param obsWSAF Observed allele frequency within sample
#'
#' @param exclusive When TRUE 0 < WSAF < 1; otherwise 0 <= WSAF <= 1.
#'
#' @param title Histogram title
#'
#' @return histogram
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390$altCount, PG0390$refCount )
#' histWSAF(obsWSAF)
#' myhist = histWSAF(obsWSAF, FALSE)
#'
histWSAF <- function ( obsWSAF, exclusive = TRUE, title ="Histogram 0<WSAF<1" ){
    tmpWSAF_index = 1:length(obsWSAF)
    if ( exclusive ){
        tmpWSAF_index = which(((obsWSAF<1) * (obsWSAF>0) ) == 1)
    }
    return (hist(obsWSAF[tmpWSAF_index], main=title, breaks = seq(0, 1, by =0.1), xlab = "WSAF"))
}


plot.llk <- function (llkTable, ref, alt, expWSAF, title = "" ){
    llk = llkTable$V2
    llkEvent = llkTable$V1
#    llk_sd = sd(llk)
#    llk_range = range(llk)

#    dic.by.var = fun.dic.by.llk.var (llk)
#    dic.by.theta = fun.dic.by.theta ( llk, fun.llk(ref, alt, expWSAF))

    plot(llk, lty=2, type="l", col="black", xlab="Iteration", ylab="LLK", main=title);
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
                  ", dic.by.var: ",round(dic.by.var),
                  ", dic.by.theta: ",round(dic.by.theta)))
}


#' @title Plot WSAF vs PLAF
#'
#' @description Plot allele frequencies within sample against population level.
#'
#' @param plaf Numeric array of population level allele frequency.
#'
#' @param obsWSAF Numeric array of observed altenative allele frequencies within sample.
#'
#' @param expWSAF Numeric array of expected WSAF from model.
#'
#' @param title Figure title, "WSAF vs PLAF" by default
#'
#' @return
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390$altCount, PG0390$refCount )
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' plaf = extractPLAF(plafFile)
#' plotWSAFvsPLAF(plaf, obsWSAF)
#'
plotWSAFvsPLAF <- function ( plaf, obsWSAF, expWSAF = c(), title = "WSAF vs PLAF" ){
    plot ( plaf, obsWSAF, cex = 0.5, xlim = c(0,1), ylim = c(0,1), col = "red", main = title, xlab = "PLAF", ylab = "WSAF" )
    if ( length(expWSAF) > 0 ){
        points ( plaf, expWSAF, cex = 0.5, col = "blue")
    }
}


fun.getWSAF.corr <- function( obsWSAF, expWSAF, dicLogFileName ){
    currentWSAFcov  = cov(obsWSAF, expWSAF)
    currentWSAFcorr = cor(obsWSAF, expWSAF)
    cat ( "corr: ",  currentWSAFcorr, "\n", file = dicLogFileName)

    return (paste("Sample Freq (cov =",format(currentWSAFcov,digits=4), "corr = ", format(currentWSAFcorr,digits=4),")"))
}


#' @title Plot WSAF
#'
#' @description Plot observed alternative allele frequency within sample against expected WSAF.
#'
#' @param obsWSAF Numeric array of observed WSAF.
#'
#' @param expWSAF Numeric array of expected WSAF.
#'
#' @param title Figure title.
#'
#' @return
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390$altCount, PG0390$refCount )
#' PG0390.deconv = dEploid(paste("-vcf", vcfFile, "-plaf", plafFile, "-noPanel"))
#' prop = PG0390.deconv$Proportions[dim(PG0390.deconv$Proportions)[1],]
#' expWSAF = t(PG0390.deconv$Haps) %*% prop
#' plotObsExpWSAF(obsWSAF, expWSAF)
#'
plotObsExpWSAF <- function (obsWSAF, expWSAF, title = "WSAF(observed vs expected)"){
    plot(obsWSAF, expWSAF, pch=19, col="blue", xlab="Observed WSAF (ALT/(ALT+REF))", ylab="Expected WSAF (h%*%p)",
         main=title,
         xlim = c(-0.05, 1.05), cex = 0.5, ylim = c(-0.05, 1.05));
    abline(0,1,lty="dotted");

}


#' @title Compute observed WSAF
#'
#' @description Compute observed allele frequency within sample from the allele counts.
#'
#' @param ref Numeric array of reference allele count.
#'
#' @param alt Numeric array of alternative allele count.
#'
#' @return Numeric array of observed allele frequency within sample.
#'
#' @seealso \code{\link{histWSAF}} for histogram, \code{\link{plotWSAF}} for WSAFs plotted against indices.
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390 = extractCoverageFromVcf(vcfFile)
#' computeObsWSAF( PG0390$altCount, PG0390$refCount )
#'
computeObsWSAF <- function (alt, ref) {
    return ( alt / (ref + alt + 0.00000001) )
}


fun.dataExplore <- function (coverage, PLAF, prefix = "") {
#    PLAF = plafInfo$PLAF
    ref = coverage$refCount
    alt = coverage$altCount

    png ( paste ( prefix, "altVsRefAndWSAFvsPLAF.png", sep = "" ), width = 1800, height = 600)
    par( mfrow = c(1,3) )

    plotAltVsRef ( ref, alt )

    obsWSAF = computeObsWSAF ( alt, ref )

    histWSAF ( obsWSAF )

    plotWSAFvsPLAF ( PLAF, obsWSAF )

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
    par( mfrow = c(2,3) )
    plotAltVsRef ( ref, alt )

    obsWSAF = computeObsWSAF ( alt, ref )
    histWSAF ( obsWSAF )

    if (exclude$excludeBool){
        excludeLogic = ( paste(coverage$CHROM, coverage$POS) %in% paste(exclude$excludeTable$CHROM, exclude$excludeTable$POS) )
        excludeindex = which(excludeLogic)
        includeindex = which(!excludeLogic)
        obsWSAF = obsWSAF[includeindex]
        PLAF = PLAF[includeindex]
        ref = ref[includeindex]
        alt = alt[includeindex]
    }
    plotWSAFvsPLAF ( PLAF, obsWSAF, expWSAF )

    plotProportions( tmpProp )

    tmpTitle = fun.getWSAF.corr (obsWSAF, expWSAF, dEploidOutput$dicLogFileName)
    plotObsExpWSAF ( obsWSAF, expWSAF, tmpTitle )

    tmpTitle = fun.getllk.dic (llkTable, ref, alt, expWSAF, dEploidOutput$dicLogFileName )
    plot.llk( llkTable, ref, alt, expWSAF, tmpTitle )

    dev.off()
}


plot.wsaf.vs.index <- function ( coverage, expWSAF = c(), expWSAFChrom = c(), exclude, titlePrefix = "" ){
    chromList = levels(coverage$CHROM)
    ref = coverage$refCount
    alt = coverage$altCount
    obsWSAF = computeObsWSAF ( alt, ref )

    for ( chromI in chromList ){
        plot( obsWSAF[coverage$CHROM==chromI], col="red", ylim=c(0,1), main = paste(titlePrefix, chromI, "WSAF"))

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
    par(mfrow = c(ncol,length(chromName)/ncol))
    plot.wsaf.vs.index ( coverage, expWSAF, hapChrom, exclude)
    dev.off()

}

#' @title Painting haplotype according the reference panel
#'
#' @description Plot the posterior probabilities of a haplotype given the refernece panel.
#'
#' @param posteriorProbabilities Posterior probabilities matrix with the size of number of loci by the number of reference strain.
#'
#' @param title Figure title.
#'
#' @return
#'
#' @export
#'
#' @examples
#' to do ...
#'
haplotypePainter <-function (posteriorProbabilities, title = ""){
    rainbowColorBin = 16
    barplot(t(posteriorProbabilities), beside=F, border=NA, col=rainbow(rainbowColorBin), space=0, xlab="SNP index", ylab="Posterior probabilities", main=title)
}


plot.postProb.ofCase <- function ( inPrefix, outPrefix, case ){
    png(paste(outPrefix, ".", case, ".png", sep = ""), width = 3500, height = 2000)
    obj = read.table( paste(inPrefix, ".", case, sep = ""), header=T)
    chromName = levels(obj$CHROM)
    ncol = ceiling(length(chromName)/2)
    par(mfrow = c(ncol,length(chromName)/ncol))
    for ( chromI in chromName ){
        haplotypePainter ( obj[which( chromI == obj$CHROM),c(3:dim(obj)[2])], "")
    }
    dev.off()
}


fun.interpretDEploid.3 <- function ( inPrefix, outPrefix = "" ){
    strainI = 0
    while ( file.exists(paste(inPrefix, ".single", strainI, sep="")) ){
        plot.postProb.ofCase( inPrefix, outPrefix, paste("single", strainI, sep=""))
        strainI = strainI+1
    }
}
