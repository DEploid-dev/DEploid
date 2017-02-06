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
    ref <- read.table(refFileName, header = TRUE, comment.char = "")
    alt <- read.table(altFileName, header = TRUE, comment.char = "")
    return ( data.frame( CHROM = ref[, 1],
                         POS = ref[, 2],
                         refCount = ref[, 3],
                         altCount = alt[, 3] )
           )
}


#' @title Extract read counts from VCF
#'
#' @description Extract read counts from VCF file of a single sample.
#'
#' @note The VCF file should only contain one sample. If more samples present in the VCF, it only returns coverage for of the first sample.
#'
#' @param vcfFileName Path of the VCF file.
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
extractCoverageFromVcf <- function ( vcfFileName, ADFieldIndex = 2 ){
    # Assume that AD is the second field
    h <- function(w){
         if ( any( grepl( "gzfile connection", w) ) )
         invokeRestart( "muffleWarning" )
    }

    gzf <- gzfile(vcfFileName, open = "rb")
    numberOfHeaderLines <- 0
    line <- withCallingHandlers( readLines(gzf, n = 1), warning = h)
    while ( length(line) > 0 ){
        if (grepl("##", line )){
            numberOfHeaderLines <- numberOfHeaderLines + 1
        } else {
            break
        }
        line <- withCallingHandlers( readLines(gzf, n = 1), warning = h)
    }
    close(gzf)

    vcf <- read.table( gzfile(vcfFileName), skip = numberOfHeaderLines,
        header = T, comment.char = "", stringsAsFactors = FALSE,
        check.names = FALSE)

    sampleName <- names(vcf)[10]

    tmp <- vcf[[sampleName]]
    field <- strsplit(as.character(tmp), ":")

    tmpCovStr <- unlist(lapply(field, `[[`, ADFieldIndex))
    tmpCov <- strsplit(as.character(tmpCovStr), ",")

    refCount <- as.numeric(unlist(lapply(tmpCov, `[[`, 1)))
    altCount <- as.numeric(unlist(lapply(tmpCov, `[[`, 2)))

    return ( data.frame( CHROM = vcf[, 1],
                         POS = vcf[, 2],
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
#' @param plafFileName Path of the PLAF text file.
#'
#' @return A numeric array of PLAF
#'
#' @export
#'
#' @examples
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' plaf = extractPLAF(plafFile)
#'
extractPLAF <- function ( plafFileName ){
    return ( read.table(plafFileName, header = T)$PLAF )
}


#' @title Plot proportions
#'
#' @description Plot the MCMC samples of the proportion, indexed by the MCMC chain.
#'
#' @param proportions Matrix of the MCMC proportion samples. The matrix size is number of the MCMC samples by the number of strains.
#'
#' @param title Figure title.
#'
#' @param cex.lab Label size.
#'
#' @param cex.main Title size.
#'
#' @param cex.axis Axis text size.
#'
#' @export
#'
#' @examples
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' panelFile = system.file("extdata", "labStrains.test.panel.txt", package = "DEploid")
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390CoverageTxt = extractCoverageFromTxt(refFile, altFile)
#' PG0390CoverageTxt.deconv = dEploid(paste("-ref", refFile, "-alt", altFile,
#'     "-plaf", plafFile, "-noPanel"))
#' plotProportions( PG0390CoverageTxt.deconv$Proportions, "PG0390-C proportions" )
#'
plotProportions <- function (proportions, title = "Components",
                       cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    rainbowColorBin <- 16
    barplot(t(proportions), beside = F, border = NA,
        col = rainbow(rainbowColorBin), space = 0, xlab = "Iteration",
        ylab = "Component proportion", main = title,
        cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
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
#' @param cex.lab Label size.
#'
#' @param cex.main Title size.
#'
#' @param cex.axis Axis text size.
#'
#' @export
#'
#' @examples
#' # Example 1
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390CoverageTxt = extractCoverageFromTxt(refFile, altFile)
#' plotAltVsRef( PG0390CoverageTxt$refCount, PG0390CoverageTxt$altCount )
#'
#' # Example 2
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390CoverageVcf = extractCoverageFromVcf(vcfFile)
#' plotAltVsRef( PG0390CoverageVcf$refCount, PG0390CoverageVcf$altCount )
#'
plotAltVsRef <- function ( ref, alt, title = "Alt vs Ref",
                    exclude.ref = c(), exclude.alt = c(),
                    cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    cr <- colorRampPalette(colors = c("#de2d26", "#2b8cbe"))
    colors <- cr(31)
    ratios <- ref / (ref + alt + 0.0000001)
    tmpRange <- 1.1 * mean(max(alt), max(ref))
    plot ( ref, alt, xlim = c(0, tmpRange), ylim = c(0, tmpRange),
        pch = 20, col = scales::alpha(colors[ceiling(ratios * 30) + 1], 0.7),
        xlab = "Reference # Reads", ylab = "Alternative # Reads", main = title,
        cex = 0.5, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    legend("topright", legend = c("100% Alt", "100% Ref", "50/50"),
        fill = colors[c(1, 31, 15)], cex = cex.lab, border = NA, box.lwd = 0,
        box.col = "white", bg = NA)
    abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
    points(exclude.ref, exclude.alt, col = "red")
    abline(v = 50, untf = FALSE, lty = 2)
    abline(h = 50, untf = FALSE, lty = 2)

    abline(h = 150, untf = FALSE, lty = 2)
    abline(v = 150, untf = FALSE, lty = 2)
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
#' @param cex.lab Label size.
#'
#' @param cex.main Title size.
#'
#' @param cex.axis Axis text size.
#'
#' @return histogram
#'
#' @export
#'
#' @examples
#' # Example 1
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390CoverageTxt = extractCoverageFromTxt(refFile, altFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageTxt$altCount, PG0390CoverageTxt$refCount )
#' histWSAF(obsWSAF)
#' myhist = histWSAF(obsWSAF, FALSE)
#'
#' # Example 2
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390CoverageVcf = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageVcf$altCount, PG0390CoverageVcf$refCount )
#' histWSAF(obsWSAF)
#' myhist = histWSAF(obsWSAF, FALSE)
#'
histWSAF <- function ( obsWSAF, exclusive = TRUE,
                title ="Histogram 0<WSAF<1",
                cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    tmpWSAFIndex <- 1:length(obsWSAF)
    if ( exclusive ){
        tmpWSAFIndex <- which( ( (obsWSAF < 1) * (obsWSAF > 0) ) == 1)
    }
    return (hist(obsWSAF[tmpWSAFIndex], main = title,
        breaks = seq(0, 1, by = 0.1), xlab = "WSAF", col = "gray",
        cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis))
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
#' @param cex.lab Label size.
#'
#' @param cex.main Title size.
#'
#' @param cex.axis Axis text size.
#'
#' @export
#'
#' @examples
#' # Example 1
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390CoverageTxt = extractCoverageFromTxt(refFile, altFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageTxt$altCount, PG0390CoverageTxt$refCount )
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' plaf = extractPLAF(plafFile)
#' plotWSAFvsPLAF(plaf, obsWSAF)
#'
#' # Example 2
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390CoverageVcf = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageVcf$altCount, PG0390CoverageVcf$refCount )
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' plaf = extractPLAF(plafFile)
#' plotWSAFvsPLAF(plaf, obsWSAF)
#'
plotWSAFvsPLAF <- function ( plaf, obsWSAF, expWSAF = c(),
                      title = "WSAF vs PLAF",
                      cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    plot ( plaf, obsWSAF, cex = 0.5, xlim = c(0, 1), ylim = c(0, 1),
        col = "red", main = title, xlab = "PLAF", ylab = "WSAF",
        cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    if ( length(expWSAF) > 0 ){
        points ( plaf, expWSAF, cex = 0.5, col = "blue")
    }
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
#' @param cex.lab Label size.
#'
#' @param cex.main Title size.
#'
#' @param cex.axis Axis text size.
#'
#' @export
#'
#' @examples
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390CoverageVcf = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageVcf$altCount, PG0390CoverageVcf$refCount )
#' plafFile = system.file("extdata", "labStrains.test.PLAF.txt", package = "DEploid")
#' PG0390CoverageVcf.deconv = dEploid(paste("-vcf", vcfFile, "-plaf", plafFile, "-noPanel"))
#' prop = PG0390CoverageVcf.deconv$Proportions[dim(PG0390CoverageVcf.deconv$Proportions)[1],]
#' expWSAF = t(PG0390CoverageVcf.deconv$Haps) %*% prop
#' plotObsExpWSAF(obsWSAF, expWSAF)
#'
plotObsExpWSAF <- function (obsWSAF, expWSAF,
                      title = "WSAF(observed vs expected)",
                      cex.lab = 1, cex.main = 1, cex.axis = 1 ){
    plot(obsWSAF, expWSAF, pch = 19, col = "blue",
        xlab = "Observed WSAF (ALT/(ALT+REF))", ylab = "Expected WSAF (h%*%p)",
        main = title, xlim = c(-0.05, 1.05), cex = 0.5, ylim = c(-0.05, 1.05),
        cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis)
    abline(0, 1, lty = "dotted");
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
#' @seealso \code{\link{histWSAF}} for histogram.
#'
#' @export
#'
#' @examples
#' # Example 1
#' refFile = system.file("extdata", "PG0390-C.test.ref", package = "DEploid")
#' altFile = system.file("extdata", "PG0390-C.test.alt", package = "DEploid")
#' PG0390CoverageTxt = extractCoverageFromTxt(refFile, altFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageTxt$altCount, PG0390CoverageTxt$refCount )
#'
#' # Example 2
#' vcfFile = system.file("extdata", "PG0390-C.test.vcf.gz", package = "DEploid")
#' PG0390CoverageVcf = extractCoverageFromVcf(vcfFile)
#' obsWSAF = computeObsWSAF( PG0390CoverageVcf$altCount, PG0390CoverageVcf$refCount )
#'
computeObsWSAF <- function (alt, ref) {
    return ( alt / (ref + alt + 0.00000001) )
}


#' @title Painting haplotype according the reference panel
#'
#' @description Plot the posterior probabilities of a haplotype given the refernece panel.
#'
#' @param posteriorProbabilities Posterior probabilities matrix with the size of number of loci by the number of reference strain.
#'
#' @param title Figure title.
#'
#' @param labelScaling Scaling parameter for plotting.
#'
#' @export
#'
haplotypePainter <- function (posteriorProbabilities, title = "", labelScaling,
                        numberOfInbreeding = 0){
    rainbowColorBin <- 16
    rainbowColors = rainbow(rainbowColorBin)
    if ( numberOfInbreeding > 0 ){
        panelSize <- dim(posteriorProbabilities)[2]-numberOfInbreeding
        rainbowColors <- c(rep("#46a8e1", panelSize),
                           rep("#f34747", numberOfInbreeding))
    }
    barplot(t(posteriorProbabilities), beside = F, border = NA,
        col = rainbowColors, space = 0, xlab = "SNP index",
        ylab = "", main = title, cex.axis = labelScaling / 5,
        cex.lab = labelScaling / 6, cex.main = labelScaling / 5,
        xaxt = "n", yaxt = "n")
    newXaxt = round(seq(1, dim(posteriorProbabilities)[1], length.out = 6))
    axis(1, at = newXaxt, labels = as.character(newXaxt),
        cex.axis= labelScaling / 7)
    newYaxt = seq(0, 1, length.out = 3)
    axis(2, at = newYaxt, labels = as.character(newYaxt),
        cex.axis= labelScaling / 7)
}
