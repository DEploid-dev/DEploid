# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    R --slave "--args -vcf tests/testData/PG0389-C.vcf -plaf tests/testData/labStrains_samples_PLAF.txt -o PG0389-C " < utilities/dataExplore.r

rm(list= ls())
source (


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


fun.makingPlots <- function (vcfInfo, plafInfo) {
    PLAF = plafInfo$PLAF
    ref = vcfInfo$refCount
    alt = vcfInfo$altCount

    png ( paste ( myInput$outPrefix, "altVsRefAndWSAFvsPLAF.png", sep = "" ), width = 1800, height = 600)
    par( mfrow = c(1,3) )
    tmp.range = 1.1*min(max(alt), max(ref))

    plot(ref, alt, xlim=c(0, tmp.range), ylim=c(0,tmp.range), cex = 0.5, main = paste("Alt vs Ref"), xlab = "REF", ylab = "ALT")
    abline(v =50, untf = FALSE, lty = 2)
    abline(h =50, untf = FALSE, lty = 2)

    abline(h =150, untf = FALSE, lty = 2)
    abline(v =150, untf = FALSE, lty = 2)

    WSAF = alt / (ref + alt + 0.00000001)
    tmpWSAF_index = which(((WSAF<1) * (WSAF>0) ) == 1)

    hist(WSAF[tmpWSAF_index], main="Histogram 0<WSAF<1", xlab = "WSAF")

    print(length(alt))
    print(length(WSAF))
        print(length(PLAF))
    plot ( PLAF, WSAF, cex = 0.5, main = paste("PLAF vs WSAF"), xlab = "PLAF", ylab = "WSAF" )

    dev.off()
}

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myVcfInfo = fun.extract.vcf ( myInput$vcfFileName )

myPlafInfo = fun.extract.plaf ( myInput$plafFileName )

fun.makingPlots (myVcfInfo, myPlafInfo)

ls()
