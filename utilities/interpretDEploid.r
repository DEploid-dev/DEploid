#!/usr/bin/env Rscript
rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid-r/.DEploid"
# DESCRIPTION:
#
# USAGE:
#    ./interpretDEploid.r -vcf FILE -plaf FILE -dEprefix STRING -o STRING
#    R --slave "--args -vcf FILE -plaf FILE -dEprefix STRING -o STRING " < utilities/interpretDEploid.r
#
# EXAMPLE:
#    ./interpretDEploid.r -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel " < utilities/interpretDEploid.r
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanelExclude -o PG0390-CNopanelExclude -exclude data/testData/labStrains.test.exclude.txt " < utilities/interpretDEploid.r

if (!exists("dEploidRootDir")){
    print("dEploidRootDir undefined, try make dEploid again!")
}

library(methods)
source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

if (myInput$helpBool){
    fun.print.help.interpret()
}

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

myExcludeInfo = fun.extract.exclude (myInput$excludeFileName, myInput$excludeBool)

if ( myInput$skip1Bool == FALSE ){
    fun.interpretDEploid.1 (myCoverageInfo, myPlafInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool)
}

fun.interpretDEploid.2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool)

fun.interpretDEploid.3 (myInput$dEploidPrefix, myInput$outPrefix, myInput$pdfBool, myInput$inbreedingBool)

if (myInput$ringBool == TRUE){
    fun.interpretDEploid.2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix, myExcludeInfo, myInput$pdfBool, myInput$ringBool)
    fun.interpretDEploid.3.ring (myInput$dEploidPrefix, myInput$outPrefix, myInput$pdfBool, myInput$inbreedingBool, myCoverageInfo, myExcludeInfo, myInput$ringDecreasingOrder, myInput$trackHeight, myInput$transformP)
}

if (myInput$ibdBool == TRUE){
    if ( myInput$skip1Bool == FALSE ){
        fun.interpretDEploid.1 (myCoverageInfo, myPlafInfo, paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myExcludeInfo, myInput$pdfBool)
    }

    fun.interpretDEploid.2 (myCoverageInfo, paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myExcludeInfo, myInput$pdfBool)

    fun.interpretDEploid.3 (paste(myInput$dEploidPrefix, ".ibd", sep=""), paste(myInput$outPrefix, ".ibd", sep=""), myInput$pdfBool, myInput$inbreedingBool)

}
