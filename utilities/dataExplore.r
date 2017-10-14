#!/usr/bin/env Rscript
rm(list=ls()); dEploidRootDir="/users/mcvean/joezhu/DEploid"
# DESCRIPTION:
#
# USAGE:
#    ./dataExplore.r -vcf FILE -plaf FILE -o STRING
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    utilities/dataExplore.r -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r
#    R --slave "--args -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r

if (!exists("dEploidRootDir")){
    print("dEploidRootDir undefined, try make dEploid again!")
}

source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

if (myInput$helpBool){
    fun.print.help.explore()
}

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

fun.dataExplore (myCoverageInfo, myPlafInfo, myInput$outPrefix, myInput$pdfBool, myInput$filter.threshold, myInput$filter.window)
