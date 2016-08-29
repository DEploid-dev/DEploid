rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid"
# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -dEprefix STRING -o STRING " < utilities/interpretDEploid.r
#
# EXAMPLE:
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel " < utilities/interpretDEploid.r
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf -plaf data/testData/labStrains.test.PLAF.txt -dEprefix PG0390-CNopanel -o PG0390-CNopanel -exclude data/testData/labStrains.test.exclude.txt " < utilities/interpretDEploid.r

if (!exists("dEploidRootDir")){
    print("dEploidRootDir undefined, try make dEploid again!")
}

source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = fun.extract.plaf ( myInput$plafFileName )

fun.interpretDEploid.1 (myCoverageInfo, myPlafInfo, myInput$dEploidPrefix, myInput$outPrefix)

fun.interpretDEploid.2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix)

fun.interpretDEploid.3 (myInput$dEploidPrefix, myInput$outPrefix)
