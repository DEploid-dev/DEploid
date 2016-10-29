rm(list=ls()); dEploidRootDir="/home/joezhu/DEploid-r/.DEploid"
# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r
#    R --slave "--args -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r

if (!exists("dEploidRootDir")){
    print("dEploidRootDir undefined, try make dEploid again!")
}

source(paste(dEploidRootDir, "/utilities/dEploidTools.r", sep=""))
source(paste(dEploidRootDir, "/utilities/DEploidR.R", sep=""))

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = extractPLAF( myInput$plafFileName )

fun.dataExplore (myCoverageInfo, myPlafInfo, myInput$outPrefix)

#ls()
