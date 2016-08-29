# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    R --slave "--args -vcf data/testData/PG0390-C.test.vcf.gz -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r
#    R --slave "--args -ref data/testData/PG0390-C.test.ref -alt data/testData/PG0390-C.test.alt -plaf data/testData/labStrains.test.PLAF.txt -o PG0390-C " < utilities/dataExplore.r

rm(list= ls())
source ("utilities/dEploidTools.r")

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = fun.extract.plaf ( myInput$plafFileName )

fun.dataExplore (myCoverageInfo, myPlafInfo, myInput$outPrefix)

#ls()
