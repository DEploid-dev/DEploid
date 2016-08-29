# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -dEprefix STRING -o STRING " < utilities/interpretDEploid.r
#
# EXAMPLE:
#    R --slave "--args -vcf data/testData/PG0389-C.vcf -plaf data/testData/labStrains_samples_PLAF.txt -dEprefix PG0389-Cpanel -o PG0389-Cpanel " < utilities/interpretDEploid.r

rm(list= ls())
source ("utilities/dEploidTools.r")

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myCoverageInfo = fun.extract.coverage ( myInput )

myPlafInfo = fun.extract.plaf ( myInput$plafFileName )

print(myInput$dEploidPrefix)

fun.interpretDEploid.1 (myCoverageInfo, myPlafInfo, myInput$dEploidPrefix, myInput$outPrefix)

fun.interpretDEploid.2 (myCoverageInfo, myInput$dEploidPrefix, myInput$outPrefix)

fun.interpretDEploid.3 (myInput$dEploidPrefix, myInput$outPrefix)
