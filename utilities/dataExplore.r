# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE -o STRING " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    R --slave "--args -vcf tests/testData/PG0389-C.vcf -plaf tests/testData/labStrains_samples_PLAF.txt -o PG0389-C " < utilities/dataExplore.r

rm(list= ls())
source ("utilities/dEploidTools.r")

args = (commandArgs(TRUE))

myInput = fun.parse ( args )

myVcfInfo = fun.extract.vcf ( myInput$vcfFileName )

myPlafInfo = fun.extract.plaf ( myInput$plafFileName )

fun.dataExplore (myVcfInfo, myPlafInfo, myInput$outPrefix)

#ls()
