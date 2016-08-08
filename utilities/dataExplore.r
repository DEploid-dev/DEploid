# DESCRIPTION:
#
# USAGE:
#    R --slave "--args -vcf FILE -plaf FILE " < dataExplore.r > dataExplore.rout
#
# EXAMPLE:
#    R --slave "--args -vcf tests/testData/PG0389-C100.vcf -plaf tests/testData/plafForTesting.csv " < utilities/dataExplore.r > dataExplore.rout
./dEploid -vcf   -plaf

rm(list= ls())

checkAndIncreaseArgI <- function ( argi ){
    return (argi+1)
}

vcfFileName = ""
plafFileName = ""

args=(commandArgs(TRUE))

arg_i = 1
while ( arg_i < length(args) ){
    argv = args[arg_i]
    if ( argv == "-vcf" ){
        arg_i = checkAndIncreaseArgI ( arg_i )
        vcfFileName = args[arg_i]
    } else if ( argv == "-plaf" ){
        arg_i = checkAndIncreaseArgI ( arg_i )
        plafFileName = args[arg_i]
    } else {
        cat ("Unknow flag: ", argv, "\n")
    }

    arg_i = arg_i + 1
}

cat ("vcfFileName: ", vcfFileName, "\n")

cat ("plafFileName: ", plafFileName, "\n")


