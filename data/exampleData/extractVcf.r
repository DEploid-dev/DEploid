# USAGE:
#    R --slave "--args STRING " < extractVcf.r
#
# EXAMPLE:
#    R --slave "--args PG0389-C " < extractVcf.r

rm(list= ls())

args = (commandArgs(TRUE))

vcfName = "labStrains_samples.vcf"

skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), T))
vcf  = read.table( vcfName, skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)

sampleName = args[1]

tmp = vcf[[sampleName]]

tmpVcf = cbind(vcf[,1:9], tmp )


wg.newVcfFilename = paste( sampleName, ".wg.vcf", sep="" )

system ( paste("grep \"##\"", vcfName, ">", wg.newVcfFilename) )
names(tmpVcf)[1] = "#CHROM"
names(tmpVcf)[10] = sampleName



write.table(tmpVcf, file = wg.newVcfFilename, append = T, sep = "\t", quote = F, row.names = F)



warnings()
