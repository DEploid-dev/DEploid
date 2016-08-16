The testing data will be working around sample `PG0390-C` from the Pf3k release 5.1 data set.
To generate the data, assume there is the file `labStrains_samples.vcf` `labStrains.wg.PLAF.txt` `labStrains.wg.panel.txt` in place, these files are obtained from Pf3k analysis. Then

```
R --slave "--args PG0390-C " < extractVcf.r
```
get `PG0390-C.wg.vcf`

```R

plaftab=read.table("labStrains.wg.PLAF.txt",header=T)
nonZeroIndex = which(plaftab$PLAF>0)
write.table (plaftab[nonZeroIndex,], file = "labStrains.eg.PLAF.txt", quote=F, sep="\t", row.names=F)

panel = read.table("labStrains.wg.panel.txt", header=T)
write.table (panel[nonZeroIndex,], file = "labStrains.eg.panel.txt", quote=F, sep="\t", row.names=F)

vcfName = "PG0390-C.wg.vcf"
skipNum = as.numeric(system(paste("cat ", vcfName, " | head -500 | grep \"##\" | wc -l"), T))
vcf  = read.table( vcfName, skip=skipNum, header=T, comment.char="", stringsAsFactors = FALSE, check.names=FALSE)
names(vcf)[1] = "#CHROM"

system ( paste("grep \"##\"", vcfName, ">", "PG0390-C.eg.vcf") )
write.table(vcf[nonZeroIndex,], file = "PG0390-C.eg.vcf", append = T, sep = "\t", quote = F, row.names = F)

set.seed(1)
testIndex = sort(unique(c(1:200,round(runif(400, 1, length(nonZeroIndex))))))
write.table (plaftab[nonZeroIndex,][testIndex,], file = "../testData/labStrains.test.PLAF.txt", quote=F, sep="\t", row.names=F)
write.table (panel[nonZeroIndex,][testIndex,], file = "../testData/labStrains.test.panel.txt", quote=F, sep="\t", row.names=F)
system ( paste("grep \"##\"", vcfName, ">", "../testData/PG0390-C.test.vcf") )
write.table(vcf[nonZeroIndex,][testIndex,], file = "../testData/PG0390-C.test.vcf", append = T, sep = "\t", quote = F, row.names = F)

```


