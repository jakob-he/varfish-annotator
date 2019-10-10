#!/bin/bash

zcat $1 | awk -F "\t" '{if($0~/^##.*/){next} else if($0~/^#C.*/){print "##fileformat=VCFv4.3\n##CADD\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} else {printf "chr%s\t%s\t.\t%s\t%s\t.\t.\tRawScore=%s;PHRED=%s\n",$1,$2,$3,$3,$4,$5,$6}}' > "${1%.tsv.gz}.vcf"
bgzip "${1%.tsv.gz}.vcf"
tabix -p vcf "${1%.tsv.gz}.vcf.gz"
