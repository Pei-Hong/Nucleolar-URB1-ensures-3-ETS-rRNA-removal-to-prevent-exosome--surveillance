#!/bin/bash

file=$1

# Get barcode files and remove barcodes from raw data
zcat $file | paste - - - - | awk '{print $1"\t"substr($3,1,9)}' > ${file%.fq.gz}.barcode
zcat $file | awk 'BEGIN{n=0}{n=n+1; if(n==2){$1=substr($1,10,(length($1)-9)); print $1} else if(n==4){$1=substr($1,10,(length($1)-9)); print $1; n=0} else {print $0}}' | gzip  > ${file%.fq.gz}.rmbar.fq.gz 

# Remove adaptors
java -jar trimmomatic-0.38.jar SE -threads 10 -phred33 ${file%.fq.gz}.rmbar.fq.gz ./${file%.fq.gz}.trim.fq.gz ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq3_2-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

