#!/bin/bash


#get command line parameters
ARGS=`getopt -o i: --long type:,index_path: -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-i) inputfile=$2 ; shift 2;;
		--index_path) index_path=$2 ; shift 2;;
		--type) type=$2 ; shift 2;;
		--)
			shift
			break
			;;
		*) echo "unknown parameter:" $1; helpdoc; exit 1;;
	esac
done

if [ $type == "rRNA" ]; then
    bowtie2 -q --sensitive -a -p 8 --no-mixed --reorder -x ./${index_path}/RNA45S5 --no-unal -U ${inputfile} -S ${inputfile%.trim.fq.gz}_rRNA.sam > ${inputfile%.trim.fq.gz}_rRNA.log 2>&1
    samtools view ${inputfile%.trim.fq.gz}_rRNA.sam -F256 -bS | samtools sort - > ${inputfile%.trim.fq.gz}_rRNA.bam
    samtools index ${inputfile%.trim.fq.gz}_rRNA.bam
fi

if [ $type == "genome" ]; then
    hisat2 -k 1 -t --dta --max-seeds 20 --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --no-softclip --known-splicesite-infile ./${index_path}/ref_all_spsites.txt -p 10 -x ./${index_path}/hg38_all -U ${i}  -S ./${inputfile%.trim.fq.gz}.sam > ./${inputfile%.trim.fq.gz}.log 2>&1
    samtools view ./${inputfile%.trim.fq.gz}.sam -F256 -bS | samtools sort - > ./${inputfile%.trim.fq.gz}.bam
    samtools index ./${inputfile%.trim.fq.gz}.bam
fi

