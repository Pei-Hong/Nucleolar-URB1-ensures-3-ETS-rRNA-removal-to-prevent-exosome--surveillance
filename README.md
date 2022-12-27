Nucleolar-URB1-ensures-3-ETS-rRNA-removal-to-prevent-exosome--surveillance
=======
Editing date: 20221226\
Author: Pei-Hong Zhang\
Email: zhangpeihong2018@sibs.ac.cn

## 1. URB1 bind on 45S pre-rRNA

### 1.1 Requirements
```
Bowtie2 v2.3.5
Hisat2 v2.1.0
samtools v1.9
bedtools v2.29.0
```
### 1.2 Data pre-processing
* Using `pre-processing.sh` script can get barcode file `URB1_iCLIP_rep1.barcode` and clean data without barcode and adaptors `URB1_iCLIP_rep1.trim.fq.gz`
```
pre-processing.sh URB1_iCLIP_rep1.fq.gz
pre-processing.sh URB1_iCLIP_rep2.fq.gz
```
### 1.3 Mapping to Genome
* Pre-processed data was mapped to hg38 Genome by hisat2 and get `URB1_iCLIP_rep1.bam`
```
mapping.sh -i URB1_iCLIP_rep1.trim.fq.gz --index_path /path/to/index_path/ --type genome
mapping.sh -i URB1_iCLIP_rep2.trim.fq.gz --index_path /path/to/index_path/ --type genome
```
* Remove PCR duplicates via `rm_pcr.py` and get `URB1_iCLIP_rep1_combine.bam`
```
rm_pcr.py URB1_iCLIP_rep1.barcode URB1_iCLIP_rep1.bam
rm_pcr.py URB1_iCLIP_rep2.barcode URB1_iCLIP_rep2.bam

bamToBed -i URB1_iCLIP_rep1_combine.bam > URB1_iCLIP_rep1_combine.bed
bamToBed -i URB1_iCLIP_rep2_combine.bam > URB1_iCLIP_rep2_combine.bed
```

### 1.4 Mapping to 45S pre-rRNA
* Pre-processed data was mapped to 45S pre-rRNA sequence by bowtie2 and get `URB1_iCLIP_rep1_rRNA.bam`
```
mapping.sh -i URB1_iCLIP_rep1.trim.fq.gz --index_path /path/to/index_path/ --type rRNA
mapping.sh -i URB1_iCLIP_rep2.trim.fq.gz --index_path /path/to/index_path/ --type rRNA
```
* Remove PCR duplicates via `rm_pcr.py` and get `URB1_iCLIP_rep1_rRNA_combine.bam`
```
rm_pcr.py URB1_iCLIP_rep1.barcode URB1_iCLIP_rep1_rRNA.bam
rm_pcr.py URB1_iCLIP_rep2.barcode URB1_iCLIP_rep2_rRNA.bam
```

### 1.5 URB1 binding region
* Get enrich-ratio file `URB1_iCLIP_rep1_rRNA_combine.srt.bg`
```
samtools sort URB1_iCLIP_rep1_rRNA_combine.bam > URB1_iCLIP_rep1_rRNA_combine.srt.bam
samtools index URB1_iCLIP_rep1_rRNA_combine.srt.bam
total=`cat URB1_iCLIP_rep1_combine.bed | wc -l`
genomeCoverageBed -split -d -ibam URB1_iCLIP_rep1_rRNA_combine.srt.bam| awk -v t=${total} '{per=$3/t; print $0"\t"per}' > URB1_iCLIP_rep1_rRNA_combine.srt.bg

paste URB1_iCLIP_rep1_rRNA_combine.srt.bg URB1_iCLIP_rep2_rRNA_combine.srt.bg | cut -f 1-4,7-8 > combine2.bg
Rscript bindingProfile.r combine2.bg
```

## 2. SHAPE-Map

### 2.1 Requirements
```
cutadapt v1.18
ShapeMapper v2.1.3
RNAfold v2.4.14
RNAstructure v6.4
deltaSHAPE v1.0
```
### 2.2 pre-processing
* Remove PCR primers from `${sample}-F1R1_F3R3.fq.gz`
```
for i in *_1P.fq.gz; do cutadapt -m 20 -g GGGTTTTAAGCAGGAGGTGT -g \
 GTGTTGTTGCCATGGTAATCCTGC -g GGGTTTAGACCGTCGTGAGAC -g GCGGGCCGCCCCCCCCTCCA \
 -G GGGTTTTAAGCAGGAGGTGT -G GTGTTGTTGCCATGGTAATCCTGC -G \
 GGGTTTAGACCGTCGTGAGAC -G GCGGGCCGCCCCCCCCTCCA \
 -o ${i%_1P.fq.gz}_1P.rmPrimer.fq.gz \
 -p ${i%_1P.fq.gz}_2P.rmPrimer.fq.gz ${i}  ${i%_1P.fq.gz}_2P.fq.gz \
 > ${i%_1P.fq.gz}.rmPrimer.log 2>&1; done
```
* Remove PCR primers from `${sample}-F4R4.fq.gz`
```
for i in *_1P.fq.gz; do cutadapt -g GGCCTCGGATAGCCGGTCCC -g \
GGCCCGGCGGGCGTGCGCGT -G GGCCTCGGATAGCCGGTCCC -G GGCCCGGCGGGCGTGCGCGT \
-o ${i%_1P.fq.gz}_1P.rmPrimer.fq.gz -p ${i%_1P.fq.gz}_2P.rmPrimer.fq.gz \
${i} ${i%_1P.fq.gz}_2P.fq.gz > ${i%_1P.fq.gz}.trimprimer.log 2>&1 ; done
```
* Combine `${sample}-F1R1_F3R3.fq.gz` `${sample}-F4R4.fq.gz`
```
zcat ${sample}-F1R1_F3R3.fq.gz ${sample}-F4R4.fq.gz | gzip > ${sample}.comb.fq.gz
```
* shapemapper
```
for treat in SCR KD1 KD4; do 
    for rep in NAI-1 NAI-2 ; do 
        hapemapper --overwrite --name ${rep}-${treat} --target ../F1R1_F3R3_F4R4.fa \
        --out ${rep}-${treat} --verbose --serial  --nproc 25 --min-depth 1000 \
        --modified --R1 ${rep}-${treat}_1P.rmPrimer.comb.fq.gz --R2 ${rep}-${treat}_2P.rmPrimer.comb.fq.gz \
        --untreated --R1 DMSO-${treat}_1P.rmPrimer.comb.fq.gz --R2 DMSO-${treat}_2P.rmPrimer.comb.fq.gz \
        --denatured --R1 DC-${treat}_1P.rmPrimer.comb.fq.gz --R2 DC-${treat}_2P.rmPrimer.comb.fq.gz \
        > ${rep}-${treat}.log 2>&1
    done
done
```
* Predict RNA secondary structure
```
for i in */*_F1R1_F3R3_F4R4.shape ; do 
    shape2struc.py -s ${i} -f ../F1R1_F3R3_F4R4.fa -o ${i%.shape}
done > log 2>&1 &
```
* deltaSHAPE
```
for rep in 1 2 ; do 
    for kd in KD1 KD4; do 
        echo "deltaSHAPE.py NAI-${rep}-${kd}_F1R1_F3R3_F4R4.map NAI-${rep}-SCR_F1R1_F3R3_F4R4.map -o p5/NAI-${rep}-${kd}.txt --all --colorfill --noshow --pdf -p 5 -z 0 -s 0 "
    done
done > deltaSHAPE.sh
cat deltaSHAPE.sh | parallel -j 3
```
## 3.Distribution analysis of proteins of interest (POIs)
* Each POI expression plasmid was transfected and co-imaged with mRuby3-DKC1 and mTagBFP2-NPM1 by epifluorescence microscope (DeltaVision Elite). The 3D stacks were imported into Fiji/ImageJ and analyzed by an in-house ImageJ script `20221226_Nuc_Classifier.ijm`.

## 4. License
Copyright (C) 2022 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use.