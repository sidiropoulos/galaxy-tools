#!/bin/bash

#Preprocessing workflow - trimming adapters, barcodes etc

# Barcode in the oligo used: NWTRYSNNNN

# Which means that each proper read must begin with: NNNN(C|G)(A|G)(C|T)A(A|T)N

# As regex:

# ^[ACGT][ACGT][ACGT][ACGT][CG][AG][CT][A][AT][ACGT]


READ1=${1}
READ2=${2}
BARCODE=$3
ADAPTER1=$4
ADAPTER2=$5
CUTOFF=$6
CUTADAPT_LOG=${7}
TRIM_LENGTH=$8
OVERLAP_LENGTH=$9
OUT_R1=${10}
OUT_R2=${11}
OUT_BARCODES=${12}

BAR_LENGTH=`eval echo ${#BARCODE}`

#########################################################################################

# Reverse complement the barcode sequence and create a regular expression.

REV_COMPLEMENT=$(perl -0777ne's/\n //g; tr/ATGCatgcNnYyRrKkMmBbVvDdHh/TACGtacgNnRrYyMmKkVvBbHhDd/; print scalar reverse $_;' <(echo $BARCODE))

REG_EXP=$(sed 's/A/[A]/g;s/C/[C]/g;s/G/[G]/g;s/T/[T]/g;s/R/[AG]/g;s/Y/[CT]/g;s/S/[GC]/g;s/W/[AT]/g;s/K/[GT]/g;s/M/[AC]/g;s/B/[CGT]/g;s/D/[AGT]/g;s/H/[ACT]/g;s/V/[ACG]/g;s/N/[AGCT]/g;' <<< $REV_COMPLEMENT)

########################################################################
#1. Remove all the reads that do not start with the signature 

# (first awk removes them, second awk removes corresponding quality strings) followed by cutadapt. 

# After cutadapt, remove last 15 nt - may be derived from the random primer [should be optional, as the user may: 1) use different random primer length, 2) have short reasd and would lose too much info] and remove all the reads that are shorter than 30 nt 

# (10 nt barcode + 20 nt for mapping)

echo -e '------------------------Read1------------------------\n' > $CUTADAPT_LOG

awk '{if(NR%4==2){if(/^'"$REG_EXP"'/){print}else{print ""}}else{print}}' $READ1 | 

awk 'BEGIN{trim_flag=0; trimming_stats=0; all_processed=0}
{
    if(NR%4==1){print; all_processed++}
    if(NR%4==2){if(length($1)==0){trim_flag=1;trimming_stats++}else{trim_flag=0};print}
    if(NR%4==3){print}
    if(NR%4==0){if(trim_flag==1){print ""}else{print $0}}
}END{print(trimming_stats, all_processed) > "trimming_stats.error"}' | 

cutadapt -a $ADAPTER1 -q $CUTOFF --format=fastq -O $OVERLAP_LENGTH - 2>>$CUTADAPT_LOG |

awk -v len="${TRIM_LENGTH}" '{if(NR%2==0){print(substr($1,0,length($1)-len))}else{print}}' | 

awk -v len="${BAR_LENGTH}" '{if(NR%2==0 && length($1)<20+len){printf("\n")}else{print}}' | gzip  > R1.fastq.gz &

wait

########################################################################
#2. Trim the adapter, primer and possible random barcode from the second read

echo -e '------------------------Read2------------------------\n' >> $CUTADAPT_LOG

cutadapt -a $ADAPTER2 -q $CUTOFF --format=fastq -O $OVERLAP_LENGTH $READ2 2>>$CUTADAPT_LOG | 

awk -v len1="${TRIM_LENGTH}" -v len2="${BAR_LENGTH}" '{if(NR%2==0){print(substr($0,len1+1,(length($0)-len1-len2)))}else{print($0)}}' - |

awk '{if(NR%2==0 && length($1)<20){printf("\n")}else{print}}' | gzip > R2.fastq.gz &

wait

########################################################################
#3. Remove empty reads - remove each pair from for which at least one read of the pair got removed (they are problematic for tophat mapping)

#first define which lines to keep from both fastq files (k for keep, d for discard in the lines_to_keep file)

paste <(zcat R1.fastq.gz) <(zcat R2.fastq.gz) | awk 'BEGIN{OFS="\n"}{if(NR%4==2 && NF==2){print("k","k","k","k")}else{if(NR%4==2 && NF<2){print("d","d","d","d")}}}' > lines_to_keep

paste lines_to_keep <(zcat R1.fastq.gz) | awk '{if($1=="k")print($2,$3)}' | gzip > R1_readsANDbarcodes.fastq.gz &

paste lines_to_keep <(zcat R2.fastq.gz) | awk '{if($1=="k")print($2,$3)}' > $OUT_R2 &

wait

########################################################################
#4. Extract the barcode sequence from the first read:
zcat R1_readsANDbarcodes.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%2==0 && length($1)<20+len){printf("\n")}else{if(NR%2==0){print(substr($0,len+1,length($0)))}else{print($0)}}}' > $OUT_R1 &

zcat R1_readsANDbarcodes.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%4==1){print($1)}else{if(NR%4==2){print(substr($0,0,len))}}}' | paste - - > ${OUT_BARCODES} &

wait 

########################################################################
#6. Remove temporary fastq files"

rm R1.fastq.gz R2.fastq.gz R1_readsANDbarcodes.fastq.gz
rm lines_to_keep

########################################################################
#7. problem! Spaces added at the end of the strings in fastq files (AWK induced). I will remove them:

mv $OUT_R1 R1.temp.fastq
mv $OUT_R2 R2.temp.fastq

awk '{print($1)}' R1.temp.fastq > $OUT_R1 &
awk '{print($1)}' R2.temp.fastq > $OUT_R2 &

wait
rm R?.temp.fastq
