#!/bin/bash

####################################################################################################
#Copyright (C) 2014 Lukasz Kielpinski, Nikos Sidiropoulos

#This program is free software: you can redistribute it and/or modify it under the terms of the
#GNU General Public License as published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
#even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details (http://www.gnu.org/licenses/).
####################################################################################################

#Preprocessing workflow - Read debarcoding and trimming.

function print_help {
cat <<End-of-message
RNA probing data preprocessing.
Read trimming and debarcoding.
-------------------------------------
Input arguments:
-h: Help
-1: Read1 (FASTQ) - Required
-2: Read2 (FASTQ) - Optional
-b: Barcode signature
-t: Trimming length
-o: Output folder (default: "output_dir")
-------------------------------------
Usage : preprocessing.sh -f <READ1> -r <READ2> -b <BARCODE_SEQ> -t <TRIM_LENGTH> -o <output_dir>
End-of-message
exit
}

#defaults
output_dir="./output_dir"
read2=""
trim_length=0

#parse input
while getopts h1:2:b:t:o: myarg
do      case "$myarg" in
        h)      print_help
                exit ;;
        1)      read1="$OPTARG" ;; #required
        2)      read2="$OPTARG" ;; #optional
        b)      barcode="$OPTARG" ;; #optional
        t)      trim_length="$OPTARG" ;; #optional
        o)      output_dir="$OPTARG" ;; #optional
        [?])    echo "ERROR: Unknown parameter"
                print_help
                exit 1 ;;
        esac
done

if [[ -z $read1 ]]; then
  echo "ERROR: Read1 file required!"
  print_help
  exit 1
fi

# Create the output folder
mkdir -p $output_dir

if [ -z "$barcode" ]; then
    barcode=''
fi
BAR_LENGTH=`eval echo ${#barcode}`

###################################################################################################
# Convert the Barcode signature to a regular exression.

# Example:
#
# Barcode in the oligo used: NWTRYSNNNN
# Which means that each proper read must begin with: NNNN(C|G)(A|G)(C|T)A(A|T)N
# As regex:
# ^[ACGT][ACGT][ACGT][ACGT][CG][AG][CT][A][AT][ACGT]
# Reverse complement the barcode sequence and create a regular expression.

REV_COMPLEMENT=$(perl -0777ne's/\n //g; tr/ATGCatgcNnYyRrKkMmBbVvDdHh/TACGtacgNnRrYyMmKkVvBbHhDd/; print scalar reverse $_;' <(echo $barcode))

REG_EXP=$(sed 's/A/[A]/g;s/C/[C]/g;s/G/[G]/g;s/T/[T]/g;s/R/[AG]/g;s/Y/[CT]/g;s/S/[GC]/g;s/W/[AT]/g;s/K/[GT]/g;s/M/[AC]/g;s/B/[CGT]/g;s/D/[AGT]/g;s/H/[ACT]/g;s/V/[ACG]/g;s/N/[AGCTN]/g;' <<< $REV_COMPLEMENT)

###################################################################################################
#Remove reads that do not start with the signature (Read1)

#First awk removes them, second awk removes corresponding quality strings.
#Remove last N nt - may be derived from the random primers.

awk '{if(NR%4==2){if(/^'"$REG_EXP"'/){print}else{print ""}}else{print}}' $read1 |

awk 'BEGIN{trim_flag=0; trimming_stats=0; all_processed=0}
{
    if(NR%4==1){print; all_processed++}
    if(NR%4==2){if(length($1)==0){trim_flag=1;trimming_stats++}else{trim_flag=0};print}
    if(NR%4==3){print}
    if(NR%4==0){if(trim_flag==1){print ""}else{print $0}}
}END{print(trimming_stats, all_processed) > "trimming_stats.error"}' |

awk -v len="${trim_length}" '{if(NR%2==0){print(substr($1,0,length($1)-len))}else{print}}' |

awk -v len="${BAR_LENGTH}" '{if(NR%2==0 && length($1)<20+len){printf("\n")}else{print}}' | gzip  > R1.fastq.gz &

wait

mv trimming_stats.error $output_dir/trimming_stats.txt

###################################################################################################
#Trim primers and possible random barcode

if [ -z "$read2" ]; then
    #single-end

    #Extract the barcode sequence from the first read:
    zcat R1.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%2==0 && length($1)<20+len){printf("\n")}else{if(NR%2==0){print(substr($0,len+1,length($0)))}else{print($0)}}}' | awk '{print($1)}' > $output_dir/read1.fastq &

    zcat R1.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%4==1){print($1)}else{if(NR%4==2){print(substr($0,0,len))}}}' | paste - - > $output_dir/barcodes.txt &

    wait

    #Remove temp files
    rm R1.fastq.gz

else
    #paired-end

    #Trim primers (Read2)
    awk -v len1="${trim_length}" -v len2="${BAR_LENGTH}" '{if(NR%2==0){print(substr($0,len1+1,(length($0)-len1-len2)))}else{print($0)}}' $read2 |
    awk '{if(NR%2==0 && length($1)<20){printf("\n")}else{print}}' | gzip > R2.fastq.gz &

    wait

    #Remove empty reads - remove each pair from for which at least one read of the pair got removed (they are problematic when mapping)

    #First define which lines to keep from both fastq files (k for keep, d for discard in the lines_to_keep file)
    paste <(zcat R1.fastq.gz) <(zcat R2.fastq.gz) | awk 'BEGIN{OFS="\n"}{if(NR%4==2 && NF==2){print("k","k","k","k")}else{if(NR%4==2 && NF<2){print("d","d","d","d")}}}' > lines_to_keep

    paste lines_to_keep <(zcat R1.fastq.gz) | awk '{if($1=="k")print($2,$3)}' | gzip > R1_readsANDbarcodes.fastq.gz &

    paste lines_to_keep <(zcat R2.fastq.gz) | awk '{if($1=="k")print($2,$3)}' | awk '{print($1)}' > $output_dir/read2.fastq &

    wait

    ########################################################################
    #Extract the barcode sequence from the first read:
    zcat R1_readsANDbarcodes.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%2==0 && length($1)<20+len){printf("\n")}else{if(NR%2==0){print(substr($0,len+1,length($0)))}else{print($0)}}}' | awk '{print($1)}' > $output_dir/read1.fastq &

    zcat R1_readsANDbarcodes.fastq.gz | awk -v len="${BAR_LENGTH}" '{if(NR%4==1){print($1)}else{if(NR%4==2){print(substr($0,0,len))}}}' | paste - - > $output_dir/barcodes.txt &

    wait

    #Remove temp files
    rm R1_readsANDbarcodes.fastq.gz R1.fastq.gz R2.fastq.gz lines_to_keep

fi
