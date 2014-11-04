#!/bin/bash

############################################################
#Create table of fragments containing coordinates[RNA molecule, start,end], number of mapped reads, number of unique barcodes and 3'most nucleotide of cDNA (on that was ligated to)
############################################################
#Remove untemplated nucleotides and create positions_temp file

function print_help {
cat <<End-of-message
Summarize Unique Barcodes.
Counts the number of unique random barcodes and reads associated with each sequenced fragment.
-------------------------------------
Input arguments:
-h: Help
-f: Aligned reads in BAM format.
-b: Barcode file (optional).
-p: Set priming position to a fixed value.
-t: Trim untemplated nucleotides.
-k: Produce k2n file. Warning: Can be sloooow!
-o: Output folder (default: "output_dir")
-------------------------------------
Usage : summarize_unique_barcodes.sh -f <BAM_file> -b <BARCODES> -p <PRIMING_POSITION> -t -k
End-of-message
exit
}


#defaults
output_dir="output_dir"
trim_flag="False"

#parse input
while getopts hf:b::p:o:ktr: myarg
do  case "$myarg" in
    h)  print_help
        exit ;;
    f)  bamfile="$OPTARG" ;; #required
    b)  barcodes="$OPTARG" ;; #optional
    t)  trim_flag="True" ;;
    k)  k2n="True" ;;
    p)  priming_pos="$OPTARG" ;;
    o)  output_dir="$OPTARG" ;;
    r)  R_SCRIPT_PATH="$OPTARG" ;; #required
    [?])  echo "ERROR: Unknown parameter"
        print_help
        exit 1 ;;
    esac
done

###### Sanity checks ######
if [ -z $bamfile ]; then
    echo "Error: Aligned reads file is missing!"
    print_help
    exit 1
fi

if [ "$barcodes" == "None" ] && [ "$k2n" == "True" ]; then
    echo "Error: k2n file cannot be produced without a barcode file!"
    exit 1
fi


mkdir -p $output_dir

#Check if bamfile contains single or paired-end data
samtools view -f 0x1 $bamfile | head -n 1 > paired

if [ -s paired ]; then
    #paired-end
    samtools view $bamfile | awk 'BEGIN{OFS="\t"}{if(substr($0,1,1)!="@"){print}}' - | awk -v out="${output_dir}/trimming_stats.txt" -v flag="${trim_flag}" 'BEGIN{OFS="\t";counter[0]=0;counter[1]=0;counter[2]=0;counter[3]=0;counter[4]=0}
    function abs(value){return(value<0?-value:value)}
    function return_offset(local_offset){print($1, $3, $4+local_offset, $4+abs($9)-1);counter[local_offset]++}
    ($2 != 99) {next};
    (flag == "False") {return_offset(0);next};
    (/[\s\t]MD:Z:/ && !/MD:Z:([012][ACGT])/)  {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/ && substr($10,1,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/) {return_offset(1);next};

    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/) {return_offset(2);next};

    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,2)=="NN") {return_offset(0);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,1)=="N") {return_offset(2);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};

    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,2)=="NN") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(3);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,3,1)=="N") {return_offset(2);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/) {return_offset(3);next};

    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N" && substr($10,1,1)=="N") {return_offset(0);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,1,1)=="N") {return_offset(3);next};
    (/MD:Z:0[ACGT]1[ACGT]/) {return_offset(3);next};

    (/MD:Z:2[ACGT]/ && substr($10,3,1)=="N")  {return_offset(0);next};
    (/MD:Z:2[ACGT]/)  {return_offset(3);next};

    (substr($10,3,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(3);next};
    (substr($10,2,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};
    (substr($10,1,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(1);next};

    {return_offset(0);counter[4]++}
    END{print("No trimming:",counter[0],", out of which not recognized MD field for:",counter[4],"; 1 nt trimmed:", counter[1],"; 2 nt trimmed:", counter[2],"; 3 nt trimmed:",counter[3]) > out}' | sort -S1G -k1,1 | gzip > positions_temp_sorted.gz
else
    #single-end
    samtools view $bamfile | awk 'BEGIN{OFS="\t"}{if(substr($0,1,1)!="@"){print}}' - | awk -v out="${output_dir}/trimming_stats.txt" -v flag="${trim_flag}" 'BEGIN{OFS="\t";counter[0]=0;counter[1]=0;counter[2]=0;counter[3]=0;counter[4]=0}
    function abs(value){return(value<0?-value:value)}
    function return_offset(local_offset){print($1, $3, $4+local_offset, $4+abs($9)-1);counter[local_offset]++}
    ($2 != 0) {next};
    (flag == "False") {return_offset(0);next};
    (/[\s\t]MD:Z:/ && !/MD:Z:([012][ACGT])/)  {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/ && substr($10,1,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/) {return_offset(1);next};

    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/) {return_offset(2);next};

    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,2)=="NN") {return_offset(0);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,1)=="N") {return_offset(2);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};

    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,2)=="NN") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(3);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,3,1)=="N") {return_offset(2);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/) {return_offset(3);next};

    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N" && substr($10,1,1)=="N") {return_offset(0);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,1,1)=="N") {return_offset(3);next};
    (/MD:Z:0[ACGT]1[ACGT]/) {return_offset(3);next};

    (/MD:Z:2[ACGT]/ && substr($10,3,1)=="N")  {return_offset(0);next};
    (/MD:Z:2[ACGT]/)  {return_offset(3);next};

    (substr($10,3,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(3);next};
    (substr($10,2,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};
    (substr($10,1,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(1);next};

    {return_offset(0);counter[4]++}
    END{print("No trimming:",counter[0],", out of which not recognized MD field for:",counter[4],"; 1 nt trimmed:", counter[1],"; 2 nt trimmed:", counter[2],"; 3 nt trimmed:",counter[3]) > out}' | sort -S1G -k1,1 | gzip > positions_temp_sorted.gz

fi

# Computing barcode length (Use the first line and compute the string length of the second column
if [ "$barcodes" != "None" ]; then

    TMP=`head -1 $barcodes | awk '{print $2}'`
    BAR_LEN=`echo ${#TMP}`

    # Remove "@" from barcodes and sort them
    sed 's/^.//' $barcodes | sort -k1,1 -S1G | gzip > barcodes_temp_sorted.gz

    # Merge poistions and barcodes
    join -1 1 <(zcat positions_temp_sorted.gz) <(zcat barcodes_temp_sorted.gz) | cut -f 2,3,4,5 -d " " | awk '{if($4 !~ /N/){print}}' | awk -v bar_len="${BAR_LEN}" '{if(length($4)==bar_len){print}}' | gzip > merged_temp.gz

    rm barcodes_temp_sorted.gz

else
    zcat positions_temp_sorted.gz | cut -f 2,3,4 | gzip > merged_temp.gz
fi

# If the experiment is single-end
if [ ! -s paired ]; then
    zcat merged_temp.gz | awk '{print $1, $2, "NA", $4, $5}' - > merged_temp2
    cat merged_temp2 | gzip > merged_temp.gz
fi

#### If priming flag is set....
if [ ! -z $priming_pos ]; then
    zcat merged_temp.gz | awk -v pos="${priming_pos}" '{print $1, $2, pos, $4, $5}' - > merged_temp2
    cat merged_temp2 | gzip > merged_temp.gz
fi

#File summary.gz columns: RNA_ID, Start, End, barcode sequence, sequenced_count[=number of sequenced fragments fulfilling previous requiremnts]

zcat merged_temp.gz | awk '{barcode[$1][$2][$3][$4]++}END{
for(RNA in barcode){
for(start_position in barcode[RNA]){
for(end_position in barcode[RNA][start_position]){
for(barseq in barcode[RNA][start_position][end_position]){print RNA,start_position,end_position,barseq,barcode[RNA][start_position][end_position][barseq]}}}}}' > $output_dir/summary.txt

#File unique_barcodes columns: RNA_ID, Start, End, number of unique barcodes observed for this fragment [PROBLEM: How to treat the different 3cdns for the same fragment? if the template was homogenous then it should be always the same]

awk '{barcode[$1][$2][$3]++}END{
for(RNA in barcode){
for(start_position in barcode[RNA]){
for(end_position in barcode[RNA][start_position]){print RNA "\t" start_position "\t" end_position "\t" barcode[RNA][start_position][end_position]}}}}' $output_dir/summary.txt > $output_dir/unique_barcodes.txt &

#read_counts.gz colums: RNA_ID, Start, End, sequenced_count

zcat merged_temp.gz | awk '{barcode[$1][$2][$3]++}END{
for(RNA in barcode){
for(start_position in barcode[RNA]){
for(end_position in barcode[RNA][start_position]){print RNA "\t" start_position "\t" end_position "\t" barcode[RNA][start_position][end_position]}}}}' > $output_dir/read_counts.txt &

wait

if [ "$barcodes" != "None" ]; then

    cut -f 4 $output_dir/unique_barcodes.txt | sort -S1G -rn > sorted_bars
    max_observed_barcodes=`head -n 1 sorted_bars`

    echo "Maximum observed Barcodes = ${max_observed_barcodes}"
    rm sorted_bars
fi

if [ "$k2n" == "True" ]; then
    bar_length=$(head -n 1 $barcodes | cut -f 2 | xargs expr length)
    Rscript $R_SCRIPT_PATH/k2n.R -m merged_temp.gz -c $output_dir/read_counts.txt -b ${max_observed_barcodes} -l $bar_length -o $output_dir/k2n.txt
fi

#Remove temp files
rm paired merged_temp.gz positions_temp_sorted.gz
