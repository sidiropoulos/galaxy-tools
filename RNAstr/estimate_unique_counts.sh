#!/bin/bash

############################################################
#Create table of fragments containing coordinates[RNA molecule, start,end], number of mapped reads, number of unique barcodes and 3'most nucleotide of cDNA (on that was ligated to)
############################################################
#Remove untemplated nucleotides and create positions_temp file 

#defaults
output_dir="output_dir"
priming_pos=-1

#parse input
while getopts hf:b:p:o: myarg
do	case "$myarg" in
	h)	echo "Usage: estimate_unique_counts.sh -f <bam_file> -b <barcodes_file> -o <output_dir>"
		exit ;;
	f)	bamfile="$OPTARG" ;; #required
	b)	barcodes="$OPTARG" ;; #required
	p)	priming_pos="$OPTARG" ;; #optional
	o)	output_dir="$OPTARG" ;;	#optional
	[?])	echo "Usage: estimate_unique_counts.sh -f <myfile.bam> -b <barcodes.txt>"
		exit 1 ;;
	esac
done

mkdir $output_dir

samtools view $bamfile |  awk 'BEGIN{OFS="\t"}{if(substr($0,1,1)!="@"){print}}' - | awk -v out="${output_dir}/trimming_stats.txt" 'BEGIN{OFS="\t";counter[0]=0;counter[1]=0;counter[2]=0;counter[3]=0}
function abs(value){return(value<0?-value:value)}
($2==99 && /[\s\t]MD:Z:/ && !/MD:Z:([012][ACGT])/)  {print($1, $3, $4+0, $4+abs($9)-1);counter[0]++;next};
($2==99 && /[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/) {print($1, $3, $4+1, $4+abs($9)-1);counter[1]++;next};
($2==99 && (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/ || (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/))) {print($1, $3, $4+2, $4+abs($9)-1);counter[2]++;next};
($2==99 && ((/[\s\t]MD:Z:1[ACGT]0[ACGT]/)||(/MD:Z:0[ACGT]1[ACGT]/)||(/MD:Z:0[ACGT]0[ACGT]0[ACGT]/)||(/MD:Z:2[ACGT]/))) {print($1, $3, $4+3, $4+abs($9)-1);counter[3]++;next}
END{print("No trimming:",counter[0],"1 nt trimmed:", counter[1],"2 nt trimmed:", counter[2],"3 nt trimmed:",counter[3]) > out}' | sort -k1,1 | gzip > positions_temp_sorted.gz &

# Computing barcode length (Use the first line and compute the string length of the second column

TMP=`head -1 $barcodes | awk '{print $2}'`
BAR_LEN=`echo ${#TMP}`

# Remove "@" from barcodes and sort them
sed 's/^.//' $barcodes | sort -k1,1 | gzip > barcodes_temp_sorted.gz &

wait 

# Merge poistions and barcodes
join -1 1 <(zcat positions_temp_sorted.gz) <(zcat barcodes_temp_sorted.gz) | cut -f 2,3,4,5 -d " " | awk '{if($4 !~ /N/){print}}' | awk -v bar_len="${BAR_LEN}" '{if(length($4)==bar_len){print}}' | gzip > merged_temp.gz

#### If priming flag is set....
if [ $priming_pos != -1 ]; then
    zcat merged_temp.gz | awk -v pos="${priming_pos}" '{print $1, $2, pos, $4}' - > merged_temp2
    cat merged_temp2 | gzip > merged_temp.gz
    rm merged_temp2
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
for(end_position in barcode[RNA][start_position]){print RNA,start_position,end_position,barcode[RNA][start_position][end_position]}}}}' $output_dir/summary.txt > $output_dir/unique_barcodes.txt &

#read_counts.gz colums: RNA_ID, Start, End, sequenced_count

zcat merged_temp.gz | awk '{barcode[$1][$2][$3]++}END{
for(RNA in barcode){
for(start_position in barcode[RNA]){
for(end_position in barcode[RNA][start_position]){print RNA,start_position,end_position,barcode[RNA][start_position][end_position]}}}}' > $output_dir/read_counts.txt &

wait

##Remove temporary files - didn't do it, in case debugging needed.

rm positions_temp_sorted.gz
rm barcodes_temp_sorted.gz
#rm merged_temp.gz

##End of remove temp files
