#!/bin/bash
# $Id: bigwig2summary_stdout 34 2014-02-20 08:31:20Z jens $

#USE IN GALAXY
#This script extract summary values (mean, min, max, std or coverage) from a bigwig file for a number of equal sized bins across genomic regions given in bed file.
#If the bed file contains 6 columns (or more), column 6 is expected to contain strand information. Summary values from a negative strand will be reversed.

########################################
#bigwigSummary has 3 non-standard outputs:
#1) n/a #(no data in bin)
#2) no data #(no data in entire region)
#3) <number> is not a valid option (typically if negative coordinate)

#Default settings & input parsing. "" indicates required user input.
nbins=1
rm_header_line=0
summary_type=mean

#parse input
while getopts hef:b:o:n:t: myarg
do	case "$myarg" in
	h)	echo "Usage: bigwig2summary_stdout -f <BIGWIG_FILE> -b <BED_FILE> -o <ASSEMBLY> -n <NUMBER OF BINS> -t <SUMMARY TYPE> -e"
	        echo "Corrects for strand if bed-file has 6 columns or more. Col 6 assumed to contain strand infomation"
		exit ;;
	f)	bigwig_file="$OPTARG" ;;
	b)	bed_file="$OPTARG" ;; #must be tab separated without header
	o)	org_assembly="$OPTARG" ;;
	n)	nbins="$OPTARG" ;;
	t)	summary_type="$OPTARG" ;;
	e)	rm_header_line=2 ;; #flag. if -e then first line is removed
	[?])	echo "Usage: bigwig2summary_stdout -f <BIGWIG_FILE> -b <BED_FILE> -o <ASSEMBLY> -n <NUMBER OF BINS> -t <SUMMARY TYPE> -e"
	        echo "Corrects for strand if bed-file has 6 columns or more. Col 6 assumed to contain strand infomation"
		exit 1 ;;
	esac
done

###################################################
###VALIDATE INPUT
###################################################

#get chromosome sizes from bigwig file. bigwig file does not contain name of genome assembly.
org_assembly_file=`mktemp -u`
fetchChromSizes $org_assembly 2>/dev/null > $org_assembly_file
if [ $? -ne 0 ]; then 
  echo "ERROR: Organism genome assembly does not exist" 
  rm $org_assembly_file
  exit
fi

#check input bed_file. bedClip only checks first 3 columns!
if [ `bedClip -verbose=2 <(tail -n +${rm_header_line} $bed_file) $org_assembly_file /dev/null 2>&1 | wc -l` -gt 0 ]; then
  echo -e "ERROR: Input bed file is not in proper format!\nTry 'bedClip' to find lines causing error"
  echo "Make sure that bigwig and bed files are using the same genome assembly"
  exit 1
fi

#make string of "nbins" 0's to insert in regions, where no reads are found
if [ $nbins -gt 1 ]; then
  seq_string=`seq 1 $nbins`
  zero_string=`printf "0\t%.s" {$seq_string} | perl -pe "s/\t$//"` 
fi

#make sure the given summary type exists
if [ `echo $summary_type | egrep "(mean|max|min|std|coverage)" | wc -l` -ne 1 ]; then
  echo "ERROR: Summary type must be: mean, max, min, std or coverage. Default is 'mean'"
  exit 1
fi

#determine number of fields in bed_file
if [ `tail -n +${rm_header_line} $bed_file | awk '{print NF}' | uniq | wc -l` -eq 1 ]; then
  nfields_bed=`tail -n +${rm_header_line} $bed_file | awk '{print NF}' | uniq`
else
  echo "ERROR: Bed file does not have constant number of line columns"
  exit 1
fi  

if [[ $nbins -gt 1 && $nfields_bed -ge 6 ]]; then
  strand_uniq_chars=`tail -n +${rm_header_line} $bed_file | cut -f6 | sort -u | perl -pe "s/\n//"`
  if [[ $strand_uniq_chars != "+-" && $strand_uniq_chars != "-+" ]] ; then
    echo "ERROR: Column 6 in bed file must only contain '+' or '-' characters"
    exit 1
  fi
fi 

###################################################
###EXTRACT DENSITIES FROM NORMALIZED BIGWIG FILE
###################################################


#if more than 1 bin AND >= 6 fields (i.e. has strand column)
if [[ $nbins -gt 1 && $nfields_bed -ge 6 ]]; then

  #cut columns 1-3+6 | rm header if flag set
  cut -f1-3,6 $bed_file | tail -n +${rm_header_line} | while read -r line; do
  
    #read 4 fields into variables
    read -r cur_chr cur_start cur_end cur_strand <<<"$line"
    #run bigWigSummary. Combine stdout and stderr | treat exceptions and errors after 'done'
    bigWigSummary $bigwig_file $cur_chr $cur_start $cur_end $nbins -type=${summary_type} 2>&1 | perl -pe "s/no data.+$/${zero_string}/" | awk 'BEGIN{OFS="\t"}{ if("'"$cur_strand"'"=="-") { for (i=NF; i>0; i--) { printf("%s\t",$i) } printf("\n") } else { print $0 } }'
  done | perl -pe "s/n\/a/0/g" | perl -pe "s/^\t?0\t?$/${zero_string}/" | perl -pe "s/ +/\t/g"  | perl -pe "s/\t$//" |  sed '/^$/d'

#if more than 1 bin AND less than 6 fields (i.e. no strand column)
elif [[ $nbins -gt 1 && $nfields_bed -lt 6 ]]; then

  #cut columns 1-3 | rm header if flag set
  cut -f1-3 $bed_file | tail -n +${rm_header_line} | while read -r line; do
  
    #read 3 fields into variables
    read -r cur_chr cur_start cur_end <<<"$line"
    #run bigWigSummary. Combine stdout and stderr | treat exceptions and errors after 'done'
    bigWigSummary $bigwig_file $cur_chr $cur_start $cur_end $nbins -type=${summary_type} 2>&1 
  
  done | perl -pe "s/n\/a/0/g" | perl -pe "s/no data.+$/${zero_string}/" | perl -pe "s/^\t?0\t?$/${zero_string}/" | perl -pe "s/ +/\t/g" | sed '/^$/d'


#if 1 bin. Strand column irrelevant 
else

  cut -f1-3 $bed_file | tail -n +${rm_header_line} | while read -r line; do
  
    read -r cur_chr cur_start cur_end <<<"$line"
    bigWigSummary $bigwig_file $cur_chr $cur_start $cur_end 1 -type=${summary_type} 2>&1
    
  done | perl -pe "s/no data.+$/0/" | perl -pe "s/n\/a/0/g" | sed '/^$/d'
  
fi 
  


  
  
