#!/bin/bash

# This tool creates a UCSC track line given a bigWig or bedGraph input file, a track name and genome assembly.

## Define URL prefix

URL="http://data.bric.dk/galaxy/UCSC"

#parse input
while getopts hi:t:n:o: myarg
do      case "$myarg" in
        h)      echo "Usage: create UCSC track line -i <INPUT_FILE> -t <FILE_TYPE> -n <TRACK_NAME> -o <ASSEMBLY>"
                exit ;;
        i)      input_file="$OPTARG" ;;
	t)	file_type="$OPTARG" ;;
        n)      track_name="$OPTARG" ;; 
        o)      org_assembly="$OPTARG" ;;
        [?])    echo "Usage: create UCSC track line -i <INPUT_FILE> -t <FILE_TYPE> -n <TRACK_NAME> -o <ASSEMBLY>"

                exit 1 ;;
        esac
done

####################################################################
# Check if a link to that file exists in the UCSC webfolder
flag=0 

for link in `ls -1 $UCSCPATH` 
do 
    if [ `readlink -f $input_file` == `readlink -f $UCSCPATH/$link` ] ; then
        link_path=/$link # use already existing link
	flag=1 
    fi
done

# If link to input file doesn't exist, create anonymous soft-link for input file in a folder accessible from the web.

if [ $flag == 0 ]; then
    
    link_path=`mktemp -u --tmpdir=/`
    ln -s $input_file $UCSCPATH$link_path
fi

####################################################################
##Convert 4th character of file_type to uppercase (e.g. bigwig -> bigWig). UCSC browser doesn't recognize the file type otherwise.

##Split file_type string to prefix (big) and suffix (wig, bed) and uppercase the first character of the latter.
prefix=`echo ${file_type:0:3}`
suffix=`echo ${file_type:3} | sed -e "s/\b\(.\)/\u\1/g"`

file_type=`echo ${prefix}${suffix}`

####################################################################
#Contruct URL
bigURL=`echo ${URL}${link_path}`

####################################################################
## Print track line
echo "track name='$track_name' db=$org_assembly type=$file_type visibility=full alwaysZero=ON bigDataUrl=$bigURL"
