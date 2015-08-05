# Setup R error handling to go to stderr
options( show.error.messages = FALSE, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# Declare functions

k2n_calc <- function(merged_file, unique_barcode_file, output_file, paired = TRUE)
{

    # Read in inputs:
    merged <- read.table( merged_file )
    colnames(merged) <- c("RNAid", "Start", "End", "Barcodes")
    barcode_length <- max(nchar(as.character(merged$Barcodes)))

    #Single end
    if( !paired)
	read_counts <- as.data.frame(table(subset(merged, select = - (End:Barcodes))))
    #Paired end
    else
        read_counts <- as.data.frame(table(subset(merged, select = - Barcodes)))
    read_counts <- read_counts[read_counts$Freq != 0, ]

    max_observed <- max(read.table(unique_barcode_file)[,4])

    #Check if max_observed equals max possible complexity. If yes - subtract 1
    #(otherwise 'while' loop below never ends)
    if( max_observed == 4**barcode_length ) {
        max_observed <- max_observed - 1
        barcodes_saturated <- TRUE
    } else
        barcodes_saturated <- FALSE

    # Remove top-quartile of reads:
    #Single end
    if (!paired){
	barcodes_nt <- merged[ do.call(
            paste, as.list(subset(merged, select = - (End:Barcodes)))) %in%
                do.call(paste, as.list(read_counts[(
                    read_counts$Freq ) <= quantile(read_counts$Freq, 0.75),
                    c("RNAid", "Start") ] ) ) , "Barcodes" ]
    } else {
    barcodes_nt <- merged[ do.call(
        paste, as.list(subset(merged, select = - Barcodes))) %in%
            do.call(paste, as.list(read_counts[(
                read_counts$Freq ) <= quantile(read_counts$Freq, 0.75),
                c("RNAid", "Start", "End") ] ) ) , "Barcodes" ]
    }

    # make the matrix with the nucleotide freqs per position:
    nt_counts <- matrix( nrow = 4, ncol = barcode_length )
    for( h in 1:ncol( nt_counts ) ){
        j <- 1
        for( nt_local in c( "A","C","G","T" ) ) {
            nt_counts[ j, h ] <- sum( substr( as.character( barcodes_nt ), h,
                                              h) == nt_local )
            j <- j + 1
        }
    }

    # Calculate frequencies
    nt_freqs <- nt_counts / colSums( nt_counts )

    nt_values <- split(nt_freqs, rep(1:ncol(nt_freqs), each = nrow(nt_freqs)))

    all_posible_comb <- expand.grid( nt_values )

    probs <- apply( all_posible_comb, 1, prod )

    ###Create Mf_to_Sf:
    results <- c()

    i <- 1
    results[ i ] <- sum( 1 - ( 1 - probs )**i )
    j <- 1
    while( results[ i ] <= max_observed ) {
        i <- i + 1
        results[ i ] <- sum( 1 - ( 1 - probs )**i )   #Mf to Sf
        if(results[ i ] == results[ i - 1]){
            if(j > 2)
                stop("Too low complexity to estimate k2n distribution")
            else
                j <- j + 1
        }else
            j <- 1
    }

    #assign molecules number to unique barcode number:
    k2n <- c()
    for( i in 1:floor( max( results ) ) ) {
        abs( results - i ) -> difference

        #if you want to know how many molecules underlie n unique barcodes
        #ask k2n[n]
        k2n[ i ] <- which( difference == min( difference ) )
    }

    #Assign +Inf for fragments which have saturated the barcodes
    #(see correct_oversaturation function):
    if( barcodes_saturated )
        k2n[max_observed + 1] <- Inf

    if(!missing(output_file)) {
        write(k2n, file = output_file )
    } else
        k2n
}

# Read inputs
merged <- args[1]
unique_barcodes <- args[2]
output <- args[3]
paired_check <- as.logical(args[4])

k2n_calc( merged, unique_barcodes, output, paired = paired_check)
