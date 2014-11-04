suppressMessages(library('getopt'))

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix( c(
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical",
  'merged', 'm', 1, "character",
  'read_counts', 'c', 1, "character",
  'max_observed', 'b', 1, "integer",
  'barcode_length', 'l', 1, "integer",
  'output', 'o', 1, "character"
  ), byrow = TRUE, ncol = 4 );

opt = getopt( spec );

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null( opt$help ) ) {
    cat( getopt( spec, usage=TRUE ) ) ;
    q( status = 1 );
}

# Read inputs
merged <- read.table( opt$merged )
read_counts <- read.table( opt$read_counts )

barcodes_nt <- merged[ do.call( paste, as.list( merged[ ,1:3 ] ) ) %in% do.call( paste, as.list( read_counts[ ( read_counts[ , 4 ] ) <= summary( read_counts[ , 4 ] )[ 5 ], 1:3 ] ) ) , 4 ]

##make the matrix with the nucleotide freqs per position:
nt_counts <- matrix( nrow = 4, ncol = opt$barcode_length )
for( h in 1:ncol( nt_counts ) ){
    j <- 1
    for( nt_local in c( "A","C","G","T" ) ) {
        nt_counts[ j, h ] <- sum( substr( as.character( barcodes_nt ), h, h) == nt_local )
	j <- j+1
    }
}

# Calculate frequencies
nt_freqs <- nt_counts / colSums( nt_counts )

nt_values <- list()
for( i in 1:ncol( nt_freqs ) ) {
    nt_values[[ i ]] <- nt_freqs[ , i ]
}


all_posible_comb <- expand.grid( nt_values )

probs <- apply( all_posible_comb, 1, prod )

###Create Mf_to_Sf:
results <- c()

i <- 1
results[ i ] <- sum( 1 - ( 1 - probs )**i )

while( results[ i ] <= opt$max_observed ) {
    i <- i + 1
    results[ i ] <- sum( 1 - ( 1 - probs )**i )   #Mf to Sf
}

#assign molecules number to unique barcode number:
Uf_to_Mf <- c()
for( i in 1:floor( max( results ) ) ) {
abs( results - i ) -> difference
Uf_to_Mf[ i ] <- which( difference == min( difference ) ) #if you want to know how many molecules underlie n unique barcodes ask Uf_to_Mf[n]
}

write( Uf_to_Mf, file = opt$output )
