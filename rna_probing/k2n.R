options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# Read inputs
merged <- read.table( args[1] )
read_counts <- read.table( args[2] )

barcodes_nt <- merged[ do.call( paste, as.list( merged[ ,1:3 ] ) ) %in% do.call( paste, as.list( read_counts[ ( read_counts[ , 4 ] ) <= summary( read_counts[ , 4 ] )[ 5 ], 1:3 ] ) ) , 4 ]

##make the matrix with the nucleotide freqs per position:
nt_counts <- matrix( nrow = 4, ncol = as.numeric(args[4]) )
for( h in 1:ncol( nt_counts ) ){
    j <- 1
    for( nt_local in c( "A","C","G","T" ) ) {
        nt_counts[ j, h ] <- sum( substr( as.character( barcodes_nt ), h, h) == nt_local )
        j <- j + 1
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

while( results[ i ] <= as.numeric(args[3]) ) {
    i <- i + 1
    results[ i ] <- sum( 1 - ( 1 - probs )**i )   #Mf to Sf
}

#assign molecules number to unique barcode number:
Uf_to_Mf <- c()
for( i in 1:floor( max( results ) ) ) {
abs( results - i ) -> difference
Uf_to_Mf[ i ] <- which( difference == min( difference ) ) #if you want to know how many molecules underlie n unique barcodes ask Uf_to_Mf[n]
}

write( Uf_to_Mf, file = args[5] )
