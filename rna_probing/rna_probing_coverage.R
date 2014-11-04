## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
#Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressMessages(library('getopt'))

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical",
  'sample', 's', 1, "character",
  'control', 'c', 2, "character",
  'window', 'w', 1, "integer",
  'ball', 'b', 1, "double",
  'peak', 'p', 1, "integer",
  'output', 'o', 1, "character",
  'counts', 'x', 1, "character",
  'ratio', 'y', 1, "character",
  'priming', 'z', 1, "character",
  'plots', 'l', 2, "character"
  ), byrow=TRUE, ncol=4);

opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}


#This is streamlined protocol that works only with reads mapped to one RNA molecule (16S rRNA) #For plot used in NAR preparation (without two bottom plots) see under many #'s

suppressMessages(library('GenomicRanges'));
#setwd("/home/nikos/rna_probing_galaxy")

#Inputs:
#Options:
window_size <- opt$window
ball_size <- opt$ball
peak <- opt$peak

#Initialize variables:
counts <- list() 							
gene_coverage <- c()
gene_counts <- c()
gene_ratio <- c()
gene_priming <- c()
#End of initializing

####Read in datasets:
counts[[1]] <- read.table( opt$sample )

#if ( !is.null( opt$control ) ) {
#    counts[[2]] <- read.table( opt$control )
#}
#####

#Filter out the short inserts not to deal with size selection correctors.
for(i in seq(1,length(counts))){
  counts[[i]] <- subset(counts[[i]], (counts[[i]][,3]-counts[[i]][,2])>=peak) 
}

############Initilize matrices for storing data:
counter <-1
for(i in seq(1,length(counts))){
  temp_storage <- counts[[i]]
  gene_coverage <- matrix(nrow=max(c(temp_storage[,3], nrow(gene_coverage))), ncol=counter)
  gene_counts <- matrix(nrow=max(c(temp_storage[,3], nrow(gene_counts))), ncol=counter)
  gene_priming <- matrix(nrow=max(c(temp_storage[,3], nrow(gene_priming))), ncol=counter)
  gene_ratio <- matrix(nrow=max(c(temp_storage[,3], nrow(gene_ratio))), ncol=counter)
  counter <- counter+1
}
###########End of Initilize matrices for storing data

##################
##Fill in the matrices with data 
counter <-1
for(i in seq(1,length(counts))){
  temp_storage <- counts[[i]]
  cover <- coverage(IRanges(start=temp_storage[,2], end=(temp_storage[,3]-peak)), weight=temp_storage[,4])
  gene_coverage[1:length((rep(runValue(cover), runLength(cover)))),counter] <- c((rep(runValue(cover), runLength(cover))))
  stopping_reads <- aggregate(temp_storage[,4]~temp_storage[,2],temp_storage, sum)
  gene_counts[stopping_reads[,1],counter] <- stopping_reads[,2]
  gene_counts[is.na(gene_counts)] <- 0
  priming_reads <- aggregate(temp_storage[,4]~temp_storage[,3],temp_storage, sum)
  gene_priming[priming_reads[,1],counter] <- priming_reads[,2]
  gene_priming[is.na(gene_priming)] <- 0
  counter <- counter+1
}

#End of filling in
gene_ratio <- gene_counts/gene_coverage

###Export to Galaxy
data <- data.frame(gene_counts, gene_coverage, gene_priming, gene_ratio )
colnames(data) <- c("Read counts", "Effective Coverage EUC", "Termination EUC", "TCR")

write.table( data, opt$output, sep = "\t", quote = F, row.names = F)


#Return plots
if ( !is.null(opt$plots) ) {
  pdf(opt$plots)
 
  # Termination signal
  plot(gene_counts, type= 'l', main = "Termination signal", ylab = "", xlab = "Position")

  # Priming signal
  plot(gene_priming, type= 'l', main = "Priming signal", ylab = "", xlab = "Position")
  # Effective Coverage

  y <- gene_coverage
  y[is.na(y)] <- 0
  x <- seq_along(y)
  y2 <- rep(y, each=2)
  y2 <- y2[-length(y2)]
  x2 <- rep(x, each=2)[-1]
  x3 <- c(min(x2), x2, max(x2))
  y3 <- c(0, y2, 0)
  
  # because polygon() is dumb and wants a pre-existing plot
  plot(x, y, ylim=c(0, max(y)), type="n", main = "Effective Coverage", ylab = "", xlab = "Position")
  
  polygon(x3, y3, border=NA, col="black")
  lines(x2, y2)
  
  # Termination Coverage Ratio
  plot(gene_ratio, type= 'l', main = "Termination Coverage Ratio", ylab = "", xlab = "Position")
  dump <- dev.off()
}
