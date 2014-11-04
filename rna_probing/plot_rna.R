#Setup R error handling to go to stderr
options( show.error.messages = FALSE, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressMessages(library('getopt'))

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'help' , 'h', 0, "logical",
    'input', 'i', 1, "character",
    'transcript', 't', 1, "character",
    'method', 'm', 2, "character",
    'cutoff', 'c', 2, "double",
    'type', 'p', 2, "character",
    'output', 'o', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

suppressMessages(require(RNAstr))

#Read and convert input to GRanges object
data <- read.table(opt$input, header = TRUE)
dataGR <- norm_df2GR(data)

#Check if given transcript exists in input file
if ( ! opt$transcript %in% data$RNAid ) { stop("Transript not found. Check input file.") }

pdf(opt$output)
for (method in strsplit(opt$method, ",")[[1]] ) {

    #Check if columns exists in data file
    if (! method %in% colnames(data)) { next } 

    plotRNA(dataGR, opt$transcript, method, stat_cutoff = opt$cutoff, type = opt$type,
        main=paste(opt$transcript,": ", method, sep=""))
}
dev.off()
