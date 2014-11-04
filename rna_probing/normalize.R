##!/usr/bin/Rscript

## Setup R error handling to go to stderr
options( show.error.messages = FALSE, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressMessages(library('getopt'))

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'treated', 't', 1, "character",
  'control', 'c', 2, "character",
  'method', 'm', 1, "character",
  'fuComplexity', 'b', 2, "integer",
  'k2nTreated', 'k', 2, "character",
  'k2nControl', 'a', 2, "character",
  'cutoff', 'g', 1, "integer",
  'reference', 'r', 2, "character",
  'compdata', 'h', 2, "character",
  'dtcr', 'd', 0, "logical",
  'dtcrWindow', 'e', 2, "integer",
  'dtcrToZero', 'f', 2, "character",
  'slograt', 's', 0, "logical",
  'slogratWindow', 'p', 2, "integer",
  'depthCorrection', 'q', 2, "character",
  'pseudocount', 'l', 2, "integer",
  'swinsor', 'w', 0, "logical",
  'swinsorWindow', 'x', 2, "integer",
  'winsorLevel', 'y', 2, "double",
  'fixQuantile', 'z', 2, "character",
  'ntOffset', 'n', 2, "integer",
  'outputDir', 'o', 1, "character",
  'bedgraph', 'v', 2, "character",
  'bed', 'j', 2, "character",
  'genome', 'i', 2, "character",
  'trackName', 'u', 2, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

suppressMessages(require(RNAprobR))


#Create output dir
dir.create(opt$outputDir, showWarnings = FALSE, recursive = TRUE)


if (opt$method=="counts"){
    treated_euc <- readsamples(opt$treated, euc=opt$method)
} else if (opt$method=="Fu") {
    treated_euc <- readsamples(opt$treated, euc=opt$method, m = opt$fuComplexity)
} else {
    treated_euc <- readsamples(opt$treated, euc=opt$method, k2n_files=opt$k2nTreated)
}


if ( !is.null(opt$reference)){
    comp_treated <- comp(treated_euc, cutoff=opt$cutoff,
                         fasta_file=opt$reference)
} else {
    comp_treated <- comp(treated_euc, cutoff=opt$cutoff)
}


#If present, read control file
if ( !is.null(opt$control) && opt$control != 'None'){
    if (opt$method=="counts"){
        control_euc <- readsamples(opt$control, euc=opt$method)
    } else if (opt$method=="Fu") {
        control_euc <- readsamples(opt$control, euc=opt$method, m = opt$fuComplexity)
    } else {
        control_euc <- readsamples(opt$control, euc=opt$method, k2n_files=opt$k2nControl)
    }
    if ( !is.null(opt$reference)){
        comp_control <- comp(control_euc, cutoff=opt$cutoff,
                         fasta_file=opt$reference)
    } else {
        comp_control <- comp(control_euc, cutoff=opt$cutoff)
    }
}

#compdata
if ( !is.null(opt$compdata) ) {

    normalized <- compdata(comp_treated, nt_offset=opt$ntOffset)
    colnames(mcols(normalized))[1:4] <- paste(colnames(mcols(normalized))[1:4], ".treated", sep="")

    if ( !is.null(opt$control) && opt$control != 'None') {
        normalized <- compdata(comp_control, nt_offset=opt$ntOffset, add_to=normalized)
        colnames(mcols(normalized))[6:9] <- c("TC.control", "TCR.control", "Cover.control", "PC.control")
    }
}

#dtcr
if ( !is.null(opt$dtcr) ) {
    if ( exists("normalized")) {
        normalized <- dtcr(comp_control, comp_treated, window_size=opt$dtcrWindow,
                           nt_offset=opt$ntOffset, bring_to_zero=opt$dtcrToZero,
                           add_to=normalized)
    }
    else {
        normalized <- dtcr(comp_control, comp_treated, window_size=opt$dtcrWindow,
                           nt_offset=opt$ntOffset, bring_to_zero=opt$dtcrToZero)
    }
}

#slograt
if ( !is.null(opt$slograt) ) {
    if ( exists("normalized")) {
        normalized <- slograt(comp_control, comp_treated, window_size=opt$slogratWindow,
                           nt_offset=opt$ntOffset, depth_correction = opt$depthCorrection,
                           pseudocount=opt$pseudocount, add_to=normalized)
    }
    else {
        normalized <- slograt(comp_control, comp_treated, window_size=opt$slogratWindow,
                              nt_offset=opt$ntOffset, depth_correction = opt$depthCorrection,
                              pseudocount=opt$pseudocount)
    }
}

#swinsor
if ( !is.null(opt$swinsor) ) {
    if ( exists("normalized")) {
        normalized <- swinsor(comp_treated, winsor_level = opt$winsorLevel,
                              window_size=opt$swinsorWindow,
                              only_top=opt$fixQuantile, nt_offset=opt$ntOffset,
                              add_to=normalized)
    }
    else {
        normalized <- swinsor(comp_treated, winsor_level = opt$winsorLevel,
                              window_size=opt$swinsorWindow,
                              only_top=opt$fixQuantile, nt_offset=opt$ntOffset)
    }
}

#bedgraph output
if ( !is.null(opt$bedgraph)) {
    if ( !is.null(opt$dtcr)){
        norm2bedgraph(normalized, bed_file = opt$bed, norm_method = "dtcr",
                      genome_build = opt$genome,
                      bedgraph_out_file = paste(opt$outputDir, "/dtcr",
                                                sep=""),
                      track_name = opt$trackName,
                      track_description = opt$trackDesc)
    }
    if ( !is.null(opt$slograt)){
        norm2bedgraph(normalized, bed_file = opt$bed, norm_method = "slograt",
                      genome_build = opt$genome,
                      bedgraph_out_file = paste(opt$outputDir, "/slograt",
                                                sep=""),
                      track_name = opt$trackName,
                      track_description = opt$trackDesc)
    }
    if ( !is.null(opt$swinsor)){
        norm2bedgraph(normalized, bed_file = opt$bed, norm_method = "swinsor",
                      genome_build = opt$genome,
                      bedgraph_out_file = paste(opt$outputDir, "/swinsor",
                                                sep=""),
                      track_name = opt$trackName,
                      track_description = opt$trackDesc)
    }
}

output <- GR2norm_df(normalized)
write.table( output, paste(opt$outputDir, "/norm_df.txt", sep = ""), sep = "\t",
             quote = F, row.names = F)
