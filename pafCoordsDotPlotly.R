#!/usr/bin/env Rscript

## Make Dot Plot with Percent Divergence on color scale
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))

option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="coords file from mummer program 'show.coords' [default %default]",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]"),
  make_option(c("-q", "--min-query-length"), type="numeric", default=400000,
              help="filter queries with total alignments less than cutoff X bp [default %default]",
              dest="min_query_aln"),
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=10000,
              help="filter alignments less than cutoff X bp [default %default]",
              dest="min_align"),
  make_option(c("-p","--plot-size"), type="numeric", default=15,
              help="plot size X by X inches [default %default]",
              dest="plot_size"),
  make_option(c("-l", "--show-horizontal-lines"), action="store_true", default=FALSE,
              help="turn on horizontal lines on plot for separating scaffolds  [default %default]",
              dest="h_lines"),
  make_option(c("-k", "--number-ref-chromosomes"), type="numeric", default=NULL,
              help="number of sorted reference chromosomes to keep [default all chromosmes]",
              dest="keep_ref"),
  make_option(c("-s", "--identity"), action="store_true", default=FALSE,
              help="turn on color alignments by % identity [default %default]",
              dest="similarity"),
  make_option(c("-t", "--identity-on-target"), action="store_true", default=FALSE,
              help="turn on calculation of % identity for on-target alignments only [default %default]",
              dest="on_target"),
  make_option(c("-x", "--interactive-plot-off"), action="store_false", default=TRUE,
              help="turn off production of interactive plotly [default %default]",
              dest="interactive"),
  make_option(c("-r", "--reference-ids"), type="character", default=NULL,
              help="comma-separated list of reference IDs to keep [default %default]",
              dest="refIDs")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i alignments.coords -o out [options]",option_list=option_list)
opt = parse_args(parser)

# rm(list=ls())
# setwd("~/GitHubLocal/dotPlotly/testing/minimap_paf/")
# opt = list(input_filename="susie_ggo_grch38.minimap0.txt",
#            output_filename="testOut.allRef",
#            min_align = 500, min_query_aln = 500000,
#            keep_ref=23, 
#            similarity=T, h_lines=T, interactive=F, plot_size=15, on_target = T, v=FALSE)


if(opt$v){
  cat(paste0("PARAMETERS:\ninput (-i): ", opt$input_filename,"\n"))
  cat(paste0("output (-o): ", opt$output_filename,"\n"))
  cat(paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln,"\n"))
  cat(paste0("minimum alignment length (-m): ", opt$min_align,"\n"))
  cat(paste0("plot size (-p): ", opt$plot_size,"\n"))
  cat(paste0("show horizontal lines (-l): ", opt$h_lines,"\n"))
  cat(paste0("number of reference chromosomes to keep (-k): ", opt$keep_ref,"\n"))
  cat(paste0("show % identity (-s): ", opt$similarity,"\n"))
  cat(paste0("show % identity for on-target alignments only (-t): ", opt$similarity,"\n"))
  cat(paste0("produce interactive plot (-x): ", opt$interactive,"\n"))
  cat(paste0("reference IDs to keep (-r): ", opt$refIDs,"\n"))
}
opt$output_filename = unlist(strsplit(opt$output_filename, "/"))[length(unlist(strsplit(opt$output_filename, "/")))]

# read in alignments
alignments = read.table(opt$input_filename, stringsAsFactors = F, fill = T)

# set column names
# PAF IS ZERO-BASED - CHECK HOW CODE WORKS
colnames(alignments)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")

# Fixes for PAF
# Some measure of similarity - need to check on this
alignments$percentID = alignments$numResidueMatches / alignments$lenAln
queryStartTemp = alignments$queryStart
# Flip starts, ends for negative strand alignments
alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
rm(queryStartTemp)

cat(paste0("\nNumber of alignments: ", nrow(alignments),"\n"))
cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)),"\n"))

# sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
if(is.null(opt$refIDs)){
  chromMax = tapply(alignments$refEnd, alignments$refID, max)
  if(is.null(opt$keep_ref)){
    opt$keep_ref = length(chromMax)
  }
  refIDsToKeepOrdered = names(sort(chromMax, decreasing = T)[1:opt$keep_ref])
  alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
  
} else {
  refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
  alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
}

# filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]

# filter alignment by length
alignments = alignments[which(alignments$lenAln > opt$min_align),]

# re-filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]

cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))

# sort df on ref
alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
alignments = alignments[with(alignments,order(refID,refStart)),]
chromMax = tapply(alignments$refEnd, alignments$refID, max)

# make new ref alignments for dot plot
if(length(levels(alignments$refID)) > 1){
  alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
} else {
  alignments$refStart2 = alignments$refStart
  alignments$refEnd2 = alignments$refEnd
  
}

## queryID sorting step 1/2
# sort levels of factor 'queryID' based on longest alignment
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))) 
queryMaxAlnIndex = tapply(alignments$lenAln,
                          alignments$queryID,
                          which.max,
                          simplify = F)
alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
  function(x, i)
    alignments$refStart2[which(i == alignments$queryID)][x],
  queryMaxAlnIndex,
  names(queryMaxAlnIndex)
))])

## queryID sorting step 2/2
## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
# per query ID, get aggregrate alignment length to each refID 
queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
if(length(levels(alignments$refID)) > 1){
  queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
} else {queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])}
# set order for queryID
alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])

#  flip query starts stops to forward if most align are in reverse complement
queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
queryRevComp = names(queryRevComp)[which(queryRevComp)]
queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1

## make new query alignments for dot plot
# subtract queryStart and Ends by the minimum alignment coordinate + 1
queryMin = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), min)
names(queryMin) = levels(alignments$queryID)
alignments$queryStart = as.numeric(alignments$queryStart - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)
alignments$queryEnd = as.numeric(alignments$queryEnd - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)

queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
alignments$queryEnd2 = alignments$queryEnd +     sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

# get mean percent ID per contig
#   calc percent ID based on on-target alignments only
if(opt$on_target & length(levels(alignments$refID)) > 1){
  alignments$queryTarget = queryID_Ref[match(as.character(alignments$queryID), names(queryID_Ref))]
  alignmentsOnTarget = alignments[which(as.character(alignments$refID) == alignments$queryTarget),]
  scaffoldIDmean = tapply(alignmentsOnTarget$percentID, alignmentsOnTarget$queryID, mean)
  alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
  alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
} else{
  scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
  alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
}

# plot
yTickMarks = tapply(alignments$queryEnd2, alignments$queryID, max)
options(warn = -1) # turn off warnings
if (opt$similarity) {
  gp = ggplot(alignments) +
    geom_point(
      mapping = aes(x = refStart2, y = queryStart2, color = percentIDmean),
      size = 0.009
    ) +
    geom_point(
      mapping = aes(x = refEnd2, y = queryEnd2, color = percentIDmean),
      size = 0.009
    ) +
    geom_segment(
      aes(
        x = refStart2,
        xend = refEnd2,
        y = queryStart2,
        yend = queryEnd2,
        color = percentIDmean,
        text = sprintf(
          'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
          queryID,
          queryStart,
          queryEnd,
          refID,
          refStart,
          refEnd,
          round(lenAln / 1000, 1)
        )
      )
    ) +
    scale_x_continuous(breaks = cumsum(as.numeric(chromMax)),
                       labels = levels(alignments$refID)) +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 4, angle = 15)
    ) +
    scale_y_continuous(breaks = yTickMarks, labels = substr(levels(alignments$queryID), start = 1, stop = 20)) +
    { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                  color = "grey60",
                                  size = .1) }} +
    scale_color_distiller(palette = "Spectral") +
    labs(color = "Mean Percent Identity (per query)", 
         title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                           paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                           paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                           paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
         )) +
    xlab("Target") +
    ylab("Query")
} else {
  gp = ggplot(alignments) +
    geom_point(mapping = aes(x = refStart2, y = queryStart2),
               size = 0.009) +
    geom_point(mapping = aes(x = refEnd2, y = queryEnd2),
               size = 0.009) +
    geom_segment(aes(
      x = refStart2,
      xend = refEnd2,
      y = queryStart2,
      yend = queryEnd2,
      text = sprintf(
        'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
        queryID,
        queryStart,
        queryEnd,
        refID,
        refStart,
        refEnd,
        round(lenAln / 1000, 1)
      )
    )) +
    scale_x_continuous(breaks = cumsum(chromMax),
                       labels = levels(alignments$refID)) +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 4, angle = 15)
    ) +
    scale_y_continuous(breaks = yTickMarks, labels = substr(levels(alignments$queryID), start = 1, stop = 20)) +
    { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                  color = "grey60",
                                  size = .1) }} +
    labs(color = "Mean Percent Identity (per query)", 
         title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                           paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                           paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                           paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
         )) +
    xlab("Target") +
    ylab("Query")
}
# gp
ggsave(filename = paste0(opt$output_filename, ".png"), width = opt$plot_size, height = opt$plot_size, units = "in", dpi = 300, limitsize = F)

if(opt$interactive){
  pdf(NULL)
  gply = ggplotly(gp, tooltip = "text")
  htmlwidgets::saveWidget(as.widget(gply), file = paste0(opt$output_filename, ".html"))
}

options(warn=0) # turn on warnings
#
