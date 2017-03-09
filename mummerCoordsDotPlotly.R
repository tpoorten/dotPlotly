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
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=10000,
              help="filter alignments less than cutoff X bp [default %default]",
              dest="min_align"),
  make_option(c("-k", "--number-ref-chromosomes"), type="numeric", default=NULL,
              help="number of sorted reference chromosomes to keep [default all chromosmes]",
              dest="keep_ref"),
  make_option(c("-s", "--similarity"), action="store_true", default=FALSE,
              help="turn on color alignments by percent similarity [default %default]",
              dest="similarity")
  )

options(error=traceback)

parser <- OptionParser(usage = "%prog -i alignments.coords -o out [options]",option_list=option_list)
opt = parse_args(parser)
# opt = list(input_filename="FvescaV4_v_p_ctg.nucmer.coords", output_filename="out", min_align=15000, keep_ref=7, similarity=TRUE)
# print(opt)

# read in alignments
alignments = read.table(opt$input_filename, stringsAsFactors = F, skip = 5)
alignments = alignments[,-c(3,6,9,11,14)]

# set column names
colnames(alignments) = c("refStart","refEnd","queryStart","queryEnd","lenAlnRef","lenAlnQuery","percentID","percentAlnCovRef","percentAlnCovQuery","refID","queryID")

print(paste("Number of alignments:", nrow(alignments)))
print(paste("Number of query sequences:", length(unique(alignments$queryID))))

# filter alignment by length
alignments = alignments[which(alignments$lenAlnQuery > opt$min_align),]

# sort by ref chromosome sizes, keep top X chromosomes
chromMax = tapply(alignments$refEnd, alignments$refID, max)
if(is.null(opt$keep_ref)){
  opt$keep_ref = length(chromMax)
  alignments = alignments[which(alignments$refID %in% names(sort(chromMax, decreasing = T)[1:opt$keep_ref])),]
} else{
  alignments = alignments[which(alignments$refID %in% names(sort(chromMax, decreasing = T)[1:opt$keep_ref])),]
}

print(paste("After filtering... Number of alignments:", nrow(alignments)))
print(paste("After filtering... Number of query sequences:", length(unique(alignments$queryID))))

# sort df on ref
alignments$refID = factor(alignments$refID, levels = names(sort(chromMax, decreasing = T)[1:opt$keep_ref])) # set order of refID
alignments = alignments[with(alignments,order(refID,refStart)),]

# get mean percent ID per contig
scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
alignments$percentIDmean = scaffoldIDmean[match(alignments$queryID, names(scaffoldIDmean))]

# make new ref alignments for dot plot
chromMax = tapply(alignments$refEnd, alignments$refID, max)
alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names(sort(chromMax, decreasing = T))[1], 0, cumsum(chromMax)[match(x, names(sort(chromMax, decreasing = T))) - 1]) )
alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x) ifelse(x == names(sort(chromMax, decreasing = T))[1], 0, cumsum(chromMax)[match(x, names(sort(chromMax, decreasing = T))) - 1]) )

# sort levels of factor 'refId' based on longest alignment
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))) 
queryMaxAlnIndex = tapply(alignments$lenAlnQuery, alignments$queryID, which.max, simplify = F)
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))[order(mapply(function(x, i) alignments$refStart2[which(i==alignments$queryID)][x], queryMaxAlnIndex, names(queryMaxAlnIndex)))]) 

#  flip query starts stops to forward if most align are in reverse complement
queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
queryRevComp = names(queryRevComp)[which(queryRevComp)]
queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1

# make new query alignments for dot plot
queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
names(queryMax) = levels(alignments$queryID)
alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
alignments$queryEnd2 = alignments$queryEnd +     sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

# plot
options(warn=-1) # turn off warnings
if(opt$similarity){
  gp = ggplot(alignments) + 
    geom_point(mapping = aes(x=refStart2, y=queryStart2, colour=percentIDmean), size=0.5) +
    geom_point(mapping = aes(x=refEnd2, y=queryEnd2, colour=percentIDmean), size=0.5) +
    geom_segment(aes(x = refStart2, xend = refEnd2, y = queryStart2, yend = queryEnd2, 
                     color = percentIDmean, 
                     text=sprintf('Query: %s<br>Target: %s<br>Length: %s kb', queryID, refID, round(lenAlnQuery/1000,1)))) +  
    scale_x_continuous(breaks = cumsum(chromMax), labels = levels(alignments$refID)) +
    theme(text = element_text(size = 8)) +
    theme_bw() +
    scale_color_distiller(palette = "Spectral") +
    xlab("Reference") +
    ylab("Query")
} else {
  gp = ggplot(alignments) + 
    geom_point(mapping = aes(x=refStart2, y=queryStart2), size=0.5) +
    geom_point(mapping = aes(x=refEnd2, y=queryEnd2), size=0.5) +
    geom_segment(aes(x = refStart2, xend = refEnd2, y = queryStart2, yend = queryEnd2, 
                     text=sprintf('Query: %s<br>Target: %s<br>Length: %s kb', queryID, refID, round(lenAlnQuery/1000,1)))) +  
    scale_x_continuous(breaks = cumsum(chromMax), labels = levels(alignments$refID)) +
    theme(text = element_text(size = 8)) +
    scale_color_distiller(palette = "Spectral") +
    xlab("Reference") +
    ylab("Query")
}
# gp
ggsave(filename=paste0(opt$output_filename, ".pdf"), width=10, height=10)

pdf(NULL)
gply = ggplotly(gp, tooltip = "text")
htmlwidgets::saveWidget(as.widget(gply), file = paste0(opt$output_filename, ".html"))

options(warn=0) # turn on warnings

#
#
# ////////////////
