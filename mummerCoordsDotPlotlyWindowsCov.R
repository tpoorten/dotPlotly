#!/usr/bin/env Rscript

## Make Dot Plot with Percent Divergence on color scale
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

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
  make_option(c("-w", "--coverage_win"), type="numeric", default=100000,
              help="coverage window size (bp) [default %default]",
              dest="coverage_win")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i alignments.coords -o out [options]",option_list=option_list)
opt = parse_args(parser)
# rm(list=ls())
# setwd("~/Google Drive/Experiments/Genome_Assembly/PacBio/runFalcon13/2-asm-falcon12k/")
# opt = list(input_filename="FvescaV4_falcon13_12k.nucmerSet2RM.delta.filter.coords", output_filename="FvescaV4_falcon13_12k.nucmerSet2RM.dotplot", 
#            min_align = 0, min_query_aln = 0, 
#            keep_ref=7, similarity=TRUE, h_lines=F, interactive=FALSE, plot_size=24, on_target = T, v=FALSE,
#            coverage_win=100000)
# rm(list=ls())
# setwd("~/Genome_Assembly/NRGene/scaffolds_v2/Analysis/")
# opt = list(input_filename="FvescaV4_NRGene.nucmerSet2RM.delta.filter.coords",
#            output_filename="FvescaV4_NRGene.nucmerSet2RM.dotplot",
#            min_align = 10000, min_query_aln = 10000,
#            keep_ref=7, similarity=TRUE, h_lines=F, interactive=FALSE, plot_size=24, on_target = T, v=FALSE,
#            coverage_win=200000)

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
  cat(paste0("coverage window size (-w): ", opt$coverage_win,"\n\n"))
}
# read in alignments
alignments = read.table(opt$input_filename, stringsAsFactors = F, skip = 5)
alignments = alignments[,-c(3,6,9,11,14)]

# set column names
colnames(alignments) = c("refStart","refEnd","queryStart","queryEnd","lenAlnRef","lenAlnQuery","percentID","percentAlnCovRef","percentAlnCovQuery","refID","queryID")

cat(paste0("Number of alignments: ", nrow(alignments),"\n"))
cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)),"\n"))

# filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAlnQuery, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]

# filter alignment by length
alignments = alignments[which(alignments$lenAlnQuery > opt$min_align),]

# re-filter queries by alignment length, for now include overlapping intervals
queryLenAgg = tapply(alignments$lenAlnQuery, alignments$queryID, sum)
alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]

# sort by ref chromosome sizes, keep top X chromosomes
chromMax = tapply(alignments$refEnd, alignments$refID, max)
if(is.null(opt$keep_ref)){
  opt$keep_ref = length(chromMax)
  alignments = alignments[which(alignments$refID %in% names(sort(chromMax, decreasing = T)[1:opt$keep_ref])),]
} else{
  alignments = alignments[which(alignments$refID %in% names(sort(chromMax, decreasing = T)[1:opt$keep_ref])),]
}

cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))

# sort df on ref
alignments$refID = factor(alignments$refID, levels = names(sort(chromMax, decreasing = T)[1:opt$keep_ref])) # set order of refID
alignments = alignments[with(alignments,order(refID,refStart)),]

# make new ref alignments for dot plot
chromMax = tapply(alignments$refEnd, alignments$refID, max)
alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names(sort(chromMax, decreasing = T))[1], 0, cumsum(chromMax)[match(x, names(sort(chromMax, decreasing = T))) - 1]) )
alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x) ifelse(x == names(sort(chromMax, decreasing = T))[1], 0, cumsum(chromMax)[match(x, names(sort(chromMax, decreasing = T))) - 1]) )

## queryID sorting step 1/2
# sort levels of factor 'queryID' based on longest alignment
alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))) 
queryMaxAlnIndex = tapply(alignments$lenAlnQuery,
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
queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAlnQuery[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
# set order for queryID
alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])

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

# get mean percent ID per contig
#   calc percent ID based on on-target alignments only
if(opt$on_target){
  alignments$queryTarget = queryID_Ref[match(as.character(alignments$queryID), names(queryID_Ref))]
  alignmentsOnTarget = alignments[which(as.character(alignments$refID) == alignments$queryTarget),]
  scaffoldIDmean = tapply(alignmentsOnTarget$percentID, alignmentsOnTarget$queryID, mean)
  alignments$percentIDmean = scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))]
  alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
} else{
  scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
  alignments$percentIDmean = scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))]
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
          'Query: %s<br>Target: %s<br>Length: %s kb',
          queryID,
          refID,
          round(lenAlnQuery / 1000, 1)
        )
      )
    ) +
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
    scale_y_continuous(breaks = yTickMarks, labels = levels(alignments$queryID)) +
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
    xlab("Reference") +
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
        'Query: %s<br>Target: %s<br>Length: %s kb',
        queryID,
        refID,
        round(lenAlnQuery / 1000, 1)
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
    scale_y_continuous(breaks = yTickMarks, labels = levels(alignments$queryID)) +
    { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                  color = "grey60",
                                  size = .1) }} +
    labs(color = "Mean Percent Identity (per query)", 
         title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                           paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                           paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                           paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
         )) +
    xlab("Reference") +
    ylab("Query")
}
# gp
ggsave(filename = paste0(opt$output_filename, ".png"), width = opt$plot_size, height = opt$plot_size*1.25, units = "in")

if(opt$interactive){
  pdf(NULL)
  gply = ggplotly(gp, tooltip = "text")
  htmlwidgets::saveWidget(as.widget(gply), file = paste0(opt$output_filename, ".html"))
}

options(warn=0) # turn on warnings
#
#
# ////////////////


# ##################
# # Run window analysis - calc. percentID across windows
# library(GenomicRanges)
# binnedAverage2 = 
#   function (bins, numvar, varname) 
#   {
#     if (!is(bins, "GRanges")) 
#       stop("'x' must be a GRanges object")
#     if (!is(numvar, "RleList")) 
#       stop("'numvar' must be an RleList object")
#     if (!identical(seqlevels(bins), names(numvar))) 
#       stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
#     viewMeans2 <- function(v) {
#       means <- viewMeans(v, na.rm = T)
#       w0 <- width(v)
#       w1 <- width(trim(v))
#       means <- means * w1/w0
#       means[w0 != 0L & w1 == 0L] <- 0
#       means
#     }
#     bins_per_chrom <- split(ranges(bins), seqnames(bins))
#     means_list <- lapply(names(numvar), function(seqname) {
#       v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
#       viewMeans2(v)
#     })
#     new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
#     mcols(bins)[[varname]] <- new_mcol
#     bins
#   }
# 
# 
# asmSeqInfo = Seqinfo(names(queryMax), seqlengths=as.numeric(queryMax), isCircular=NA, genome="runFalcon13")
# # add sorted starts and ends
# alignments$queryStartSort = NA
# alignments$queryEndSort = NA
# alignments[,c("queryStartSort", "queryEndSort")] = t(apply(alignments[,c("queryStart","queryEnd")], 1, function(x) sort(x)))
# # alignments$strand = ifelse(alignments$queryStart < alignments$queryEnd, "+", "-")
# # alnGRanges = makeGRangesFromDataFrame(alignments[,c("queryID","queryStartSort","queryEndSort","percentID","strand")], 
# #                                       keep.extra.columns = T, ignore.strand = F, seqinfo = asmSeqInfo, seqnames.field = "queryID", start.field = "queryStartSort", end.field = "queryEndSort", strand.field = "strand")
# # on-target alignments
# alnGRanges = makeGRangesFromDataFrame(alignments[which(!is.na(alignments$percentIDmean)),c("queryID","queryStartSort","queryEndSort","percentID")], 
#                                       keep.extra.columns = T, ignore.strand = T, seqinfo = asmSeqInfo, seqnames.field = "queryID", start.field = "queryStartSort", end.field = "queryEndSort")
# 
# # get disjoint set
# alnGRangesDisjoin = disjoin(alnGRanges, with.revmap=TRUE)
# hist(unlist(lapply(alnGRangesDisjoin$revmap, length)))
# # alnGRangesDisjoin2 = alnGRangesDisjoin[unlist(alnGRangesDisjoin$revmap[which(unlist(lapply(alnGRangesDisjoin$revmap, length)) == 1)])]
# # alnGRangesDisjoin$percentID = NA
# alnGRangesDisjoin$percentID[which(unlist(lapply(alnGRangesDisjoin$revmap, length)) == 1)] = alnGRanges$percentID[unlist(alnGRangesDisjoin$revmap[which(unlist(lapply(alnGRangesDisjoin$revmap, length)) == 1)])]
# # isDisjoint(alnGRangesDisjoin)
# score2 <- mcolAsRleList(alnGRangesDisjoin, "percentID")
# #
# gBins <- tileGenome(asmSeqInfo, tilewidth = 1000000, cut.last.tile.in.chrom = T)
# #
# gBinsAve = binnedAverage2(bins=gBins, numvar=score2, varname="binned_score")
# alignmentsWin = as.data.frame(gBinsAve)
# table(is.na(alignmentsWin$binned_score))
# 
# ## connect to alignments
# pairs <- findOverlapPairs(gBinsAve, alnGRanges, ignore.strand = TRUE)
# ans <- pintersect(pairs, ignore.strand = TRUE)
# ans2 = as.data.frame(ans, row.names = NULL)
# #
# # pairs <- findOverlapPairs(alnGRanges, gBinsAve, ignore.strand = TRUE)
# ###
# # alignments = alignments[which(!is.na(alignments$percentIDmean)),]
# alignments$binned_score = NA
# alignments$binned_score = ans2$binned_score[match(paste(alignments$queryID, alignments$queryStartSort, sep = "_"), paste(ans2$seqnames, ans2$start, sep = "_"))]
# hist(alignments$binned_score, breaks = 100)
# # ////////
# #########
# #
# 
# ggplot(alignments) +
#   geom_point(
#     mapping = aes(x = refStart2, y = queryStart2, color = binned_score),
#     size = 0.009
#   ) +
#   geom_point(
#     mapping = aes(x = refEnd2, y = queryEnd2, color = binned_score),
#     size = 0.009
#   ) +
#   geom_segment(
#     aes(
#       x = refStart2,
#       xend = refEnd2,
#       y = queryStart2,
#       yend = queryEnd2,
#       color = binned_score,
#       text = sprintf(
#         'Query: %s<br>Target: %s<br>Length: %s kb',
#         queryID,
#         refID,
#         round(lenAlnQuery / 1000, 1)
#       )
#     )
#   ) +
#   scale_x_continuous(breaks = cumsum(chromMax),
#                      labels = levels(alignments$refID)) +
#   theme_bw() +
#   theme(text = element_text(size = 8)) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.text.y = element_text(size = 4, angle = 15)
#   ) +
#   scale_y_continuous(breaks = yTickMarks, labels = levels(alignments$queryID)) +
#   { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
#                                 color = "grey60",
#                                 size = .1) }} +
#   scale_color_distiller(palette = "Spectral") +
#   labs(color = "Mean Percent Identity (per query)", 
#        title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
#                          paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
#                          paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
#                          paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
#        )) +
#   xlab("Reference") +
#   ylab("Query")
# opt$output_filename = "FvescaV4_NRGeneDovetail.dotplot.window1Mb"
# # opt$plot_size = 12
# ggsave(filename = paste0(opt$output_filename, ".png"), width = opt$plot_size, height = opt$plot_size*1, units = "in")
# 
# if(opt$interactive){
#   pdf(NULL)
#   gply = ggplotly(gp, tooltip = "text")
#   htmlwidgets::saveWidget(as.widget(gply), file = paste0(opt$output_filename, ".html"))
# }
# 
# #
# # ////////////////

########
# coverage, exclude masked regions
# refSeqInfo = Seqinfo(names(chromMax), seqlengths=as.numeric(chromMax), isCircular=NA, genome="Fvesca_v4")
genome <- import("~/Google Drive/Experiments/Genome_Assembly/Fvesca_Edger/F_vesca_V4.1b.fasta")
refSeqInfo = Seqinfo(seqnames = names(genome)[which(names(genome) %in% levels(alignments$refID))], seqlengths = width(genome)[which(names(genome) %in% levels(alignments$refID))], isCircular=NA, genome="Fvesca_v4")

alnRefGRanges = makeGRangesFromDataFrame(alignments[which(!is.na(alignments$percentIDmean)),c("refID","refStart","refEnd")],
                                         keep.extra.columns = F, ignore.strand = T, seqinfo = refSeqInfo, seqnames.field = "refID", start.field = "refStart", end.field = "refEnd")

refSeqInfoAll = Seqinfo(seqnames = names(genome), seqlengths = width(genome), isCircular=NA, genome="Fvesca_v4")
repMaskGRanges <- import("~/Google Drive/Experiments/Genome_Assembly/Fvesca_Edger/F_vesca_V4.1b.masked.bed", genome = refSeqInfoAll)
repMaskGRanges = intersect(as(refSeqInfo, "GRanges"), unstrand(repMaskGRanges))
nonRepMaskGRanges = setdiff(as(refSeqInfo, "GRanges"), unstrand(repMaskGRanges))
#
alnCoverage = coverage(alnRefGRanges)
alnCoverageGRanges = GRanges(alnCoverage)
alnCoverageGRanges <- subset(alnCoverageGRanges, score > 0)
summary(alnCoverageGRanges$score)
#
# non-repeats coverage
nonRepPairs <- findOverlapPairs(alnCoverageGRanges, nonRepMaskGRanges, ignore.strand = TRUE)
alnCoverageNonRepeatsGRanges <- pintersect(nonRepPairs, ignore.strand = TRUE)
summary(alnCoverageNonRepeatsGRanges$score)
# hist(alnCoverageNonRepeatsGRanges$score)
#
# repeats coverage
repPairs <- findOverlapPairs(alnCoverageGRanges, repMaskGRanges, ignore.strand = TRUE)
alnCoverageRepeatsGRanges <- pintersect(repPairs, ignore.strand = TRUE)
summary(alnCoverageRepeatsGRanges$score)
# hist(alnCoverageRepeatsGRanges$score)
# /////////



##################
# Calc. contig coverage across ref chromosomes

# coverage, average in windows
binnedAverage2 = 
  function (bins, numvar, varname) 
  {
    if (!is(bins, "GRanges")) 
      stop("'x' must be a GRanges object")
    if (!is(numvar, "RleList")) 
      stop("'numvar' must be an RleList object")
    if (!identical(seqlevels(bins), names(numvar))) 
      stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
    viewMeans2 <- function(v) {
      means <- viewMeans(v, na.rm = T)
      w0 <- width(v)
      w1 <- width(trim(v))
      means <- means * w1/w0
      means[w0 != 0L & w1 == 0L] <- 0
      means
    }
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar), function(seqname) {
      v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
      viewMeans2(v)
    })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[varname]] <- new_mcol
    bins
  }
#
#
#####
# coverage, average in windows, cut intervals with 0 coverage, tuned for repeat-masked ref
genome <- import("~/Google Drive/Experiments/Genome_Assembly/Fvesca_Edger/F_vesca_V4.1b.fasta")
refSeqInfo = Seqinfo(seqnames = names(genome)[which(names(genome) %in% levels(alignments$refID))], seqlengths = width(genome)[which(names(genome) %in% levels(alignments$refID))], isCircular=NA, genome="Fvesca_v4")

alnRefGRanges = makeGRangesFromDataFrame(alignments[which(!is.na(alignments$percentIDmean)),c("refID","refStart","refEnd")],
                                         keep.extra.columns = F, ignore.strand = T, seqinfo = refSeqInfo, seqnames.field = "refID", start.field = "refStart", end.field = "refEnd")

alnCoverage = coverage(alnRefGRanges)

alnCoverageGRanges = GRanges(alnCoverage)
alnCoverageGRanges <- subset(alnCoverageGRanges, score > 0)

# summary(alnCoverageGRanges$score)
# hist(alnCoverageGRanges$score)
cov2 <- mcolAsRleList(alnCoverageGRanges, "score")
# opt$coverage_win = 200000
refBins <- tileGenome(refSeqInfo, tilewidth = opt$coverage_win, cut.last.tile.in.chrom = T)
#
refBinsAve = binnedAverage2(bins=refBins, numvar=cov2, varname="binned_cov")
refBinsAveDF = as.data.frame(refBinsAve)

summary(refBinsAveDF$binned_cov)
png(paste0(opt$output_filename, "_coverage_window",opt$coverage_win/1000,"k.png"))
hist(refBinsAveDF$binned_cov[which(refBinsAveDF$binned_cov < 10)], breaks=40, main = paste0("median = ",round(median(refBinsAveDF$binned_cov, na.rm=T),3)))
dev.off()
