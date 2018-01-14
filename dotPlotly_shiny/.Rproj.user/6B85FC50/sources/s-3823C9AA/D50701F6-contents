library(shiny)
library(ggplot2)
library(plotly)

# set max upload size
options(shiny.maxRequestSize = 30*1024^2)

ui = fluidPage(
  sidebarPanel(
    fileInput("file1", "Choose PAF File",
              accept = c(
                "paf")
    ),
    tags$hr(),
    numericInput("num", 
                 h3("Number of reference scaffolds to keep"), 
                 value = 2),
    numericInput("min_query_aln", 
                 h3("Filter queries with total alignments less than cutoff X bp"), 
                 value = 500000),
    numericInput("min_align", 
                 h3("Filter alignments less than cutoff X bp"), 
                 value = 1000),
    checkboxInput("similarity", "Show % identity with color gradient", value = TRUE),
    checkboxInput("h_lines", "Show horizontal lines between scaffolds", value = TRUE)
  ),
  mainPanel(
    plotlyOutput('plot1'),
    verbatimTextOutput("event")
  )
)

server = function(input, output) {
  dat <- reactive({
    infile <- input$file1
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    } else {
      opt = list(refIDs = NULL, keep_ref = input$num, 
                 min_query_aln = input$min_query_aln, 
                 min_align = input$min_align, 
                 similarity = input$similarity,
                 on_target = TRUE,
                 h_lines = input$h_lines)
      alignments = read.table(infile$datapath, header=F, stringsAsFactors = F, fill = T, row.names = NULL)
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
        scaffoldIDmean = sapply(unique(as.character(alignmentsOnTarget$queryID)), function(x)
          weighted.mean(x = alignmentsOnTarget$percentID[which(as.character(alignmentsOnTarget$queryID) == x)], 
                        w = alignmentsOnTarget$lenAln[which(as.character(alignmentsOnTarget$queryID) == x)])
        )
        alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
        alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
      } else{
        scaffoldIDmean = sapply(unique(as.character(alignments$queryID)), function(x)
          weighted.mean(x = alignments$percentID[which(as.character(alignments$queryID) == x)], 
                        w = alignments$lenAln[which(as.character(alignments$queryID) == x)])
        )
        alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
      }
      
      results = list(aln = alignments, opt = opt)
      return(results)
    }
  })

  output$plot1 <- renderPlotly({
    if (is.null(input$file1)) {
      # User has not uploaded a file yet
      return(NULL)
    } else {
      mydataList = dat()
      mydata = mydataList$aln
      opt = mydataList$opt
      #
      chromMax = tapply(mydata$refEnd, mydata$refID, max)
      yTickMarks = tapply(mydata$queryEnd2, mydata$queryID, max)
      print(yTickMarks)
      print(length(yTickMarks))
      print(length(substr(levels(mydata$queryID), start = 1, stop = 20)))
      options(warn = -1) # turn off warnings
      if (opt$similarity) {
        gp = ggplot(mydata) +
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
                round(lenAln / 1000, 1)
              )
            )
          ) +
          scale_x_continuous(breaks = cumsum(as.numeric(chromMax)),
                             labels = levels(mydata$refID)) +
          theme_bw() +
          theme(text = element_text(size = 8)) +
          theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y = element_text(size = 4, angle = 15),
            legend.position = "top"
          ) +
          scale_y_continuous(breaks = yTickMarks, labels = substr(levels(mydata$queryID), start = 1, stop = 20)) +
          { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                        color = "grey60",
                                        size = .1) }} +
          scale_color_distiller(palette = "Spectral") +
          labs(color = "Mean Percent Identity (per query)", 
               title = paste0(   paste0("Post-filtering number of alignments: ", nrow(mydata),"\t\t\t\t"),
                                 paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                                 paste0("Post-filtering number of queries: ", length(unique(mydata$queryID)),"\t\t\t\t\t\t\t\t"),
                                 paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
               )) +
          xlab("Reference") +
          ylab("Query")
      } else {
        gp = ggplot(mydata) +
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
              round(lenAln / 1000, 1)
            )
          )) +
          scale_x_continuous(breaks = cumsum(chromMax),
                             labels = levels(mydata$refID)) +
          theme_bw() +
          theme(text = element_text(size = 8)) +
          theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y = element_text(size = 4, angle = 15)
          ) +
          scale_y_continuous(breaks = yTickMarks, labels = substr(levels(mydata$queryID), start = 1, stop = 20)) +
          { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                        color = "grey60",
                                        size = .1) }} +
          labs(color = "Mean Percent Identity (per query)", 
               title = paste0(   paste0("Post-filtering number of alignments: ", nrow(mydata),"\t\t\t\t"),
                                 paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                                 paste0("Post-filtering number of queries: ", length(unique(mydata$queryID)),"\t\t\t\t\t\t\t\t"),
                                 paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
               )) +
          xlab("Reference") +
          ylab("Query")
      }
      
      }
  })
}
# //////////////////
shinyApp(ui = ui, server = server)
