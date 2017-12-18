library(shiny)
library(DT)
library(ggplot2)
library(ggbeeswarm)
library(plyr)
library(plotly)

load("Data.RData")

autoBoxplot <-
  function(gene,
           clusList,
           tissueList,
           patientList,
           plotBy,
           plotBee,
           dotSize,
           plotColor,
           plotLabel,
           plotProg = NULL) {
    maiN = paste(Genes[gene], "Expression")
    resulT <- NULL
    try(plotProg$set(message = "Filtering data...", value = 0.2),
        silent = TRUE)
    crit <- Sample.Def$MajorCellType %in% clusList &
      Sample.Def$tissueSource %in% tissueList &
      Sample.Def$patient %in% patientList
    dataPool <- data.frame(Expression = Sample.Exp[crit, gene])
    ylaB = "Normalized expression"
    
    try(plotProg$set(message = "Combining data...", value = 0.4),
        silent = TRUE)
    if (plotBy == 1) {
      dataPool$XINDEX <- Sample.Def$MajorCellType[crit]
      xlaB = "Cluster"
    }
    if (plotBy == 2) {
      dataPool$XINDEX <- Sample.Def$tissueSource[crit]
      xlaB = "Tissue type"
    }
    
    try(plotProg$set(message = "Generating plot...", value = 0.6),
        silent = TRUE)
    resulT <-
      ggplot(data = dataPool, aes(x = XINDEX, y = Expression)) +
      geom_boxplot(
        outlier.alpha = 0.5,
        outlier.color = "#a3a3a3",
        outlier.size = dotSize
      ) +
      theme_bw() + ggtitle(maiN) + ylab(ylaB) + xlab(xlaB) +
      theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 10,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
      )
    if (plotBee == TRUE)
      if (plotColor == TRUE)
        resulT <-
      resulT + geom_quasirandom(
        cex = dotSize,
        aes(col = XINDEX),
        show.legend = FALSE,
        width = 0.25,
        alpha = 0.5
      )
    else
      resulT <-
      resulT + geom_quasirandom(
        cex = dotSize,
        col = "#a3a3a3",
        show.legend = FALSE,
        width = 0.25,
        alpha = 0.5
      )
    if (plotLabel == TRUE) {
      meds <-
        ddply(dataPool, .(XINDEX), summarise, med = median(Expression))
      resulT <-
        resulT + geom_text(
          data = meds,
          aes(
            label = round(med, 2),
            x = XINDEX,
            y = med
          ),
          size = 5,
          vjust = -0.5
        )
    }
    return(resulT)
  }

autoTsne <-
  function(gene,
           clusList,
           tissueList,
           patientList,
           dotSize,
           plotProg = NULL) {
    resulT <- NULL
    try(plotProg$set(message = "Filtering data...", value = 0.2),
        silent = TRUE)
    maiN = paste(Genes[gene], "Expression in tSNE Clusters")
    crit <- Sample.Def$MajorCellType %in% clusList &
      Sample.Def$tissueSource %in% tissueList &
      Sample.Def$patient %in% patientList
    
    dataPool <- data.frame(Expression = Sample.Exp[crit, gene],
                           x = Sample.tSNE[crit, 1],
                           y = Sample.tSNE[crit, 2])
    llaB <- "Norm. Exp."
    
    try(plotProg$set(message = "Generating plot...", value = 0.4),
        silent = TRUE)
    
    dataPool$ID <- rownames(Sample.tSNE)[crit]
    dataPool$Cluster <- Sample.Def$MajorCellType[crit]
    resulT <-
      ggplot(data = dataPool, aes(
        x = x,
        y = y,
        color = Expression,
        text = paste("Cell ID: ", ID, "\nCluster: ", Cluster)
      )) + geom_point(alpha = 0.7, size =
                        dotSize) +
      scale_color_gradientn(
        colors = c(
          "#ffffcc",
          "#ffeda0",
          "#fed976",
          "#feb24c",
          "#fd8d3c",
          "#fc4e2a",
          "#e31a1c",
          "#bd0026",
          "#800026"
        )
      ) +
      theme_bw() + ggtitle(maiN) +
      theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        #legend.position = c(0.99, 0.99),
        #legend.justification = c(1, 1),
        legend.text = element_text(size = 10),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
      ) + labs(color = llaB)
    return(resulT)
  }

autoTsneSet <-
  function(geneset,
           clusList,
           tissueList,
           patientList,
           dotSize,
           plotProg = NULL) {
    resulT <- NULL
    try(plotProg$set(message = "Filtering data...", value = 0.2),
        silent = TRUE)
    maiN = paste("Geneset Expression in tSNE Clusters")
    crit <- Sample.Def$MajorCellType %in% clusList &
      Sample.Def$tissueSource %in% tissueList &
      Sample.Def$patient %in% patientList
    
    dataPool <- data.frame(Expression = rowSums(Sample.Exp[crit, geneset]),
                           x = Sample.tSNE[crit, 1],
                           y = Sample.tSNE[crit, 2])
    llaB <- "Norm. Exp."
    
    try(plotProg$set(message = "Generating plot...", value = 0.4),
        silent = TRUE)
    
    dataPool$ID <- rownames(Sample.tSNE)[crit]
    dataPool$Cluster <- Sample.Def$MajorCellType[crit]
    resulT <-
      ggplot(data = dataPool, aes(
        x = x,
        y = y,
        color = Expression,
        text = paste("Cell ID: ", ID, "\nCluster: ", Cluster)
      )) + geom_point(alpha = 0.7, size =
                        dotSize) +
      scale_color_gradientn(
        colors = c(
          "#ffffcc",
          "#ffeda0",
          "#fed976",
          "#feb24c",
          "#fd8d3c",
          "#fc4e2a",
          "#e31a1c",
          "#bd0026",
          "#800026"
        )
      ) +
      theme_bw() + ggtitle(maiN) +
      theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        #legend.position = c(0.99, 0.99),
        #legend.justification = c(1, 1),
        legend.text = element_text(size = 10),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
      ) + labs(color = llaB)
    return(resulT)
  }





ui <- bootstrapPage(singleton(tags$head(
  tags$script(src = "message-handler.js"),
  tags$script('function getKey(){if(event.keyCode==13){$("#update").click();}}'),
  tags$script(src = "md5.js", type='text/javascript'),
  tags$script(src = "shinyBindings.js", type='text/javascript')
)),
tags$body(onload="setvalues()"),
tags$input(id = "ipid", class = "ipaddr", value='', type="text", style="display:none;"),

column(5,
       wellPanel(
         HTML(
           '<label for="gene">Gene</label>
           <div class="input-group shiny-input-container">
           <input id="gene" type="text" class="form-control" value="" placeholder="Type a gene name or gene ID." onkeypress="getKey();"/>
           <span class="input-group-btn">
           <button id="update" type="button" class="btn btn-primary action-button">Submit</button>
           </span>
           </div>'
         ),
         helpText("For batch or geneset plotting, separate genes with a comma(,)."),
         
         fluidRow(
           column(
             6,
             tags$hr(),
             checkboxGroupInput(
               "clusList",
               "Clusters",
               choices = clusterNames,
               selected = clusterNames
             ),
             
             checkboxGroupInput(
               "tissueList",
               "Tissue type",
               choices = tissueNames,
               selected = tissueNames
             ),
             
             checkboxGroupInput(
               "patientList",
               "Patient",
               choices = patientNames,
               selected = patientNames
             )
           ),
           
           column(
             6,
             tags$hr(),
             selectInput(
               "plotBy",
               "Generate Boxplot by",
               choices = plotParam,
               selected = 1
             ),
             HTML(
               '<div class="form-group shiny-input-container shiny-input-checkboxgroup" id="plotOptions">
               <label class="control-label" for="plotOptions">Plotting options</label>
               <div class="shiny-options-group">
               
               <div class="checkbox">
               <label>
               <input id="plotBee" type="checkbox" checked="checked"/>
               <span>Display scatter dots in Boxplot</span>
               </label>
               </div>
               
               <div class="checkbox">
               <label>
               <input id="plotColor" type="checkbox" checked="checked"/>
               <span>Display colored dots in Boxplot</span>
               </label>
               </div>
               
               <div class="checkbox">
               <label>
               <input id="plotLable" type="checkbox" checked="checked"/>
               <span>Display median labels in Boxplot</span>
               </label>
               </div>
               
               </div>
               </div>'
             ),
             sliderInput(
               "boxDotSize",
               "Boxplot dot size:",
               min = 0,
               max = 2,
               value = 0.4,
               step = 0.1
             ),
             sliderInput(
               "tsneDotSize",
               "tSNE plot dot size:",
               min = 0,
               max = 2,
               value = 1,
               step = 0.1
             )
             )
           
         )
         )),

column(7,
       tabsetPanel(
         tabPanel(
           "Boxplot",
           conditionalPanel(
             condition = "input.update",
             plotOutput("boxPlot",height=540),
             tags$br(),
             helpText(
               "The Boxplot shows the expression of selected gene in different populations. The thick horizontal line represents the median."
             ),
             helpText(
               "The Scatter dot is a visual presentation of data distribution. Every dot represents one cell."
             ),
             helpText(
               "For batch plotting, only the first gene will be displayed above. Please download the file for all plots."
             ),
             downloadButton('downloadBoxplot', 'Download')
           )
         ),
         tabPanel(
           "tSNE",
           conditionalPanel(
             condition = "input.update",
             plotOutput("tsnePlot",height = 540),
             tags$br(),
             helpText(
               "This plot shows the expression of selected gene in individual cells. Every dot represents one cell. The color of the dot represents the expression level."
             ),
             helpText(
               "For batch plotting, only the first gene will be displayed above. Please download the file for all plots."
             ),
             downloadButton('downloadTsne', 'Download')
           )
         ),
         tabPanel(
           "Geneset",
           conditionalPanel(
             condition = "input.update",
             plotOutput("tsneSetPlot",height = 540),
             verbatimTextOutput("geneset"),
             tags$br(),
             helpText(
               "This plot shows the expression of selected geneset in individual cells. Every dot represents one cell. The color of the dot represents the expression level."
             ),
             downloadButton('downloadTsneSet', 'Download')
           )
         ),
         tabPanel(
           "Interactive",
           conditionalPanel(
             condition = "input.update",
             plotlyOutput("tsneIntPlot",height=540),
             verbatimTextOutput("intset"),
             tags$br(),
             helpText(
               "This plot shows the expression of selected geneset in individual cells. Every dot represents one cell. The color of the dot represents the expression level."
             )
           )
         ),
         tabPanel("Cluster def.",
                  br(),
                  img(src = "tsneLegend.png", width = "100%"))
         ,
         tabPanel(
           "Dataset",
           conditionalPanel(
             condition = "input.update",
             tags$br(),
             DT::dataTableOutput("table"),
             tags$br(),
             helpText(
               "The Expression unit is normalized expression/log(TPM+1), depending on your selection. View Supplementery for more details."
             ),
             helpText(
               "For batch query, only the first gene will be displayed above. Please download the file for the complete dataset."
             ),
             downloadButton('downloadTable', 'Download')
           )
         )
       )))


server <- function(input, output, session) {
  getGene <-
    eventReactive(input$update,
                  toupper(strsplit(input$gene, " *, *")[[1]]),
                  ignoreNULL = TRUE,
                  ignoreInit = TRUE)
  
  
  observe({
    logFile <- file("genes.log", 'a')
    cat(paste(isolate(input$ipid), Sys.time(), "\n"),
        file = logFile,
        append = TRUE)
    cat(x = getGene(),
        file = logFile,
        append = TRUE)
    cat("\n----------\n\n", file = logFile , append = TRUE)
    close(logFile)
  })
  
  output$boxPlot <- renderPlot({
    p <- NULL
    progress <- shiny::Progress$new(style = "old")
    on.exit(progress$close())
    progress$set(message = "Checking parameters...", value = 0)
    gene <- getGene()[1]
    if (is.na(gene))
      session$sendCustomMessage(type = 'textMessage', message = "Empty/unrecognized gene name!")
    else if (is.null(isolate(input$clusList)) |
             is.null(isolate(input$tissueList)) |
             is.null(isolate(input$patientList)))
      session$sendCustomMessage(type = 'textMessage', message = "Empty dataset!")
    else if (!gene %in% names(Genes))
      session$sendCustomMessage(type = 'textMessage', message = "Gene not found in the dataset!")
    else {
      p <- autoBoxplot(
        gene,
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$plotBy),
        isolate(input$plotBee),
        isolate(input$boxDotSize),
        isolate(input$plotColor),
        isolate(input$plotLable),
        progress
      )
      if (length(getGene()) > 1)
        session$sendCustomMessage(type = 'textMessage', message = "Batch plot generated, press Download to download them as a .PDF file!")
    }
    progress$set(message = "Plot complete.", value = 1)
    p
  })
  output$downloadBoxplot <- downloadHandler(
    filename = function() {
      paste("Boxplot-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      if (isolate(input$gene) == "")
        return()
      if (is.null(isolate(input$clusList)) |
          is.null(isolate(input$tissueList)) |
          is.null(isolate(input$patientList)))
        return()
      
      pdf(file)
      for (m in getGene())
        if (m == "" |
            is.na(m)) {
          plot.new()
          grid.text('Empty/unrecognized gene name.')
        }
      else if (!m %in% names(Genes)) {
        plot.new()
        grid.text(paste0('Gene "', m, '" not found in the dataset.'))
      }
      else
        plot(
          autoBoxplot(
            m,
            isolate(input$clusList),
            isolate(input$tissueList),
            isolate(input$patientList),
            isolate(input$plotBy),
            isolate(input$plotBee),
            isolate(input$boxDotSize),
            isolate(input$plotColor),
            isolate(input$plotLable)
          )
        )
      dev.off()
    }
  )
  
  output$tsnePlot <- renderPlot({
    p <- NULL
    progress <- shiny::Progress$new(style = "old")
    on.exit(progress$close())
    gene <- getGene()[1]
    if (is.na(gene))
      session$sendCustomMessage(type = 'textMessage', message = "Empty/unrecognized gene name!")
    else if (is.null(isolate(input$clusList)) |
             is.null(isolate(input$tissueList)) |
             is.null(isolate(input$patientList)))
      session$sendCustomMessage(type = 'textMessage', message = "Empty dataset!")
    else if (!gene %in% names(Genes))
      session$sendCustomMessage(type = 'textMessage', message = "Gene not found in the dataset!")
    else {
      p <- autoTsne(
        gene,
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$tsneDotSize),
        progress
      )
      if (length(getGene()) > 1)
        session$sendCustomMessage(type = 'textMessage', message = "Batch plot generated, press Download to download them as a .PDF file!")
      progress$set(message = "Plot complete.", value = 1)
    }
    p
  })
  output$downloadTsne <- downloadHandler(
    filename = function() {
      paste("tSNE-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      if (isolate(input$gene) == "")
        return()
      if (is.null(isolate(input$clusList)) |
          is.null(isolate(input$tissueList)) |
          is.null(isolate(input$patientList)))
        return()
      pdf(file)
      for (m in getGene())
        if (m == "" |
            is.na(m)) {
          plot.new()
          grid.text('Empty/unrecognized gene name.')
        }
      else if (!m %in% names(genes)) {
        plot.new()
        grid.text(paste0('Gene "', m, '" not found in the dataset.'))
      }
      else
        plot(autoTsne(
          m,
          isolate(input$clusList),
          isolate(input$tissueList),
          isolate(input$patientList),
          isolate(input$tsneDotSize)
        ))
      dev.off()
    }
  )
  
  output$tsneSetPlot <- renderPlot({
    p <- NULL
    progress <- shiny::Progress$new(style = "old")
    on.exit(progress$close())
    geneset <- intersect(getGene(),names(Genes))
    output$geneset <- renderPrint(Genes[geneset])
    if (length(geneset) == 0)
      session$sendCustomMessage(type = 'textMessage', message = "Empty/unrecognized gene name!")
    else if (is.null(isolate(input$clusList)) |
             is.null(isolate(input$tissueList)) |
             is.null(isolate(input$patientList)))
      session$sendCustomMessage(type = 'textMessage', message = "Empty dataset!")
    else if (length(geneset) == 1) {
      p <- autoTsne(
        geneset[1],
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$tsneDotSize),
        progress
      )
      progress$set(message = "Plot complete.", value = 1)
    }
    else {
      p <- autoTsneSet(
        geneset,
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$tsneDotSize),
        progress
      )
      progress$set(message = "Plot complete.", value = 1)
    }
    p
  })
  output$downloadTsneSet <- downloadHandler(
    filename = function() {
      paste("Geneset-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      if (isolate(input$gene) == "")
        return()
      if (is.null(isolate(input$clusList)) |
          is.null(isolate(input$tissueList)) |
          is.null(isolate(input$patientList)))
        return()
      geneset <- intersect(getGene(),names(Genes))
      if (length(geneset) == 0)
        return()
      pdf(file)
      if (length(geneset) == 1)
        plot(
          autoTsne(
            geneset[1],
            isolate(input$clusList),
            isolate(input$tissueList),
            isolate(input$patientList),
            isolate(input$tsneDotSize),
            progress
          )
        )
      else
        plot(autoTsneSet(
          geneset,
          isolate(input$clusList),
          isolate(input$tissueList),
          isolate(input$patientList),
          isolate(input$tsneDotSize)
        ))
      dev.off()
    }
  )
  
  output$tsneIntPlot <- renderPlotly({
    p <- NULL
    progress <- shiny::Progress$new(style = "old")
    on.exit(progress$close())
    geneset <- intersect(getGene(),names(Genes))
    output$intset <- renderPrint(Genes[geneset])
    if (length(geneset) == 0)
      session$sendCustomMessage(type = 'textMessage', message = "Empty/unrecognized gene name!")
    else if (is.null(isolate(input$clusList)) |
             is.null(isolate(input$tissueList)) |
             is.null(isolate(input$patientList)))
      session$sendCustomMessage(type = 'textMessage', message = "Empty dataset!")
    else if (length(geneset) == 1) {
      p <- autoTsne(
        geneset[1],
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$tsneDotSize),
        progress
      )
      progress$set(message = "Plot complete.", value = 1)
    }
    else {
      p <- autoTsneSet(
        geneset,
        isolate(input$clusList),
        isolate(input$tissueList),
        isolate(input$patientList),
        isolate(input$tsneDotSize),
        progress
      )
      progress$set(message = "Plot complete.", value = 1)
    }
    if (is.null(p))
      ggplot() + theme_bw()
    else
      p
  })
  
  output$table <- DT::renderDataTable({
    if (length(getGene()) > 1)
      session$sendCustomMessage(type = 'textMessage', message = "Multiple-gene dataset generated, press Download to download!")
    gene <- getGene()[1]
    if (isolate(input$gene) == "" | is.na(gene))
      session$sendCustomMessage(type = 'textMessage', message = "Empty/unrecognized gene name!")
    else if (is.null(isolate(input$clusList)) |
             is.null(isolate(input$tissueList)) |
             is.null(isolate(input$patientList)))
      session$sendCustomMessage(type = 'textMessage', message = "Empty dataset!")
    else if (!gene %in% names(Genes))
      session$sendCustomMessage(type = 'textMessage', message = "Gene not found in the dataset!")
    else {
      crit <- Sample.Def$MajorCellType %in% isolate(input$clusList) &
        Sample.Def$tissueSource %in% isolate(input$tissueList) &
        Sample.Def$patient %in% isolate(input$patientList)
      Expression <- Sample.Exp[crit, gene]
      
      DT::datatable(
        colnames = c("Cell ID" = 1),
        class = "stripe table-condensed",
        filter = "bottom",
        cbind(Sample.Def[crit,], Expression),
        options = list(pageLength = 6)
      )
    }
  })
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("Dataset-", Sys.Date(), '.csv', sep = '')
    },
    content = function(file) {
      if (isolate(input$gene) == "")
        return()
      if (is.null(isolate(input$clusList)) |
          is.null(isolate(input$tissueList)) |
          is.null(isolate(input$patientList)))
        return()
      crit <- Sample.Def$MajorCellType %in% isolate(input$clusList) &
        Sample.Def$tissueSource %in% isolate(input$tissueList) &
        Sample.Def$patient %in% isolate(input$patientList)
      resulT <- Sample.Def[crit,]
      for (gene in getGene())
        if (gene %in% names(Genes)) {
          Expr <- Sample.Exp[crit, gene]
          resulT <- cbind(resulT, Expr)
          colnames(resulT)[ncol(resulT)] <- Genes[gene]
        }
      write.csv(cbind("Cell ID" = rownames(resulT), resulT), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)