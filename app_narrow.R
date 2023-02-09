###########################################################################
## Copyright 2018 Nantomics LLC                                          ##
##                                                                       ##
## Redistribution and use in source and binary forms, with or without    ##
## modification, are permitted for educational, research and non-profit  ##
## purposes, by non-profit institutions only provided that the following ##
## conditions are met:                                                   ##
## 1. Redistributions of source code must retain the above copyright     ##
## notice, this list of conditions and the following disclaimer.         ##
##                                                                       ##
## 2. Redistributions in binary form must reproduce the above copyright  ##
## notice, this list of conditions and the following disclaimer in the   ##
## documentation and/or other materials provided with the distribution.  ##
##                                                                       ##
## 3. Neither the name of the copyright holder nor the names of its      ##
## contributors may be used to endorse or promote products derived from  ##
## this software without specific prior written permission.              ##
##                                                                       ##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   ##
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ##
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ##
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT  ##
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,            ##
## INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,  ##
## BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS ##
## OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED    ##
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT           ##
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY ##
## WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           ##
## POSSIBILITY OF SUCH DAMAGE.                                           ##
###########################################################################

library(ggvis)
library(shiny)
library(Seurat)
library(Matrix)
library(feather)
library(tibble)

# metadata stored in "coldata"
#load("/home/shiny/data/kasp11_final_coldata.RData")
load("/data/main/kasp11_final_coldata.RData")

# log2 + 1 imputued data in magic.data
#load("/home/shiny/data/zinb-20_magic-t10_3UMI_100cells_1Pct_log.Rdata")
#magic.data <- read_feather("/data/main/magic.feather")
magic.data <- read_feather("/data/main/magic_selected.feather")
magic.data <- column_to_rownames(data.frame(magic.data, check.names=F))

dim(magic.data)

cluster.colors <- c('#008B00FF', '#7CAE00FF', '#FFA500FF', '#FFD700FF', '#00BFC4FF',
                    '#308ED2FF', '#7771B5FF', '#912CEEFF',  '#8B5A2BFF', '#CD853FFF','#FF5347FF')
names(cluster.colors) <- c("basal1", "basal2", "WNT1", "follicular", "channel",
                             "mitotic", "spinous", "granular", "mel1", "mel2", "immune")

tissue.colors <- cluster.colors[c(1,4,8,11)]
names(tissue.colors) <- levels(coldata$tissue)

sample.colors <- c(cluster.colors, "blue")
names(sample.colors) <- levels(coldata$sample)

#levels(coldata$cluster) <- c("basal1", "basal2", "WNT1", "follicular", "channel",
#                             "mitotic", "spinous", "granular", "mel2", "mel1", "immune")
#coldata$cluster <- factor(coldata$cluster, levels=c("immune", "mel2", "mel1", "granular", "spinous", "mitotic",
#                                                    "channel", "follicular", "WNT1", "basal2", "basal1"))

levels(coldata$cluster) <- c("basal1", "basal2", "WNT1", "follicular", "channel",
                             "mitotic", "spinous", "granular", "mel1", "mel2", "immune")
coldata$cluster  <- factor(coldata$cluster, levels=rev(levels(coldata$cluster)))


app.dat <- list(data=magic.data, data.info=coldata)

gene.list <- sort(row.names(app.dat$data))


makePlotData <- function(genes, expdata=app.dat, cells=rownames(app.dat$data.info)){

    curData <- data.frame(as.matrix(t(expdata$data[genes,])))
    curData <- cbind(curData[cells,], expdata$data.info[cells, c("sample", "tissue", "cluster")])

    return(curData)
}

scaleGeneVector <- function(gene.vec) {
    return((gene.vec/max(gene.vec))*100)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    tags$head(
        tags$style("#plot1{height:50vh !important;}"),
        tags$style("#plot2{height:50vh !important;}")
    ),

    # Application title
    titlePanel("SCARAB (Single Cell Analysis of RNA Browser)"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
	    selectizeInput('xvar', 'X Gene (Sorted Gene)', gene.list, selected="KRT10", options=list(maxOptions=7000)),
	    selectizeInput('yvar', 'Y Gene', gene.list, selected="KRT5", options=list(maxOptions=7000)),
            #checkboxGroupInput("facet", "Show", c("by.sample"), selected=c()),
	    checkboxGroupInput("selectTissue", "Tissue Select", c("foreskin", "trunk", "scalp", "psoriasis"), selected=c("foreskin", "scalp", "trunk")),
	    #radioButtons("color", "Color",
	    #             choices = list("cluster" = "cluster", "sample" = "sample", "tissue" = "tissue"), selected="cluster"),
	    #radioButtons("scale", "Scale",
	    #             choices = list("" = "cluster", "sample" = "sample", "tissue" = "tissue"), selected="cluster"),
            submitButton("Update View"),
	    tags$div(class="body", checked=NA,
	        tags$br(),
	        tags$p("SCARAB provides freely browsable transcript abundance data from the Cheng ", tags$i("et al"), ' manuscript and collaborative USCF-Nantomics project "Transcriptional programming of normal and inflamed human epidermis at single cell resolution".',
		    "Each point represents abundance of a given transcript in a single epidermal cell. Data is drawn from more than 80,000 epidermal cells profiled from 12 samples representing foreskin, trunk, scalp, and psoriatic epidermis.",
		    "2804 transcripts are visually browsable on this site based on abundance of 5 UMI in at least 100 cells or log fold change > .5 in at least one cluster.",
		    tags$br(), tags$br(), "Color-coding of cells is based on spectral clsutering of epidermal cells into lineages and differentiation states described in the manuscript.",
		    tags$br(), tags$br(), "The full underlying data will be accessible after manuscript acceptance at the European Genome-phenome Archive (EGA), which is hosted by EBI and CRG, under accession number EGAS00001002927.",
		    tags$br(), tags$br(), "Contact ", tags$a(href="mailto:Jeffrey.cheng@ucsf.edu", "Jeffrey.cheng@ucsf.edu"), " with questions related to the manuscript.",
		    tags$br(), "Contact ", tags$a(href="mailto:Andrew.Sedgewick@nantomics.com", "Andrew.Sedgewick@nantomics.com"), " with questions related to this site."
                )
	    )
        ),

        # Show a plot of the generated distribution
        mainPanel(
	    tabsetPanel(
	        id = 'type',
                tabPanel('xy-scatter', plotOutput("plot1")),
     	        tabPanel('versus sorted gene', plotOutput("plot2"), height="75%")
	    ),
	    height="75%"
        ), position="right"
    )
)
				      

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$plot1 <- renderPlot({
        cellFilter <- rownames(subset(coldata, tissue %in% input$selectTissue))
        curData <- makePlotData(c(input$xvar, input$yvar), cells=cellFilter)

        # set color palletes
        #if(input$color == "cluster"){
        #    plot.colors <- cluster.colors
        #} else if(input$color == "tissue"){
        #    plot.colors <- tissue.colors
        #} else {
        #    plot.colors <- sample.colors 
        #}
        #aes.color <- input$color

        aes.color <- "cluster"
        plot.colors <- cluster.colors

        p <- curData %>%
            ggplot(aes_string(input$xvar, input$yvar)) + geom_point(aes_string(color=aes.color), alpha=.2) +
	    theme_bw() + theme(text=element_text(color="grey30", face="bold", size=14)) + scale_color_manual(values = plot.colors) +
	    guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) + 
            labs(x=paste0(input$xvar, "   log2(imputed counts per 10k + 1)"), y=paste0(input$yvar, "   log2(imputed counts per 10k + 1)"))

        if(length(input$facet) > 0 && input$facet=="by.sample")
            p <- p + facet_wrap(~sample, ncol=4)
        p
      
    })

    output$plot2 <- renderPlot({
        cellFilter <- rownames(subset(coldata, tissue %in% input$selectTissue))
        curData <- makePlotData(c(input$xvar, input$yvar), cells=cellFilter)
        gene1.order <- order(curData[,input$xvar], decreasing = F)
        curData <- curData[gene1.order,]
        curData[,input$xvar] <- scaleGeneVector(curData[,input$xvar])
        curData[,input$yvar] <- scaleGeneVector(curData[,input$yvar])

        curData$index <- 1:nrow(curData)

        ## need block like this to calculate indexes if faceting
        #if(length(input$facet) > 0 && input$facet=="by.sample"){
        #    samps <- app.dat$data.info[curFilter, "sample", drop=F]
	#    samps <- samps[gene1.order,1,drop=F]
	#    gene1.inds <-  samps[,1] %>% levels %>%
	#                  lapply(function(x){ ou <- 1:sum(samps==x); names(ou) <- rownames(samps)[samps==x]; return(ou)}) %>%
	#		  unlist
	#		  
	#    plot.dat$index <- gene1.inds[rownames(samps)]
	#} else {
        #    plot.dat$index <- 1:sum(curFilter!=0)
        #}
	#samps %>% levels %>% lapply(function(x) 1:sum(samps==x)) %>% unlist

        #set color palettes
        #if(input$color == "cluster"){
        #    plot.colors <- cluster.colors
        #} else if(input$color == "tissue"){
        #    plot.colors <- tissue.colors
        #} else {
        #    plot.colors <- sample.colors 
        #}
        #aes.color <- input$color

	plot.colors <- cluster.colors
        aes.color <- "cluster"

        p <- ggplot(curData, aes_string('index', input$yvar)) + geom_point(aes_string(color=aes.color), alpha = .1,  size=1) +
            geom_point(color = "grey60", aes_string("index", input$xvar), size=1.5, alpha=.5) +
            theme_bw() + theme(text=element_text(color="grey30", face="bold", size=14)) + 
	    scale_color_manual(values = plot.colors) +
	    guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) + 
            labs(x=paste0("Cell index sorted relative to ", input$xvar, " expression"), y=paste0(input$xvar, " (Grey), ", input$yvar, "  expression percentile"))

        if(length(input$facet) > 0 && input$facet=="by.sample")
            p <- p + facet_wrap(~sample, ncol=4, scales="free")
        p
    })

#, height = function() {
    #      session$clientData$output_plot1_width
    #      })
}

# Run the application
shinyApp(ui = ui, server = server)