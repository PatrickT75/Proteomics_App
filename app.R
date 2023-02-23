library(shiny)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(rgl)
library(dplyr)
library(cluster)
library(rlist)
library(tibble)
library(ggvenn)
library(upsetjs)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(BiocManager)
library(msigdbr)
library(rlang)
library(tidyr)
library(enrichplot)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Proteomics"),
  h5("Version updated as of: 02-18-2023 11:22"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Expression File",
            multiple = TRUE,
            accept = c("text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")),
      
      fileInput("file2", "Metadata File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      tags$hr(),
      radioButtons("na_processing", "NA values: ",
                   choices = c("Omit", "Replace"),
                   selected = "Omit"),
      tags$hr(),
      sliderInput("sig_cutoff", "Significance",
                  0.001, 0.10, 0.05, 0.001),
      sliderInput("fc_cutoff", "Fold Change Cutoff",
                  1, 10, 2, 0.05),
      tags$hr(),
      radioButtons("adjustment", "Adjustment",
                   choices = c(None = "none",
                               BenjaminiHochberg = "BH",
                               Bonferroni = "bonferroni",
                               Holm = "holm"),
                   selected = "none"),
      
      actionButton("submit", "Submit"),
      
      actionButton("downloadoptions", "Download Report"),
      
      actionButton("demo", "Load Demo Data")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type="tab",
        tabPanel("Instructions", 
          h3("Please use the following format to upload your csv files."),
          h4("Expression file: "),
          tags$img(src='expression_format.png', height=150, width=600),
          tags$hr(),
          h4("Metadata file (no headers): "),
          tags$img(src='metadata_format.png', height=150, width=200),
          tags$hr(),
          h4("Before Uploading"),
          h5("Make sure that the sample names are spelled and ordered identically in the expression and metadata file.")
        ),
        
        tabPanel("Data", 
          h3("Data Summary"),
          radioButtons("disp", "Display",
                       choices = c(Data = "head",
                                   Statistics = "summary"),
                       selected = "head"),
          tableOutput("contents"),
          htmlOutput("length"),
          tags$hr(),
          h3("Protein Distribution Viewer"),
          textInput("protein_viewer_select", "Search Protein: "),
          actionButton("protein_viewer_run", "Run"),
          plotOutput("protein_viewer"),
          tags$hr(),
          h3("Grouping"),
          tableOutput("meta_check")
        ),
        
        tabPanel("Venn",
                 h3("Upset Plot"),
                 h5("Replicate proteins are counted once only."),
                 upsetjsOutput("upsetjs1"),
                 tags$hr(),
                 h3("Venn Diagram"),
                 h5("Replicate proteins are counted separately."),
                 checkboxGroupInput("venn_choices", "Select Up to 4 Groups: "),
                 actionButton("venn_run", "Run"),
                 plotOutput("venndiagram"),
                 downloadButton("venndiagram_download", "Download Figure"),
                 tags$hr(),
                 h3("Unique Protein Viewer"),
                 checkboxGroupInput("venn_unique", "Unique proteins from group(s) in the Venn Diagram:",
                                    choices = c(Group1 = "group1",
                                                Group2 = "group2"),
                                    selected = "group1"),
                 uiOutput("venn_descriptions"),
                 tableOutput("venndiagram_table")
                 ),
        
        tabPanel("Violin",
                 h3("Violin Plot"),
                 plotOutput("violin"),
                 downloadButton("violin_download", "Download Figure"),
                 tags$hr(),
                 h3("Aggregated Plot"),
                 plotOutput("violin_agg"),
                 downloadButton("violin_agg_download", "Download Figure")
        ),
        
        tabPanel("Correlation",
          h3("Correlation Heatmap"),
          radioButtons("corr_disable", "Show Values", choices=c("Enable", "Disable"), selected="Enable"),
          radioButtons("corr_type", "Compute correlation: ", choices=c(Pearson="pearson", Spearman="spearman"), selected="pearson"),
          plotOutput("correlation"),
          downloadButton("correlation_download", "Download Figure"),
          tags$hr(),
          h3("Scatterplot Matrix"),
          actionButton("scatterplot_run", "Run"),
          plotOutput("scatter"),
          downloadButton("scatter_download", "Download Figure")
        ),
        
        tabPanel("Volcano",
          h3("Parameters"),
          htmlOutput("volcano_parameters"),
          tags$hr(),
          h3("Volcano Plot"),
          h5("The log fold change reflects the ratio of expression in Group 2 compared to that in Group 1."),
          fixedRow(
            column(6, checkboxGroupInput("groups1", "Group 1: ")),
            column(6, checkboxGroupInput("groups2", "Group 2: "))
          ),
          tags$hr(),
          actionButton("volcano_run", "Run"),
          plotOutput("volcano"),
          downloadButton("volcano_download", "Download Figure"),
          htmlOutput("volcano_info"),
          tags$hr(),
          h3("Protein Table"),
          radioButtons("sort_volcano_list", "Sort By: ",
                       choices = c("P value",
                                  "Fold change",
                                  "Average Expression"),
                       selected = "P value"),
          tableOutput("volcano_protein_list"),
          
          radioButtons("volcanodownloadformat", "Format",
                       choices = c(List = "list",
                                   FlatFile = "info",
                                   Details = "details"),
                       selected = "details"),
          downloadButton("volcano_button", "Download Table"),
          
          tags$hr(),
          h3("MA Plot"),
          plotOutput("ma_plot"),
          tableOutput("ma_plot_info")
        ),
        
        tabPanel("Heatmap",
          h3("Expression Heatmap"),
          fixedRow(column(6, numericInput("heatmap_clusters", "Clusters: ", 
                    2, min=2, max=10)),
                   column(6, checkboxGroupInput("groups_heatmap", "Select Groups: "))
                   ),
          actionButton("heatmap_run", "Run"),
          fixedRow(plotOutput("heatmap")),
          uiOutput("heatmap_cluster_numbers"),
          tags$hr(),
          uiOutput("heatmap_info"),
          downloadButton("heatmap_download", "Download Figure"),
          tags$hr(),
          h3("Number of Clusters"),
          h5("Optimal number of clusters is the maximum silhouette score represented below."),
          plotOutput("optimal_hm"),
          downloadButton("heatmap_button", "Download Clusters"),
        ),
        
        tabPanel("PCA", 
          h3("Principal Component Analysis"),
          h5("Choose the number of proteins to include in the PCA, sorted by the highest coefficient of variation."),
          numericInput("pca_topcv", "Number of Proteins:", 100,
                       min = 1, max = 2000),
          h5("Plot the PCs as (x,y,z) in order of importance. The first 5 PCs may be plotted."),
          fixedRow(
          column(4, radioButtons("pc_x", "X:",
                                 choices = c(PC1 = "1", 
                                             PC2 = "2", 
                                             PC3 = "3", 
                                             PC4 = "4", 
                                             PC5 = "5"),
                                 selected = "1")),
          column(4, radioButtons("pc_y", "Y:",
                                 choices = c(PC2 = "2", 
                                             PC3 = "3", 
                                             PC4 = "4", 
                                             PC5 = "5"),
                                 selected = "2")),
          column(4, radioButtons("pc_z", "Z:",
                                 choices = c(PC3 = "3", 
                                             PC4 = "4", 
                                             PC5 = "5"),
                                 selected = "3"))
          ),
          actionButton("pca_run", "Run"),
          plotOutput("pca"),
          downloadButton("pca_download", "Download Figure"),
          tags$hr(),
          h3("3-D Principal Component Analysis"),
          column(12, rglwidgetOutput("pca3d", width="700px", height="512px")),
          actionButton("pca3d_snapshot", "Take snapshot"),
          tags$hr(),
          h3("Scree Plot"),
          plotOutput("pca_importance"),
          downloadButton("pca_importance_download", "Download Figure"),
          tags$hr(),
          h3("Top and Bottom Loadings"),
          fixedRow(column(6, radioButtons("pc_load", "PC:",
                                 choices = c(PC1 = "PC1", 
                                             PC2 = "PC2", 
                                             PC3 = "PC3", 
                                             PC4 = "PC4", 
                                             PC5 = "PC5"),
                                 selected = "PC1"))),
          
          fixedRow(column(6, div(style='padding:5px'), tableOutput("pca_loadings_top")),
                   column(6, div(style='padding:5px'), tableOutput("pca_loadings_bot"))),
          downloadButton("loadings_button", "Download Loadings"),
          plotOutput("pca_loadings_plot"),
          downloadButton("pca_loadings_plot_download", "Download Figure"),
          tags$hr(),
          fixedRow(
          column(6, radioButtons("select_loading", "Plot Distribution of: ", choices=c("Protein 1", "Protein 2")),
                    textInput("select_loading_text", "Search: "),
                    textInput("normalize", "Normalize to (if present): "),
                    actionButton("plot_loadings_scatter", "Plot")),
          ),
          plotOutput("loadings_scatter")
        ),
        
        tabPanel("Pathway",
                 h3("Gene Ontology"),
                 radioButtons("go_proteinlist", "Select protein list: ",
                              choices = c("Heatmap",
                                          "Top CV",
                                          "Upload File"),
                              selected = "Top CV"),
                 
                 tags$hr(),
                 conditionalPanel(condition = "input.go_proteinlist == 'Heatmap'",
                                  h5("Use proteins from the heatmap tab that have been organized by cluster."),
                                  checkboxGroupInput("go_heatmap_clusters_choices", "Select Clusters:")),
                 conditionalPanel(condition = "input.go_proteinlist == 'Top CV'",
                                  numericInput("go_topcv", "Top proteins: ",
                                                1000, min=100, max=5000)),
                 conditionalPanel(condition = "input.go_proteinlist == 'Upload File'",
                                  h5("Upload a csv file with only a list of proteins, no headers."),
                                  fileInput("go_file", "Select List to Upload: ",
                                            multiple = TRUE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv"))),
                 tags$hr(),
                 radioButtons("GO_ont", "Select Ontology: ",
                              choices = c(Biological_Process = "BP",
                                          Cellular_Compartment = "CC",
                                          Molecular_Function = "MF"),
                              selected = "CC"),
                 actionButton("go_calculate", "Run"),
                 plotOutput("go_dotplot"),
                 downloadButton("go_dotplot_download", "Donwload Figure"),
                 uiOutput("go_info"),
                 tags$hr(),
                 selectInput("go_genelist", "See genes involved in: ",
                             choices=c(1,2,3)),
                 checkboxGroupInput("go_heatmap_groups", "Show groups: "),
                 actionButton("go_heatmap_run", "View"),
                 plotOutput("go_new_heatmap"),
                 downloadButton("go_heatmap_download", "Download Heatmap")
                 ),
        
        tabPanel("GSEA",
                 h3("Gene Set Enrichment Analysis"),
                 radioButtons("gsea_proteinlist", "Select protein list: ",
                              choices = c("All",
                                          "Heatmap",
                                          "Upload File"),
                              selected = "All"),
                 tags$hr(),
                 conditionalPanel(condition = "input.gsea_proteinlist == 'Heatmap'",
                                  h5("Use proteins from the heatmap tab that have been organized by cluster."),
                                  checkboxGroupInput("gsea_heatmap_clusters_choices", "Select Clusters:",
                                    choices = c(1,2), selected = c(1,2))),
                 conditionalPanel(condition = "input.gsea_proteinlist == 'Upload File'",
                                  h5("Upload a csv with no headers: first column is protein name, second column is log FC"),
                                  fileInput("gsea_file", "Select List to Upload: ",
                                            multiple = TRUE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv"))),
                 tags$hr(),
                 radioButtons("gsea_set", "Select Molecular Signature: ",
                              choices = c(Hallmark = "H",
                                          Regulatory = "C3",
                                          Oncogenic = "C6",
                                          Immunologic = "C7"),
                              selected = "H"),
                 actionButton("gsea_calculate", "Run"),
                 tags$hr(),
                 selectInput("gsea_genelist", "See genes involved in:", 
                             choices = c(1, 2),
                             selected = 1),
                 actionButton("gsea_heatmap_run", "View"),
                 tags$hr(),
                 h3("Running Enrichment Score"),
                 plotOutput("gsea_running"),
                 downloadButton("gsea_running_download", "Download Figure"),
                 tags$hr(),
                 h3("Heatmap of Selected Proteins"),
                 checkboxGroupInput("gsea_heatmap_groups", "Show groups: "),
                 plotOutput("gsea_new_heatmap"),
                 downloadButton("gsea_heatmap_download", "Download Figure")
          ),
        selected = "Data"
      )
    )
    
  )
)


# Define server logic to read selected file ----
server <- function(input, output, session) {
  
   while (!is.null(dev.list()))  dev.off()
   
   ##################### INPUT########################
   # update checkbox groups for subsequent analyses
   observe({
    choice <- names()
    select1 <- choice[1]
    select2 <- choice[2]
    
    updateCheckboxGroupInput(session = session,
                             inputId = "groups1",
                             choices = choice,
                             selected = select1)
    updateCheckboxGroupInput(session = session,
                             inputId = "groups2",
                             choices = choice,
                             selected = select2)
    updateCheckboxGroupInput(session = session,
                             inputId = "venn_choices",
                             choices = choice,
                             selected = c(select1, select2))
    updateCheckboxGroupInput(session = session,
                             inputId = "groups_heatmap",
                             choices = choice,
                             selected = choice)
    updateCheckboxGroupInput(session = session,
                             inputId = "go_heatmap_groups",
                             choices = choice,
                             selected = choice)
    updateCheckboxGroupInput(session = session,
                             inputId = "gsea_heatmap_groups",
                             choices = choice,
                             selected = choice)
    
  })
   
   # total number of samples
   observations <- eventReactive(input$submit, {
     return(base::nrow(metadf()))
   })
    
   # show grouping of samples/replicates
   metadf <- reactive({
     req(input$file2)
     df <- read.csv(input$file2$datapath,
                    header = FALSE,
                    sep = ",",
                    quote = '"',
                    encoding = "UTF-8")
     if(is.null(input$file1) && input$demo > 0) {
       df <- read.csv(demoDataMetadata,
                      header = TRUE,
                      sep = ",",
                      quote = '"',
                      encoding = "UTF-8")
     }
     
     df$V2 <- as.factor(df$V2)
     return(df)
   })
   
   # shows how many replicates are in each group
   metadf_summary <- reactive({
     sum <- base::as.data.frame(summary(metadf()$V2))
     base::colnames(sum) <- "Replicates"
     return(sum)
   })
   
   ngroup <- reactive({
     l <- length(levels(metadf()$V2))
     return(l)
   })
   
   # group names
   names <- eventReactive(input$submit, {
     return(rownames(metadf_summary()))
   })
   
   # for each group, how many replicates are there?
   times <- reactive({
     times <- as.vector(metadf_summary())
     return(times$Replicates)
   })
   
   # expression of each replicate, ordered the same as the grouping. NA can be omitted.
   df <- reactive({
     req(input$file1)
     df <- read.csv(input$file1$datapath,
                    header = TRUE,
                    sep = ",",
                    quote = '"',
                    encoding = "UTF-8")
     
     if(is.null(input$file1) && input$demo > 0) {
     df <- read.csv(demoDataExpression,
                      header = TRUE,
                      sep = ",",
                      quote = '"',
                      encoding = "UTF-8")
     }

     if (input$na_processing == "Omit") {df <- na.omit(df)}
     if (input$na_processing == "Replace") {df[is.na(df)] <- 1}

     Protein <- df[,1]
     df <- df[,2:ncol(df)]
     df <- df[,order(metadf()$V2)]
     df <- cbind.data.frame(Protein,df)
     
     return(df)
   })
   
   # matrix only - no names of proteins
   dfres <- reactive({
     obs <- observations()
     dfr <- df()[2:(obs+1)]
     return(dfr)
   })
  
   # inverted dataframe - samples are in rows. df_trans()$group gives the grouping of each sample.
   df_trans <- reactive({
     df_t <- base::as.data.frame(base::t(dfres()))
     df_t$group <- rep(names(), times())
     return(df_t)
   })
   
   
   ################# DATA AND QUALITY CHECK #################
   # show contents for check
   contents_table <- eventReactive(input$submit, {
     means <- apply(dfres(), 1, mean)
     max <- apply(dfres(), 1, max)
     min <- apply(dfres(), 1, min)
     stdev <- apply(dfres(), 1, sd)
     pcv <- stdev/means*100
     
     summary <- cbind.data.frame(df()[,1], means, max, min, stdev, pcv)
     base::colnames(summary) <- c("Gene.Name", "Average", "Maximum", "Minimum", "SD", "%CV")
     
     p <- df()[1:(observations()+1)]
     
     if(input$disp == "head") {
       p <- utils::head(p, n=20)
     }
     if(input$disp == "summary") {
       p <- utils::head(summary, n=20)
     }
     
     return(p)
   })
   
   output$contents <- renderTable({
     contents_table()
   })
   
   # description for contents table
   output$length <- renderUI({
     len <- base::nrow(dfres())
     
     str1 <- base::paste("Number of groups: ", ngroup())
     str2 <- base::paste("Number of samples: ", observations())
     str3 <- base::paste("Number of proteins included: ", len)
     HTML(base::paste(str1, str2, str3, sep = '<br/>'))
     
   })
   
   output$meta_check <- renderTable({
     m <- metadf()
     m <- m[order(m$V2),]
     base::colnames(m) <- c("Sample", "Group")
     return(m)
   })
   
   # view distribution of individual protein
   protdist_viewer <- eventReactive(input$protein_viewer_run, {
     table <- base::subset(df(), df()[,1] == input$protein_viewer_select)
     table <- table[,2:ncol(table)]
     table <- melt(table)
     table$group <- rep(names(), times())
     p <- ggplot(table, aes(x=group, y=log(value), fill=group)) +
      geom_boxplot() + geom_point(position = position_jitter(seed = 1, width = 0.2))
     return(p)
   })
   
   output$protein_viewer <- renderPlot({
     return(protdist_viewer())
   })
   
   ################# VENN DIAGRAM #########################
   # update venn diagram unique protein viewer to only include groups in the venn diagram
   observeEvent(input$venn_run, {
     choice <- input$venn_choices
     updateCheckboxGroupInput(session = session,
                              inputId = "venn_unique",
                              choices = choice,
                              selected = choice[1])
   })
   
   output$upsetjs1 <- renderUpsetjs({
     req(input$file1)
     df <- read.csv(input$file1$datapath,
                    header = TRUE,
                    sep = ",",
                    quote = '"',
                    encoding = "UTF-8")
     protein_list <- list()

     count <- 1
     for (i in 1:ngroup())
     {
       test <- !(apply(df[,(count+1):(count+times()[i])], 1, anyNA))
       proteins <- df[,1][test]
       count <- count+times()[i]
       protein_list <- list.append(protein_list, proteins)
     }

     names(protein_list) <- names()


     upsetjs() %>% fromList(protein_list) %>% interactiveChart()
   })
   
   # table output: first column = protein names, rest of the columns = Present(T)/Absent(F) for each group
   venn_diagram_table <- eventReactive(input$venn_run, {
     req(input$file1)
     df <- read.csv(input$file1$datapath,
                    header = TRUE,
                    sep = ",",
                    quote = '"',
                    encoding = "UTF-8")
     # df <- df[rowSums(df[,2:ncol(df)])>0,]
     
     # subset dataframe into smaller dataframes
     mini_dfs <- list()
     for (i in 1:length(names())){
       s <- base::subset(metadf(), metadf()$V2 == names()[i])
       mini_dfs <- list.append(mini_dfs, df[,s$V1])
     }
     
     # each smaller dataframe is examined for NAs (T/F) in any of the replicates
     na_s <- list()
     for (i in 1:length(mini_dfs)){
       n <- !(apply(mini_dfs[[i]], 1, anyNA))
       na_s <- list.append(na_s, n)
     }
     names(na_s) <- names()
     
     na_s <- as.data.frame(do.call(cbind, na_s))
     d <- cbind.data.frame(df[,1], na_s)
     base::colnames(d) <- c("Protein", names())
     return(d)
   })
   
   venn_diagram <- eventReactive(input$venn_run, {
     if(length(input$venn_choices) == 2) {
       p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2]), show_percentage=FALSE, stroke_size=0.5)
       }
     else if(length(input$venn_choices) == 3) {
       p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2], input$venn_choices[3]), show_percentage=FALSE, stroke_size=0.5)
       }
     else if(length(input$venn_choices) == 4) {
       p <- ggvenn(venn_diagram_table(), c(input$venn_choices[1], input$venn_choices[2], input$venn_choices[3], input$venn_choices[4]), show_percentage=FALSE, stroke_size=0.5)
       }
     return(p)
   })
   
   output$venndiagram <- renderPlot({
     return(venn_diagram())
   })
   
   output$venndiagram_table_test <- renderTable({
     v <- venn_diagram_table()
     
     show <- rep(TRUE, length(names()))
     for(i in 1:length(names()))
     {if (names()[i] %in% input$venn_unique) {show[i] <- TRUE}
       else {show[i] <- FALSE}
     }
     
     for(i in 1:(length(show)))
     {
       v <- v[v[,i+1] == show[i], ]
     }
     return(v[,1])
   })
   
   # this table shows unique proteins from comparisons ONLY if they are in the venn diagram (disregards unselected groups)
   filtered_venndiagram_table <- reactive({
     v <- venn_diagram_table()
     v <- v[,c("Protein", input$venn_choices)]
     
     #replaced names with venn_choices
     show <- rep(TRUE, length(input$venn_choices))
     for(i in 1:length(input$venn_choices))
     {if (input$venn_choices[i] %in% input$venn_unique) {show[i] <- TRUE}
       else {show[i] <- FALSE}
     }
     
     for(i in 1:(length(show)))
     {
       v <- v[v[,i+1] == show[i], ]
     }
     return(v[,1])
   })
   
   output$venndiagram_table <- renderTable({
     return(filtered_venndiagram_table())
   })
   
   output$venn_descriptions <- renderUI({
     return(base::paste("There are ", length(filtered_venndiagram_table()), "unique proteins in the group(s): ", base::paste(input$venn_unique, collapse=", ", ".", sep = "")))
   })
   
   
   ######################## VIOLIN PLOT ############################
   # simple formula for transparency of points in violin plot
   alpha <- reactive({
     if(base::nrow(dfres()) < 1000)
     {return(0.2)}
     if(base::nrow(dfres()) < 2500)
     {return(0.1)}
     if(base::nrow(dfres()) < 4000)
     {return(0.03)}
     if(base::nrow(dfres()) >= 4000)
     {return(0.01)}
   })
   
   violin <- eventReactive(input$submit, {
     len <- base::nrow(dfres())
     dfres_melt <- melt(dfres())
     p <- ggplot(dfres_melt, aes(x=variable, y=log(value), fill=variable)) +
       geom_violin() + geom_boxplot() + geom_point(position = position_jitter(seed = 1, width = 0.2), alpha=alpha()) +
       theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),)
     return(p)
   })
   
   output$violin <- renderPlot({
     return(violin())
   })
   
   violin_agg <- eventReactive(input$submit, {
     len <- base::nrow(dfres())
     dfres_melt <- melt(dfres())
     
     grouping <- rep(names(), times())
     
     dfr <- cbind.data.frame(dfres_melt, grouping)
     
     p <- ggplot(dfr, aes(x=grouping, y=log(value), fill=grouping)) +
       geom_violin() + geom_boxplot() + geom_point(position = position_jitter(seed = 1, width = 0.2), alpha=alpha())
     return(p)
   })
   
   output$violin_agg <- renderPlot({
     return(violin_agg())
   })
   
   
   ################## CORRELATION AND SCATTER PLOT #####################
   # Pearson or Spearman.
   corrmap <- reactive({
     cormat <- cor(dfres(), method=input$corr_type)
     
     get_lower_tri<-function(cormat){
       cormat[upper.tri(cormat)] <- NA
       return(cormat)
     }
     
     get_upper_tri<- function(cormat){
       cormat[lower.tri(cormat)] <- NA
       return(cormat)
     }
     
     # Get lower triangular matrix
     lower_tri<- get_lower_tri(cormat)
     
     melted_cormat <- reshape2::melt(lower_tri, na.rm=TRUE)
     
     crc_p<-ggplot(data= melted_cormat, aes(x = Var2, y=reorder(Var1, dplyr::desc(Var1)), fill=value)) +
       geom_tile(color="white") +
       theme_minimal() +
       theme(axis.text.x=element_text(angle=45, vjust=1,
                                      
                                      size=12, hjust=1)) +
       scale_x_discrete(labels=rep(names(), times())) +
       scale_y_discrete(labels=rev(rep(names(), times()))) +
       coord_fixed() +
       theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(),
         legend.justification = c(1, 0),
         legend.position = c(0.8, 0.7),
         legend.direction = "horizontal") +
       guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                    title.position = "top", title.hjust = 0.5))
     
     if(input$corr_disable == "Enable") {crc_p <- crc_p + geom_text(aes(label=round(value,2)), color = "black", size = 4)}
     if(input$corr_type == "pearson") {crc_p <-  crc_p + scale_fill_gradient2(low="blue", high="red", mid="white",
                                                              midpoint=((min(unlist(melted_cormat$value)))+1)/2, limit=c(min(unlist(melted_cormat$value)),1.0), space="Lab",
                                                              name="Pearson\nCorrelation")
     }
     else if(input$corr_type == "spearman") {crc_p <- crc_p + scale_fill_gradient2(low="blue", high="red", mid="white",
                                                              midpoint=((min(unlist(melted_cormat$value)))+1)/2, limit=c(min(unlist(melted_cormat$value)),1.0), space="Lab",
                                                              name="Spearman\nCorrelation")
     }
     return(crc_p)
   })
   
   output$correlation <- renderPlot({
     print(corrmap())
   })
   
   # Scatterplot for every replicate
   scattermat <- eventReactive(input$scatterplot_run, {
     lower.panel<-function(x, y){
       points(x,y, pch = 20)
     }
     
     p <- pairs(dfres(), lower.panel = lower.panel, 
                upper.panel = NULL)
     
     return(p)
   })
   
   output$scatter <- renderPlot({
     return(scattermat())
   })
   
   
   ####################### VOLCANO PLOT ##########################
   # dataframe of all replicates in group 1
   v1 <- eventReactive(input$volcano_run, {
     df_1 <- base::subset(df_trans(), df_trans()$group %in% input$groups1)
     v1 <- t(df_1[1:(base::nrow(dfres()))])
     return(v1)
   })
   
   # dataframe of all replicates in group 2
   v2 <- eventReactive(input$volcano_run, {
     df_2 <- base::subset(df_trans(), df_trans()$group %in% input$groups2)
     v2 <- t(df_2[1:(base::nrow(dfres()))])
     return(v2)
   })
   
   # calculates log2fc between v1 and v2
   Log2FC <- eventReactive(input$volcano_run, {
     n1 <- apply(v1(), 1, mean)
     n2 <- apply(v2(), 1, mean)
     
     fc <- n2/n1
     Log2FC <- log2(fc)
     
     return(Log2FC)
   })
   
   # for each protein, calculates the t-test p-value between v1 and v2
   df_pval <- eventReactive(input$volcano_run, {
     l_1 <- ncol(v1())
     l_2 <- ncol(v2())
     
     vfull <- cbind(v1(), v2())
     
     t_testing2 <- function(x, l1, l2) {
       x_1 = x[1:l1]
       x_2 = x[(l1+1):(l1+l2)]
       if(length(unique(x_1)) == 1 | length(unique(x_2)) == 1)
       {x_1[1] <- x_1[1]+0.001
        x_2[1] <- x_2[1]+0.001}
       result <- t.test(x_1, x_2)$p.value
       return(result)
     }
     
     Pval <- base::apply(vfull, 1, t_testing2, l1=l_1, l2=l_2)
     
     if(input$adjustment != 'none') {
       Pval <- p.adjust(Pval, input$adjustment) 
     }
     
     LogPval <- -1*log10(Pval) 
     
     tt_log <- cbind(Pval, LogPval)
     
     return(tt_log)
   })
   
   # combines p-value and log2fc into one dataframe
   df_tt <- eventReactive(input$volcano_run, {
     df_tt <- base::cbind(dfres(), Log2FC(), df_pval())
     df_tt_test <- df_tt[order(df_tt$Pval),]
     print(head(df_tt_test))
     return(df_tt)
   })
   
   # combines p-value, log2fc, and significance into one dataframe
   df_sig <- eventReactive(input$volcano_run, {
     fc <- df_tt()[,(base::ncol(df_tt())-2)]
     logpval <- df_pval()[,2]
     bind <- cbind(logpval, fc)
     sig <- bind[,1] >= -1*log10(input$sig_cutoff) & abs(bind[,2]) >= log2(input$fc_cutoff) 
     df_sig <- cbind(df_tt(), sig)
     return(df_sig)
   })
   
   output$volcano_parameters <- renderUI({
     str1 <- base::paste("Significance: ", input$sig_cutoff)
     str2 <- base::paste("Adjustment method: ", input$adjustment)
     str3 <- base::paste("Fold change cutoff: ", input$fc_cutoff)
     HTML(base::paste(str1, str2, str3, sep = '<br/>'))
   })
   
   volcano <- eventReactive(input$volcano_run, {
     df_volcano <- df_sig()
     g1 <- base::paste(input$groups1, collapse='/')
     g2 <- base::paste(input$groups2, collapse='/')
     title <- base::paste(g1, "vs", g2, sep=' ')
     p <- ggplot(df_volcano) +
       geom_point(mapping=aes(x=Log2FC(),y=LogPval,col=sig))+labs(x = "log2(fold change)", y = "-log10(adjusted P-value)")+ggtitle(title)
     return(p)
   })
   
   output$volcano <- renderPlot({
     return(volcano())
   })
   
   # summary of results from volcano plot
   volcano_info <- eventReactive(input$volcano_run, {
     greater <- base::nrow(base::subset(df_sig(), df_sig()$sig == TRUE & df_sig()[,(ncol(df_sig())-3)] > 0))
     smaller <- base::nrow(base::subset(df_sig(), df_sig()$sig == TRUE & df_sig()[,(ncol(df_sig())-3)] < 0))
     str1 <- base::paste("Positive LogFC: ", greater)
     str2 <- base::paste("Negative LogFC: ", smaller)
     str3 <- base::paste("A positive LogFC corresponds to higher expression in ", base::paste(input$groups2, collapse=" / "), " compared to ", base::paste(input$groups1, collapse=" / "), ".", sep="")
     HTML(base::paste(str1, str2, str3, sep = '<br/>'))
   })
   
   output$volcano_info <- renderUI({
     return(volcano_info())
   })
   
   # means of all groups in the entire dataset, to show in volcano list
   means_test_volcano <- eventReactive(input$volcano_run, {
     means_list <- list()
     
     # apply means for each group
     for (i in 1:length(names())) {
       s <- base::subset(df_trans(), df_trans()$group == names()[i])
       s <- t(s[1:base::nrow(dfres())])
       means <- apply(s, 1, mean)
       means_list <- list.append(means_list, means)
     }
     
     means_list <- t(base::do.call(base::rbind, means_list))
     return(means_list)
   })
   
   # list of significant proteins
   volcano_list <- reactive({
     protein <- df()[,1]
     t <- base::cbind(protein, df_sig())
     t <- base::subset(t, t$sig == TRUE)
     t <- t[,-c(2:(observations()+1))]
     t <- t[,-ncol(t)]
     t <- cbind.data.frame(t, MA_plot_info())
     if(input$sort_volcano_list == "Fold change") {
       t <- t[sort(abs(t[,2]), decreasing=TRUE, index.return=TRUE)[[2]],]
     }
     else if (input$sort_volcano_list == "P value") {
       t <- t[base::order(t$LogPval, decreasing=TRUE),]
     }
     else if (input$sort_volcano_list == "Average Expression") {
       t <- t[base::order(t$Expression, decreasing=TRUE),]
     }
     return(t[-c(5,6)])
   })
   
   output$volcano_protein_list <- renderTable({
     return(utils::head(volcano_list(), n=20))
   }, digits=4)
   
   # x-axis: log expression, y-axis: log2FC
   MA_plot <- eventReactive(input$volcano_run, {
     dfs <- df_sig()
     means <- apply(dfs[,1:observations()], 1, mean)
     dfs <- cbind(dfs, means)
     return(dfs)
   })
   
   output$ma_plot <- renderPlot({
     dfs <- MA_plot()
     dfs$means <- log(dfs$means)
     p <- ggplot(dfs) +
       geom_point(mapping=aes(x=means,y=Log2FC(),col=sig))+labs(x = "log expression", y = "log2(fold change)")
     return(p)
   })
   
   # list of significant proteins ordered by highest average expression
   MA_plot_info <- eventReactive(input$volcano_run, {
     names <- as.vector(names())
     for(i in 1:length(names)){
       names[i] <- paste("Average ", names[i], sep="")
     }

     dfs <- MA_plot()
     sigs <- dfs[dfs$sig == TRUE,]
     df_selected <- df()[rownames(sigs),1]
     df_selected <- as.data.frame(df_selected)
     
     means_selected <- means_test_volcano()[rownames(sigs),]
     print(is.vector(means_selected))
     
     if(is.vector(means_selected)){
       means_selected <- t(means_selected)
     }
     
     t <- cbind.data.frame(df_selected, sigs$means)
     
     t <- cbind.data.frame(t, means_selected)
     

     
     colnames(t) <- c("Protein", "Expression", names)
     print(head(t))
     return(t)
   })
   
   output$ma_plot_info <- renderTable({
     MA_pi <- MA_plot_info()[order(-MA_plot_info()$Expression),]
     MA_pi$Expression <- log(MA_pi$Expression)
     return(utils::head(MA_pi[c(1,2)], n=20))
   })
   

   
   ############################ HEATMAP ###############################
   # collects means for groups that are selected
   means_test <- reactive({
     means_list <- list()
     
     for (i in 1:length(input$groups_heatmap)) {
       s <- base::subset(df_trans(), df_trans()$group == input$groups_heatmap[i])
       s <- t(s[1:base::nrow(dfres())])
       means <- apply(s, 1, mean)
       means_list <- list.append(means_list, means)
     }
     
     means_list <- t(base::do.call(base::rbind, means_list))
     return(means_list)
   })
   
   # chooses the maximum mean
   v_max_test <- reactive({
     max_index <- apply(means_test(), 1, base::which.max)
     return(as.vector(max_index))
   })
   
   # chooses the minimum mean
   v_min_test <- reactive({
     min_index <- apply(means_test(), 1, base::which.min)
     return(as.vector(min_index))
   })
   
   # calculates log2fc between maximum average and minimum average groups
   Log2FC_test <- reactive({
     max <- apply(means_test(), 1, max)
     min <- apply(means_test(), 1, min)
     
     Log2FC <- log2(max) - log2(min)
     return(Log2FC)
   })
   
   # calculates t-test p-value between maximum average and minimum average groups, returns in dataframe
   df_pval_test <- reactive({
     ## ALGORITHM:
     # Keep rolling count (int i)
     # Take index of max, use that to ask which group that refers to
     # Subset df_trans by that group name
     # Take column i of df_trans to get the vector of values
     # Repeat for min
     # T-test values of max vs values of min
     
     Pval <- c()
     t_testing3 <- function(i, l1, l2) {
       maxs <- df_trans()[df_trans()$group == input$groups_heatmap[l1[i]], i]
       mins <- df_trans()[df_trans()$group == input$groups_heatmap[l2[i]], i]
       
       # t-test in R doesn't accept all the same values for each group
       if (length(unique(maxs) == 1) | length(unique(mins)) == 1)
       {maxs[1] <- maxs[1] + 0.001
       mins[1] <- mins[1] + 0.001}
       result <- t.test(maxs, mins)$p.value
       return(result)
     }
     
     # ALTERNATE MEASURE IF MAXS/MINS GROUPS ARE BOTH IDENTICAL
     
     # if (length(unique(maxs) == 1) & length(unique(mins)) == 1)
     # {
         # ??
     # return(result)
     # }
     
     for (i in 1:length(v_max_test())) {
       Pval <- append(Pval, t_testing3(i, v_max_test(), v_min_test()))
     }
     
     if(input$adjustment != 'none') {
       Pval <- p.adjust(Pval, input$adjustment) 
     }
     
     LogPval <- -1*log10(Pval) 
     tt_log <- cbind(Pval, LogPval)
     
     return(tt_log)
     
   })
   
   # expression, log2fc, and p-value in one dataframe
   df_tt_test <- reactive({
     df_tt <- cbind(dfres(), Log2FC_test(), df_pval_test())
     return(df_tt)
   })
   
   # expression, log2fc, p-value, and significance in one dataframe
   df_sig_test <- reactive({
     fc <- df_tt_test()[,(ncol(df_tt_test())-2)]
     logpval <- df_pval_test()[,2]
     bind <- cbind(logpval, fc)
     sig <- bind[,1] >= -1*log10(input$sig_cutoff) & abs(bind[,2]) >= log2(input$fc_cutoff) 
     df_sig <- cbind(df_tt_test(), sig)
     return(df_sig)
   })
   
   # list of significant proteins - but from all pairwise comparisons
   volcano_list_test <- eventReactive(input$heatmap_run, {
     protein <- df()[,1]
     t <- base::cbind(protein, df_sig_test())
     t <- base::subset(t, t$sig == TRUE)
     t <- t[,-c(2:(observations()+1))]
     t <- t[,-base::ncol(t)]
     return(t)
   })
   
   # index of groups selected to be shown on the heatmap
   heatmap_index_test <- reactive({
     v <- c()
     v <- append(v, which(df_trans()$group %in% input$groups_heatmap))
     return(v)
   })
   
   # subsets the groups we want, then normalizes row-wise (per protein)
   hm_df_test <- reactive({
     h_i <- heatmap_index_test()
     df_sig <- base::subset(df_sig_test(), df_sig_test()$sig == TRUE)
     df_mat <- as.matrix(df_sig[,1:ncol(dfres())])
     df_mat_2 <- df_mat[,h_i]
     
     cal_z_score <- function(x) {
       (x - mean(x)) / sd(x)
     }
     
     hm_df <- t(apply(df_mat_2, 1, cal_z_score))
     return(hm_df)
   })
   
   # calculates distance matrix of selected groups
   datadist_test <- reactive({
     number_obs <- base::nrow(df_trans()[df_trans()$group %in% input$groups_heatmap, ])
     datacor <- cor(t(hm_df_test()[,1:number_obs]))
     datadist <- as.dist(1-datacor)
     return(datadist)
   })
   
   # output: hierarchical clustering object by distance matrix
   hclust_test <- reactive({
     my_hclust_gene2 <- hclust(datadist_test(), method="average")
     return(my_hclust_gene2)
   })
   
   # output: list of each protein and which cluster it falls under
   cluster_df_test <- reactive({
     cluster <- cutree(tree = hclust_test(), k=input$heatmap_clusters)
     cluster_df <- base::as.data.frame(cluster)
     cluster_df$cluster <- as.factor(cluster_df$cluster)
     return(cluster_df)
   })
   
   # returns a summary of the number of proteins per cluster
   heatmap_cluster_numbers <- eventReactive(input$heatmap_run, {
     amts <- c()
     for (i in 1:input$heatmap_clusters){
       amt <- length(cluster_df_test()[cluster_df_test()$cluster == i,])
       amts <- append(amts, base::paste("Proteins in cluster ", i, ": ", amt, sep=""))
     }
     return(HTML(base::paste(amts, collapse='<br/>')))
   })
   
   output$heatmap_cluster_numbers <- renderUI({
     return(heatmap_cluster_numbers())
   })
   
   # generates the heatmap
   heatmap_test <- eventReactive(input$heatmap_run, {
     num_groups <- length(input$groups_heatmap)
     number_obs <- base::nrow(df_trans()[df_trans()$group %in% input$groups_heatmap, ])
     title <- ("Heatmap")
     
     df_sample_c <- data.frame(group=rep(names(), times()))
     df_sample_c <- base::subset(df_sample_c, df_sample_c$group %in% input$groups_heatmap)
     row.names(df_sample_c) <- base::colnames(hm_df_test())
     p <- pheatmap(hm_df_test(), cluster_rows = hclust_test(), cluster_cols=FALSE, main=title, annotation_row = cluster_df_test(), annotation_col = df_sample_c, show_rownames = FALSE)
     return(p)
   })
   
   output$heatmap <- renderPlot({
     return(heatmap_test())
   })
   
   # description of heatmap results and parameters
   heatmap_info <- eventReactive(input$heatmap_run, {
     nprot <- nrow(cluster_df_test())
     str1 <- base::paste("This heatmap was generated using ", nprot, " proteins.", sep="")
     str2 <- base::paste("Significant proteins were determined as having log fold change greater than ", input$fc_cutoff, " and a t-test p-value of ", input$sig_cutoff, " after adjustment: ", input$adjustment, ".", sep="")
     str3 <- base::paste("The comparisons were done pairwise between groups ", base::paste(input$groups_heatmap, collapse=", "), ".", sep="")
     return(HTML(base::paste(str1, str2, str3, sep='<br/>')))
   })
   
   output$heatmap_info <- renderUI({
     return(heatmap_info())
   })
   
   # plot a graph of number of clusters vs silhouette width
   heatmap_silhouette <- eventReactive(input$heatmap_run, {
     silhouette_score <- function(k){
       my_hclust <- hclust(datadist_test(), method="average")
       cluster <- cutree(tree = my_hclust, k=k)
       ss <- silhouette(cluster, datadist_test())
       mean(ss[, 3])
     }
     k <- 2:20
     avg_sil <- base::sapply(k, silhouette_score)
     
     p <- plot(k, avg_sil)
     return(p)
   })
   
   output$optimal_hm <- renderPlot({
     return(heatmap_silhouette())
   })

   
  
   #################### PRINCIPAL COMPONENT ANALYSIS ####################
   # takes top/bottom loadings for PCA and uses them as radio buttons to plot distribution
   observe({
     top <- pca_loadings_top()[1:3,1]
     bot <- pca_loadings_bot()[1:3,1]
     combined <- c(top, bot, "Other")
     
     updateRadioButtons(session = session,
                        inputId = "select_loading",
                        choices = combined,
                        selected = combined[1])
   })
   
   # orders dataframe by coefficient of variation
   cv_sorted <- reactive({
     unq <- base::unique(metadf()$V2)
     
     subs <- base::subset(df_trans(), df_trans()$group == unq[1])
     subs <- apply(subs[1:(ncol(subs)-1)], 2, mean)
     
     for (i in 2:length(unq)) {
       s <- base::subset(df_trans(), df_trans()$group == unq[i])
       m <- apply(s[1:(ncol(s)-1)], 2, mean)
       subs <- base::rbind(subs, m)
     }
     
     subs <- t(subs)
     
     means <- apply(subs, 1, mean)
     stdev <- apply(subs, 1, sd)
     pcv <- stdev/means*100
     
     dfr <- dfres()
     dfr$pcv <- pcv
     
     sorted <- dfr[base::order(-dfr$pcv),]
     return(sorted)
   })
   
   # takes only the proteins that have the n largest cv's
   pca_matrix_t <- reactive({
     filtered <- cv_sorted()[(1:input$pca_topcv),]
     filtered <- base::subset(filtered, select = -c(pcv))
     
     df_t <- t(filtered)
   })
   
   pca_summary <- reactive({
     pca <- prcomp(pca_matrix_t(), center=TRUE, scale.=TRUE)
     pca_sum <- summary(pca)
     return(pca_sum)
   })

   pca <- eventReactive(input$pca_run, {
     pca <- prcomp(pca_matrix_t(), center=TRUE, scale.=TRUE)
     # takes the rotation value
     df_pca <- base::as.data.frame(pca$x)
     
     # assigns group to each replicate
     df_pca$sample <- rep(names(), times())
     df_pca$sample <- factor(df_pca$sample, levels=names())
     
     # assigns a color to each group
     length_times <- length(times())
     color <- c()
     for(i in 1:length(df_pca$sample)){
       index <- which(names() == df_pca$sample[i])
       color <- append(color, rainbow(length_times)[index])
     }
     df_pca$color <- color
     
     if(input$pc_x == 1 & input$pc_y == 2) {p <- ggplot(df_pca, aes(x=PC1, y=PC2, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times()))) 
     if(input$pc_x == 1 & input$pc_y == 3) {p <- ggplot(df_pca, aes(x=PC1, y=PC3, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 1 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC1, y=PC4, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 1 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC1, y=PC5, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 2 & input$pc_y == 3) {p <- ggplot(df_pca, aes(x=PC2, y=PC3, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 2 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC2, y=PC4, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 2 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC2, y=PC5, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 3 & input$pc_y == 4) {p <- ggplot(df_pca, aes(x=PC3, y=PC4, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 3 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC3, y=PC5, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     if(input$pc_x == 4 & input$pc_y == 5) {p <- ggplot(df_pca, aes(x=PC4, y=PC5, col=sample)) + geom_point(shape=19, size=7, alpha=0.7)} + scale_colour_manual(values = rainbow(length(times())))
     return(p)
   })
   
   # 2d plot of pca
   output$pca <- renderPlot({
     return(pca())
   })
   
   # 3d plot of pca
   output$pca3d <- renderRglwidget({
     obs <- observations()
     
     df_t <- pca_matrix_t()
     pca <- prcomp(df_t, center=TRUE, scale.=TRUE)       
     df_pca <- base::as.data.frame(pca$x)
     
     df_pca$sample <- rep(names(), times())
     
     length_times <- length(times())
     color <- c()
     for(i in 1:length(df_pca$sample)){
       index <- which(names() == df_pca$sample[i])
       color <- append(color, rainbow(length_times)[index])
     }
     df_pca$color <- color

     try(close3d())
     
     if(input$pc_x == 1 & input$pc_y == 2 & input$pc_z == 3) {plot3d(df_pca$PC1, df_pca$PC2, df_pca$PC3, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC2", zlab="PC3")}
     if(input$pc_x == 1 & input$pc_y == 2 & input$pc_z == 4) {plot3d(df_pca$PC1, df_pca$PC2, df_pca$PC4, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC2", zlab="PC4") }
     if(input$pc_x == 1 & input$pc_y == 2 & input$pc_z == 5) {plot3d(df_pca$PC1, df_pca$PC2, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC2", zlab="PC5")}
     if(input$pc_x == 1 & input$pc_y == 3 & input$pc_z == 4) {plot3d(df_pca$PC1, df_pca$PC3, df_pca$PC4, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC3", zlab="PC4")}
     if(input$pc_x == 1 & input$pc_y == 3 & input$pc_z == 5) {plot3d(df_pca$PC1, df_pca$PC3, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC3", zlab="PC5")}
     if(input$pc_x == 1 & input$pc_y == 4 & input$pc_z == 5) {plot3d(df_pca$PC1, df_pca$PC4, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC1", ylab="PC4", zlab="PC5")}
     if(input$pc_x == 2 & input$pc_y == 3 & input$pc_z == 4) {plot3d(df_pca$PC2, df_pca$PC3, df_pca$PC4, col=df_pca$color, type="s", size=2, xlab="PC2", ylab="PC3", zlab="PC4")}
     if(input$pc_x == 2 & input$pc_y == 3 & input$pc_z == 5) {plot3d(df_pca$PC2, df_pca$PC3, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC2", ylab="PC3", zlab="PC5")}
     if(input$pc_x == 2 & input$pc_y == 4 & input$pc_z == 5) {plot3d(df_pca$PC2, df_pca$PC4, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC2", ylab="PC4", zlab="PC5")}
     if(input$pc_x == 3 & input$pc_y == 4 & input$pc_z == 5) {plot3d(df_pca$PC3, df_pca$PC4, df_pca$PC5, col=df_pca$color, type="s", size=2, xlab="PC3", ylab="PC4", zlab="PC5")}
     
     M <- par3d("userMatrix")
     
     legend3d("topright", legend = unique(df_pca$sample), col=rainbow(length(times())), pch = 16, cex=0.9, inset=c(0.02))
     
     if (!rgl.useNULL())
       play3d(par3dinterp(time = (0:2)*0.75, userMatrix = list(M,
                                                               rotate3d(M, pi/2, 1, 0, 0),
                                                               rotate3d(M, pi/2, 0, 1, 0) ) ), 
              duration = 10)
     
     rglwidget()
   })
   
   observeEvent(input$pca3d_snapshot, {
     rgl.snapshot("www/pca3d.png", fmt = "png", top = TRUE)
   })
   
   # scree plot
   pca_importance <- eventReactive(input$pca_run, {
     pca_summary <- pca_summary()
     imp <- pca_summary[["importance"]][2:3,]
     imp <- t(imp)
     imp <- cbind.data.frame(PC = base::rownames(imp), imp)
     base::colnames(imp) <- c("PC", "Proportion", "Cumulative")
     imp$PC <- factor(imp$PC, levels=imp$PC)
     p <- ggplot(imp) + geom_col(aes(x=PC, y=Proportion)) + geom_line(aes(x=PC, y=Cumulative), color='red', size=1, group=1) +
       geom_point(aes(x=PC, y=Cumulative), color='red', size=3)
     return(p)
   })
   
   output$pca_importance <- renderPlot({
     return(pca_importance())
   })
   
   # shows a table of top/bottom loadings
   pca_loadings_top <- reactive({
     if(input$pc_load == "PC1") {n <- 1}
     else if(input$pc_load == "PC2") {n <- 2}
     else if(input$pc_load == "PC3") {n <- 3}
     else if(input$pc_load == "PC4") {n <- 4}
     else if(input$pc_load == "PC5") {n <- 5}
     
     pca_summary <- pca_summary()
     
     top_index <- base::sort(pca_summary[["rotation"]][,n], index.return=TRUE, decreasing=TRUE)$ix[1:10]
     top_values <- base::sort(pca_summary[["rotation"]][,n], index.return=TRUE, decreasing=TRUE)$x[1:10]
     
     top_genes <- df()[top_index,1]
     top_table <- base::as.data.frame(top_values)
     top_table$proteins <- top_genes
     top_table <- top_table[, c(2,1)]
     return(top_table)
   })
   
   pca_loadings_bot <- reactive({
     if(input$pc_load == "PC1") {n <- 1}
     else if(input$pc_load == "PC2") {n <- 2}
     else if(input$pc_load == "PC3") {n <- 3}
     else if(input$pc_load == "PC4") {n <- 4}
     else if(input$pc_load == "PC5") {n <- 5}
     
     pca_summary <- pca_summary()
     
     bot_index <- base::sort(pca_summary[["rotation"]][,n], index.return=TRUE, decreasing=FALSE)$ix[1:10]
     bottom_values <- base::sort(pca_summary[["rotation"]][,n], index.return=TRUE, decreasing=FALSE)$x[1:10]
     
     bot_genes <- df()[bot_index,1]
     bot_table <- base::as.data.frame(bottom_values)
     bot_table$proteins <- bot_genes
     bot_table <- bot_table[, c(2,1)]
     return(bot_table)
   })
   
   output$pca_loadings_bot <- renderTable({
     return(pca_loadings_bot())
   }, digits=5)
   
   output$pca_loadings_top <- renderTable({
     return(pca_loadings_top())
   }, digits=5)
   
   # for a particular PC, we can examine how the proteins are loaded
   pca_loadings_plot <- reactive({
     if(input$pc_load == "PC1") {n <- 1}
     else if(input$pc_load == "PC2") {n <- 2}
     else if(input$pc_load == "PC3") {n <- 3}
     else if(input$pc_load == "PC4") {n <- 4}
     else if(input$pc_load == "PC5") {n <- 5}
     
     pca_summary <- pca_summary()
     
     sorted <- base::sort(pca_summary[["rotation"]][,n], decreasing=TRUE)
     
     p <- plot(sorted, xaxt = "n", yaxp = round(c(min(sorted), max(sorted), 10), digits=1), las=2, ylab="", xlab="", main = input$pc_load)
     return(p)
   })
   
   output$pca_loadings_plot <- renderPlot({
     pca_loadings_plot()
   })
   
   # plot a particular protein's distribution by group
   loadings_scatter <- eventReactive(input$plot_loadings_scatter, {
     if (input$select_loading == "Other") {
       table <- base::subset(df(), df()[,1] == input$select_loading_text)
     }
     else {
       table <- base::subset(df(), df()[,1] == input$select_loading)
     }
     
     expression <- as.numeric(table[,-1])
     
     if (input$normalize %in% df()[,1]) {
       standard <- base::subset(df(), df()[,1] == input$normalize)
       standard <- as.numeric(standard[,-1])
       expression <- expression/standard
     }
     
     group <- rep(names(), times())
     table_t <- cbind.data.frame(group, expression)
     
     p <- ggplot(table_t, aes(x=group, y=expression)) + geom_point()
     if (input$select_loading != "Other") {title <- input$select_loading}
     else if (input$select_loading == "Other") {title <- input$select_loading_text}
     
     if (input$normalize %in% df()[,1]) {
       title <- base::paste(title, "normalized to", input$normalize, sep=" ")
     }
     
     p <- p + ggtitle(title)
     return(p)
   })

   output$loadings_scatter <- renderPlot({
     return(loadings_scatter())
   })
   
   ######################## GENE ONTOLOGY ##########################
   # if choosing to use heatmap as list, then specify which of the heatmap's clusters are to be used.
   observe({
     choice <- seq(1:input$heatmap_clusters)
     updateCheckboxGroupInput(session = session,
                              inputId = "go_heatmap_clusters_choices",
                              choices = choice)
   })
   
   # processing if list=top cv
   go_listcv <- reactive({
     listcv <- df()[base::rownames(cv_sorted()),1]
     l <- listcv[1:input$go_topcv]
     return(as.vector(l))
   })
   
   # processing if list=heatmap
   go_listheatmap <- reactive({
     df <- cbind.data.frame(volcano_list_test(), cluster_df_test())
     df <- subset(df, df$cluster %in% input$go_heatmap_clusters_choices)
     return(df)
   })
   
   # processing if list=uploaded file
   go_listupload <- reactive({
     df <- req(input$go_file)
     df <- read.csv(input$go_file$datapath,
                    header = TRUE,
                    sep = ",",
                    quote = '"',
                    encoding = "UTF-8")
     return(df[,1])
   })
   
   # input: list of proteins, output: enrichment object - contains pathway name, pval, gene ratio, and associated genes
   go_enrichment <- eventReactive(input$go_calculate, {
     enrichment<- function(x, y){
       plot=enrichGO(x, org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont=y,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.02,
                     minGSSize=10,
                     maxGSSize = 500,
                     readable = FALSE,
                     pool = FALSE)
       p <- plot
       print("enrich done")
       return(p)
     }
     
     if (input$go_proteinlist == "Heatmap"){
       list1 <- go_listheatmap()$protein
     }
     else if (input$go_proteinlist == "Top CV"){
       list1 <- go_listcv()
     }
     else if (input$go_proteinlist == "Upload File"){
       list1 <- go_listupload()
     }
     
     list2 <- c()
     for (i in 1:length(list1)) {
       list1[i] <- unlist(strsplit(list1[i], split="_", fixed=TRUE))[1]
     }
    
     return(enrichment(list1, input$GO_ont))
   })
   
   # enriched pathway names
   go_descriptions <- eventReactive(input$go_calculate, {
     return(go_enrichment()$Description)
   })
   
   output$go_descriptions <- renderText({
     return(go_descriptions())
   })
   
   # selection for heatmap viewing
   observeEvent(input$go_calculate, {
     choice <- go_descriptions()
     updateSelectInput(session = session,
                       inputId = "go_genelist",
                       choices = choice)
   })
   
   # for each enriched pathway, return a list of genes associated with the pathway
   go_geneID <- eventReactive(input$go_calculate, {
     p <- go_enrichment()$geneID
     p2 <- list()
     for(i in 1:length(p)){
       test <- strsplit(p[i], split="/", fixed=TRUE)
       p2 <- list.append(p2, test)
     }
     names(p2) <- go_descriptions()
     return(p2)
   })
   
   # for each enriched pathway, return the p-value
   go_pvals <- eventReactive(input$go_calculate, {
     pvals <- go_enrichment()$p.adjust
     return(pvals)
   })
   
   # after choosing ONE pathway to view, this will return the genes in the pathway
   go_sig_genes <- eventReactive(input$go_heatmap_run, {
     return(go_geneID()[[input$go_genelist]])
   })
   
   # generates heatmap of proteins in the selected pathway. Can choose which groups to display
   go_new_heatmap <- eventReactive(input$go_heatmap_run, {
     oldnames <- df()[,1]
     newnames <- c()
     
     # if gene names have "_HUMAN" in them, then take only the name
     for (i in 1:length(oldnames)) {
        newnames <- list.append(newnames, unlist(strsplit(oldnames[i], split="_", fixed=TRUE))[1])
     }
     d <- cbind(newnames, df()[,2:ncol(df())])
     d <- subset(d, d[,1] %in% go_sig_genes()[[1]])
     
     # remove duplicated rows of proteins
     rn <- d[,1]
     nodup <- !duplicated(rn)
     rn <- rn[nodup]
     d <- d[,2:ncol(d)]
     d <- d[nodup,]
     rownames(d) <- rn
     
     # normalize protein per row
     cal_z_score <- function(x) {
       (x - mean(x)) / sd(x)
     }
     
     d <- t(apply(d, 1, cal_z_score))
     
     # subset dataframe by selected groups to view
     dft <- df_trans()
     dft <- subset(dft, dft$group %in% input$go_heatmap_groups)
     obs_to_keep <- rownames(dft)

     d <- d[,obs_to_keep]

     df_sample_c <- data.frame(group=rep(names(), times()))
     df_sample_c <- subset(df_sample_c, df_sample_c$group %in% input$go_heatmap_groups)

     rownames(df_sample_c) <- colnames(d)
     
     pval_index <- which(go_descriptions() == input$go_genelist)
     title <- base::paste(input$go_genelist, ", adjusted p-value: ", go_pvals()[pval_index], sep="")
     
     return(pheatmap(d, main=title, labels_col=colnames(df()[,-1]), cluster_cols=FALSE, annotation_col = df_sample_c))
   })
   
   output$go_new_heatmap <- renderPlot({
     return(go_new_heatmap())
   })
   
   # shows dotplot: x-axis/size=gene ratio, color=p-value
   go_dotplot <- eventReactive(input$go_calculate, {
     return(enrichplot::dotplot(go_enrichment()))
     })
   
   output$go_dotplot <- renderPlot({
     go_dotplot()
     })
   
   # description of dotplot
   go_info <- eventReactive(input$go_calculate, {
     str1 <- "This list used a total of"
     if(input$go_proteinlist == "Heatmap"){
       str2 <- base::paste(base::nrow(go_listheatmap()))
       str4 <- base::paste(input$go_heatmap_clusters_choices, collapse=", ")
       str4 <- base::paste("clusters:", str4)
     }
     if(input$go_proteinlist == "Top CV"){
       str2 <- base::paste(length(go_listcv()))
       str4 <- base::paste(" the top ", input$go_topcv, " largest CV.", sep="")
     }
     if(input$go_proteinlist == "Upload File"){
       str2 <- length(go_listupload())
       str4 <- "the uploaded file."
     }
     str3 <- "proteins from"
     str5 <- base::paste("The gene ontology used was ", input$GO_ont, ".", sep="")
     sub1 <- HTML(base::paste(str1, str2, str3, str4,  sep = ' '))
     return(HTML(base::paste(sub1, str5, sep = '<br/>')))
   })
   
   output$go_info <- renderUI({
     return(go_info())
   })
   
   ############################### GSEA ##############################
   # if selected list is from heatmap, then select clusters to use
   observe({
     choice <- seq(1:input$heatmap_clusters)
     updateCheckboxGroupInput(session = session,
                              inputId = "gsea_heatmap_clusters_choices",
                              choices = choice,
                              selected = choice)
   })
   
   # processing if list=heatmap
   gsea_listheatmap <- reactive({
     df <- cbind.data.frame(volcano_list_test(), cluster_df_test())
     df <- subset(df, df$cluster %in% input$gsea_heatmap_clusters_choices)
     df <- df[order(df[,2], decreasing=TRUE),]
     return(df)
   })
   
   # processing if list=upload file
   gsea_listupload <- reactive({ 
     df <- req(input$gsea_file)
     df <- read.csv(input$gsea_file$datapath,
                  header = TRUE,
                  sep = ",",
                  quote = '"',
                  encoding = "UTF-8")
     df <- df[order(df[,2], decreasing=TRUE),]
     return(df)
   })
   
   # returns a gsea object - includes gene sets and genes in each set
   gsea <- eventReactive(input$gsea_calculate, {
     geneset <- msigdbr(species= "Homo sapiens", category=input$gsea_set)
     geneset <- dplyr::select(geneset, gs_name, gene_symbol)
     
     if(input$gsea_proteinlist == "Heatmap") {
        protlist <- gsea_listheatmap()[,2]
        protnames <- gsea_listheatmap()[,1]
     }
     else if(input$gsea_proteinlist == "Upload File"){
       protlist <- gsea_listupload()[,2]
       protnames <- gsea_listupload()[,1]
     }
     # this uses all proteins in the dataframe
     else if(input$gsea_proteinlist == "All") {
       protlist <- df_tt_test()[,(ncol(df_tt_test())-2)]
       protnames <- df()[,1]
     }
     
     for (i in 1:length(protnames)){
       protnames[i] <- unlist(strsplit(protnames[i], split="_", fixed=TRUE))[1]
     }
     names(protlist) <- protnames
     protlist[!is.finite(protlist)] <- NA
     protlist <- na.omit(protlist)
     protlist <- sort(protlist, decreasing=TRUE)
     
     test_gse <- GSEA(protlist, exponent =1, minGSSize = 10, maxGSSize = 2000, eps=1e-10,
                    pvalueCutoff =0.05, pAdjustMethod= "BH", TERM2GENE = geneset)
     return(test_gse)
   })
   
   gsea_results <- eventReactive(input$gsea_calculate, {
     return(gsea()@result)
   })
   
   # returns a description of each gene set
   gsea_descriptions <- eventReactive(input$gsea_calculate, {
     g <- gsea_results()$ID
     return(g)
   })
   
   # selection for heatmap/running enrichment score viewing
   observeEvent(input$gsea_calculate, {
     choice <- gsea_descriptions()
     updateSelectInput(session = session,
                       inputId = "gsea_genelist",
                       choices = choice)
   })
   
   # for each gene set, return the list of genes associated
   gsea_geneID <- eventReactive(input$gsea_calculate, {
     p <- gsea_results()$core_enrichment
     p2 <- list()
     for(i in 1:length(p)){
       test <- strsplit(p[i], split="/", fixed=TRUE)
       p2 <- list.append(p2, test)
     }
     names(p2) <- gsea_descriptions()
     return(p2)
   })
   
   # when ONE gene set is chosen, return that list of genes
   gsea_sig_genes <- eventReactive(input$gsea_heatmap_run, {
     return(gsea_geneID()[[input$gsea_genelist]])
   })
   
   # running enrichment score plot, input: gsea object
   gsea_running <- eventReactive(input$gsea_heatmap_run, {
     descriptions <- as.vector(gsea_descriptions())
     p <- gseaplot(gsea(), geneSetID = which(descriptions %in% input$gsea_genelist))
     return(p)
   })
   
   output$gsea_running <- renderPlot({
     return(gsea_running())
   })
   
   # heatmap of selected groups showing proteins that are in the gene set
   gsea_new_heatmap <- eventReactive(input$gsea_heatmap_run, {
     oldnames <- df()[,1]
     newnames <- c()
     for (i in 1:length(oldnames)) {
       newnames <- list.append(newnames, unlist(strsplit(oldnames[i], split="_", fixed=TRUE))[1])
     }
     
     d <- cbind(newnames, df()[,2:ncol(df())])
     d <- subset(d, d[,1] %in% gsea_sig_genes()[[1]])
     
     rn <- d[,1]
     nodup <- !duplicated(rn)
     rn <- rn[nodup]
     d <- d[,2:ncol(d)]
     d <- d[nodup,]
     rownames(d) <- rn
     
     cal_z_score <- function(x) {
       (x - mean(x)) / sd(x)
     }
     
     d <- t(apply(d, 1, cal_z_score))
     
     dft <- df_trans()
     dft <- subset(dft, dft$group %in% input$gsea_heatmap_groups)
     obs_to_keep <- rownames(dft)
     
     d <- d[,obs_to_keep]
     
     df_sample_c <- data.frame(group=rep(names(), times()))
     df_sample_c <- subset(df_sample_c, df_sample_c$group %in% input$gsea_heatmap_groups)
     rownames(df_sample_c) <- colnames(d)
     
     return(pheatmap(d, main=input$gsea_genelist, clustering_distance_rows = "correlation", cluster_cols=FALSE, annotation_col = df_sample_c))
   })
   
   output$gsea_new_heatmap <- renderPlot({
     return(gsea_new_heatmap())
   })
   
   ############################ POPUP ##############################
   dataModal <- function(failed = FALSE) {
     modalDialog(
       checkboxGroupInput("include", "Include in Report: ",
                    choices = c(Venn = "venn",
                                Violin = "violin",
                                Correlation = "correlation",
                                Volcano = "volcano",
                                Heatmap = "heatmap",
                                PCA = "pca"),
                    selected = c("venn", "violin", "correlation", "volcano", "heatmap", "pca")
                    ),
       textInput("download_title", "Title of document: "),

       footer = tagList(
         modalButton("Cancel"),
         downloadButton("downloadReport", "Confirm")
       )
     )
   }
   
   observeEvent(input$downloadoptions, {
     showModal(dataModal())
   })
   
   ########################### DOWNLOAD FIGURES #######################
   output$venndiagram_download <- downloadHandler(
     filename = function() {
       base::paste('venn-diagram-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(venn_diagram())
       dev.off()
     }
   )
   
   output$violin_download <- downloadHandler(
     filename = function() {
       base::paste('violin-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(violin())
       dev.off()
     }
   )
   
   output$violin_agg_download <- downloadHandler(
     filename = function() {
       base::paste('violin-aggregate-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(violin_agg())
       dev.off()
     }
   )
   
   output$correlation_download <- downloadHandler(
     filename = function() {
       base::paste('correlation-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(corrmap())
       dev.off()
     }
   )
   
   output$scatter_download <- downloadHandler(
     filename = function() {
       base::paste('scatterplot-matrix-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(scattermat())
       dev.off()
     }
   )
   
   output$volcano_download <- downloadHandler(
     filename = function() {
       base::paste('volcano-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(volcano())
       dev.off()
     }
   )
   
   output$heatmap_download <- downloadHandler(
     filename = function() {
       base::paste('heatmap-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(heatmap())
       dev.off()
     }
   )
   
   output$pca_download <- downloadHandler(
     filename = function() {
       base::paste('pca-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(pca())
       dev.off()
     }
   )
   
   output$pca_loadings_plot_download <- downloadHandler(
     filename = function() {
       base::paste('pca-loadings-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(pca_loadings_plot())
       dev.off()
     }
   )
   
   output$pca_importance_download <- downloadHandler(
     filename = function() {
       base::paste('scree-plot-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(pca_importance())
       dev.off()
     }
   )
   
   output$go_dotplot_download <- downloadHandler(filename = function() {
     base::paste('GO-dotplot-', Sys.Date(), '.png', sep='')
   },
   content = function(file) {
     png(file)
     print(go_dotplot())
     dev.off()
   })
   
   output$go_heatmap_download <- downloadHandler(
     filename = function() {
       base::paste('GO-', base::paste(input$go_genelist, collapse='-'), '-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(go_new_heatmap())
       dev.off()
     }
   )
   
   output$gsea_heatmap_download <- downloadHandler(
     filename = function() {
       base::paste('GSEA-', base::paste(input$gsea_genelist, collapse='-'), '-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(gsea_new_heatmap())
       dev.off()
     }
   )
   
   output$gsea_running_download <- downloadHandler(
     filename = function() {
       base::paste('GSEA-running-enrichment-score-', Sys.Date(), '.png', sep='')
     },
     content = function(file) {
       png(file)
       print(gsea_running())
       dev.off()
     }
   )
   
   ## BUTTONS
   get_group1 <- reactive({return(base::paste(input$groups1, collapse='/'))})
   get_group2 <- reactive({return(base::paste(input$groups2, collapse='/'))})
   get_sig <- reactive({return(input$sig_cutoff)})
   get_fc <- reactive({return(input$fc_cutoff)})
   get_adj <- reactive({return(input$adjustment)})
   get_topcv <- reactive({return(input$pca_topcv)})
   get_downloadtitle <- reactive({return(input$download_title)})
   get_loadings <- reactive({return(pca_loadings_plot())})
   
   output$volcano_button <- downloadHandler(
     filename = function() {
       base::paste('significant-proteins-', Sys.Date(), '.csv', sep='')
     },
     content = function(file) {
       g1 <- base::paste(input$groups1, collapse='/')
       g2 <- base::paste(input$groups2, collapse='/')
       
       meta <- base::paste("#",Sys.time(),", ",g1," vs ",g2,", alpha: ",input$sig_cutoff,", adjustment: ",input$adjustment,", fold change: ", input$fc_cutoff, sep='')
       if(input$volcanodownloadformat == "details") {
          write.table(meta, file, append=TRUE, row.names=FALSE, col.names=FALSE)
       }
       write.table(volcano_list(), file, append=TRUE, row.names=FALSE, sep=',')
      
       if(input$volcanodownloadformat == "list") {
         write.table(volcano_list()[,1], file, row.names=FALSE, col.names=FALSE, append=FALSE)
       }
     }
   )
   
   output$heatmap_button <- downloadHandler(
     filename = function() {
       base::paste('protein-clusters-', Sys.Date(), '.csv', sep='')
     },
     content = function(file) {
       cluster <- cluster_df_test()
       protein <- df()[rownames(cluster),1]
       t <- cbind.data.frame(cluster, protein)
       t <- t[order(t$cluster),]
       nclusters <- length(base::unique(t$cluster))
       
       meta <- base::paste("#",Sys.time(), ", alpha: ",input$sig_cutoff,", adjustment: ",input$adjustment,", fold change: ", input$fc_cutoff, sep='')
       clusterinfo <- base::paste("# Hierarchical clustering of proteins into ", nclusters, " clusters", sep="")
       
       write.table(meta, file, append=TRUE, row.names=FALSE, col.names=FALSE)
       write.table(clusterinfo, file, append=TRUE, row.names=FALSE, col.names=FALSE)
       write.table(t, file, append=TRUE, row.names=FALSE, sep=',')
     }
   )
   
   output$loadings_button <- downloadHandler(
     filename = function() {
       base::paste('pca-loadings-', Sys.Date(), '.csv', sep='')
     },
     content = function(file) {
       meta <- base::paste("#Loadings for ", input$pc_load, ", Top ", input$pca_topcv, " Proteins",  sep='')
       meta2 <- base::paste("#Loadings for ", input$pc_load, ", Bottom ", input$pca_topcv, " Proteins",  sep='')
       write.table(meta, file, append=TRUE, row.names=FALSE, col.names=FALSE)
       write.table(pca_loadings_top(),file, append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
       write.table(meta2, file, append=TRUE, row.names=FALSE, col.names=FALSE)
       write.table(pca_loadings_bot(),file, append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
     }
   )
   
   output$downloadReport <- downloadHandler(
     filename = function() {
         base::paste('report-', Sys.Date(), '.html', sep='')
     },
     content = function(file) {
         # Copy the report file to a temporary directory before processing it, in
         # case we don't have write permissions to the current working dir (which
         # can happen when deployed).
         tempReport <- file.path(tempdir(), "report.Rmd")
         file.copy("report.Rmd", tempReport, overwrite = TRUE)
         
         # Set up parameters to pass to Rmd document
         params <- list(metadf_summary = metadf_summary,
                        contents_table = contents_table,
                        doc_title = get_downloadtitle)

         if("correlation" %in% input$include) {
           params <- list.append(params, contents_table = contents_table,
                                    correlation = corrmap,
                                    scattermat = scattermat,
                                    show_correlation = TRUE)}
         if("venn" %in% input$include) {
           params <- list.append(params, venn_diagram = venn_diagram, 
                                 show_venn = TRUE)}
         if("violin" %in% input$include) {
           params <- list.append(params, violin = violin, 
                                    violin_agg = violin_agg,
                                    show_violin = TRUE)}
         if("volcano" %in% input$include) {
           params <- list.append(params, volcano = volcano,
                            group1 = get_group1,
                            group2 = get_group2,
                            significance = get_sig,
                            fc = get_fc,
                            adjustment = get_adj,
                            show_volcano = TRUE)}
         if ("heatmap" %in% input$include) {
           params <- list.append(params, group1 = get_group1,
                            group2 = get_group2,
                            significance = get_sig,
                            fc = get_fc,
                            adjustment = get_adj,
                            heatmap = heatmap_test,
                            show_heatmap = TRUE)}
         if("pca" %in% input$include) {
           params <- list.append(params, pca = pca,
                            pca_loadings_plot = pca_loadings_plot,
                            pca_importance = pca_importance,
                            topcv = get_topcv,
                            show_pca = TRUE)}
         
         rmarkdown::render(tempReport, output_file = file,
                           params = params,
                           envir = new.env(parent = globalenv())
         )
       }
     )
   
   # close the R session when Chrome closes
   shinyServer(function(input, output, session){
     session$onSessionEnded(function() {
       stopApp()
     })
   })
}

# Create Shiny app ----
shinyApp(ui, server)

# DEPLOY: 
# options(repos = BiocManager::repositories())
# getOption("repos")
# rsconnect::deployApp('C:/Users/patri/OneDrive/Documents/Columbia/RaiLab/App-1')
