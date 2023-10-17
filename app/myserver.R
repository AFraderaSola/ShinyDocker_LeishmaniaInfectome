############################
######## Libraries #########
############################

library(shiny)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(ggrepel)
library(kohonen)
library(effectsize)
library(tibble)
library(tidyr)
library(rMQanalysis)
library(shinyBS)
library(DT)

############################
####### Variables ##########
############################

hours <- c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h")

df_infantum <- read.csv("./InputFiles/Database_data_infantum.csv")

inf_sub <- unique(df_infantum$TimePoint)

for (i in 1:length(hours)) {
  
  df_infantum$TimePoint <- gsub(pattern = inf_sub[i], replacement = hours[i], df_infantum$TimePoint)
  
}

df_infantum$Replica <- as.character(df_infantum$Replica)

df_major <- read.csv("./InputFiles/Database_data_major.csv")

mjr_sub <- unique(df_major$TimePoint)

for (i in 1:length(hours)) {
  
  df_major$TimePoint <- gsub(pattern = mjr_sub[i], replacement = hours[i], df_major$TimePoint)
  
}

df_major$Replica <- as.character(df_major$Replica)

df_mexicana <- read.csv("./InputFiles/Database_data_mexicana.csv")

mex_sub <- unique(df_mexicana$TimePoint)

for (i in 1:length(hours)) {
  
  df_mexicana$TimePoint <- gsub(pattern = mex_sub[i], replacement = hours[i], df_mexicana$TimePoint)
  
}

df_mexicana$Replica <- as.character(df_mexicana$Replica)

test_infantum <- read.csv("./InputFiles/Test_data_infantum.csv")

for (i in 1:length(hours)) {
  
  test_infantum$group1 <- gsub(pattern = inf_sub[i], replacement = hours[i], test_infantum$group1)
  test_infantum$group2 <- gsub(pattern = inf_sub[i], replacement = hours[i], test_infantum$group2)
  
}

test_major <- read.csv("./InputFiles/Test_data_major.csv")

for (i in 1:length(hours)) {
  
  test_major$group1 <- gsub(pattern = mjr_sub[i], replacement = hours[i], test_major$group1)
  test_major$group2 <- gsub(pattern = mjr_sub[i], replacement = hours[i], test_major$group2)
  
}

test_mexicana <- read.csv("./InputFiles/Test_data_mexicana.csv")

for (i in 1:length(hours)) {
  
  test_mexicana$group1 <- gsub(pattern = mex_sub[i], replacement = hours[i], test_mexicana$group1)
  test_mexicana$group2 <- gsub(pattern = mex_sub[i], replacement = hours[i], test_mexicana$group2)
  
}

full_infantum <- read.csv("./InputFiles/QuantifiedProteins_infantum.csv")

full_major <- read.csv("./InputFiles/QuantifiedProteins_major.csv")

full_mexicana <- read.csv("./InputFiles/QuantifiedProteins_mexicana.csv")

colors_infantum <- colorRampPalette(c("#000000", "#00A651"))(8)

colors_major <- colorRampPalette(c("#000000", "#EE2A7B"))(8)

colors_mexicana <- colorRampPalette(c("#000000", "#2B3990"))(8)

# hm_pal <- c(rep("white",each=5), brewer.pal(9,"Blues"))

hm_pal <- colorRampPalette(brewer.pal(9,"Blues"))(15)

orth_df <- read.csv("./InputFiles/Leishmania_Orthologs.csv")

colnames(orth_df)[2] <- "#Species"

colnames(orth_df)[4] <- "Alg. Conn."

############################
######### Script ###########
############################

server <- function(input, output, session) {
  
  #### Data functions ####
  
  # Modify the different data sets in function of
  # user's input
  
  
  
  ##### L. infantum ####
  
  ## Main data table
  
  infantum_MainDataTable1 <- reactive({
    my_main_table_data <- df_infantum 
    if(input$infantum_ui_datasetID != 'All') {
      my_main_table_data <- df_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    my_main_table_data
  })
  
  infantum_MainDataTable2 <- reactive({
    my_main_table_data <- infantum_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs) 
    my_main_table_data
  })
  
  ## Test data table
  
  infantum_TestDataTable1 <- reactive({
    my_test_table_data <- test_infantum 
    if(input$infantum_ui_datasetID != 'All') {
      my_test_table_data <- test_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    my_test_table_data
  })
  
  infantum_TestDataTable2 <- reactive({
    my_test_table_data <- infantum_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  infantum_FullDataTable1 <- reactive({
    my_full_table_data <- full_infantum 
    if(input$infantum_ui_datasetID != 'All') {
      my_full_table_data <- full_infantum %>% filter(Species == input$infantum_ui_datasetID)
    }
    my_full_table_data
  })

  infantum_FullDataTable2 <- reactive({
    full_infantum2 <- infantum_FullDataTable1()
    if (length(input$infantum_ui_selectedIDs) != 0) {
      full_infantum2 <- full_infantum2 %>% filter(Majority.protein.IDs %in% input$infantum_ui_selectedIDs)
      full_infantum2
    }else{
      full_infantum2 <- data.frame()
    }
  })
  
  ##### L. major ####
  
  ## Main data table
  
  major_MainDataTable1 <- reactive({
    my_main_table_data <- df_major 
    if(input$major_ui_datasetID != 'All') {
      my_main_table_data <- df_major %>% filter(Species == input$major_ui_datasetID)
    }
    my_main_table_data
  })
  
  major_MainDataTable2 <- reactive({
    my_main_table_data <- major_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs) 
    my_main_table_data
  })
  
  ## Test data table
  
  major_TestDataTable1 <- reactive({
    my_test_table_data <- test_major 
    if(input$major_ui_datasetID != 'All') {
      my_test_table_data <- test_major %>% filter(Species == input$major_ui_datasetID)
    }
    my_test_table_data
  })
  
  major_TestDataTable2 <- reactive({
    my_test_table_data <- major_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  major_FullDataTable1 <- reactive({
    my_full_table_data <- full_major 
    if(input$major_ui_datasetID != 'All') {
      my_full_table_data <- full_major %>% filter(Species == input$major_ui_datasetID)
    }
    my_full_table_data
  })
  
  major_FullDataTable2 <- reactive({
    full_major2 <- major_FullDataTable1()
    if (length(input$major_ui_selectedIDs) != 0) {
      full_major2 <- full_major2 %>% filter(Majority.protein.IDs %in% input$major_ui_selectedIDs)
      full_major2
    }else{
      full_major2 <- data.frame()
    }
  })
  
  ##### L. mexicana ####
  
  ## Main data table
  
  mexicana_MainDataTable1 <- reactive({
    my_main_table_data <- df_mexicana 
    if(input$mexicana_ui_datasetID != 'All') {
      my_main_table_data <- df_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    my_main_table_data
  })
  
  mexicana_MainDataTable2 <- reactive({
    my_main_table_data <- mexicana_MainDataTable1()
    my_main_table_data <- my_main_table_data %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs) 
    my_main_table_data
  })
  
  ## Test data table
  
  mexicana_TestDataTable1 <- reactive({
    my_test_table_data <- test_mexicana 
    if(input$mexicana_ui_datasetID != 'All') {
      my_test_table_data <- test_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    my_test_table_data
  })
  
  mexicana_TestDataTable2 <- reactive({
    my_test_table_data <- mexicana_TestDataTable1()
    my_test_table_data <- my_test_table_data %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs) 
    my_test_table_data
  })
  
  ## Exploratory data
  
  mexicana_FullDataTable1 <- reactive({
    my_full_table_data <- full_mexicana 
    if(input$mexicana_ui_datasetID != 'All') {
      my_full_table_data <- full_mexicana %>% filter(Species == input$mexicana_ui_datasetID)
    }
    my_full_table_data
  })
  
  mexicana_FullDataTable2 <- reactive({
    full_mexicana2 <- mexicana_FullDataTable1()
    if (length(input$mexicana_ui_selectedIDs) != 0) {
      full_mexicana2 <- full_mexicana2 %>% filter(Majority.protein.IDs %in% input$mexicana_ui_selectedIDs)
      full_mexicana2
    }else{
      full_mexicana2 <- data.frame()
    }
  })
  
  #### Output functions ####
  
  #### L. infantum #### 
  
  #### Main table output  
  
  output$infantum_MainDataTable <- renderDataTable({
    infantum_MainDataTable1() %>%
      datatable(options = list(pageLength = 4, autoWidth = TRUE), 
                rownames = F) %>%
      formatRound('Value', digits = 2)
  })
  
  infantum_MainDataTable_proxy <- dataTableProxy('infantum_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$infantum_PCA <- renderPlot({
    
    # Load and format the data
    
    plot_df_infantum <- infantum_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_infantum))
    pca <- prcomp(t(na.omit(plot_df_infantum[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2 * -1
    scores$TimePoint <- factor(rep(c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"), each = 4), levels = c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"))
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    # Input plot
    
    infantum_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_point(size=3.5)+
      stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = TimePoint))+
      geom_text_repel(mapping = aes(label = Label),
                      size = 3,
                      fontface = 'bold',
                      color = 'black')+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_color_manual(values = colors_infantum)+
      scale_fill_manual(values = colors_infantum)+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))
    
    print(infantum_PCA)
    
  })
  
  output$infantum_HeatMap <- renderPlot({
    
    # Load and format the data
    
    plot_df_infantum <- infantum_FullDataTable1()
    lfq <- plot_df_infantum[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_infantum))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(hours)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = hours[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    # Input plot
    
    infantum_HeatMap <- gplots::heatmap.2(corr,
                                         trace="none",
                                         Colv=T,
                                         Rowv=T,
                                         dendrogram="column",
                                         density.info="density",
                                         srtCol=45,
                                         margins=c(10, 10), 
                                         key.xlab="Pearson's correlation coefficient",
                                         key.title="",
                                         keysize=1.5,
                                         col=hm_pal)
    
    print(infantum_HeatMap)
    
  })
  
  #### Clustering analysis output 
  
  output$infantum_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "inf"
    
    pg_LFQ <- infantum_FullDataTable1()
    
    hl_df <- infantum_FullDataTable2()
    
    # User input variables
    
    xdim <- input$infantum_xDIM
    
    ydim <- input$infantum_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#00A651"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    infantum_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
              
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      infantum_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(infantum_ClusterLinePlot)
    
  })
  
  output$infantum_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "inf"
    
    pg_LFQ <- infantum_FullDataTable1()
    
    hl_df <- infantum_FullDataTable2()
    
    # User input variables
    
    xdim <- input$infantum_xDIM
    
    ydim <- input$infantum_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_inf <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_inf$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#00A651"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_inf$data[[1]][which(som_model_inf$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LINF.*", data$ID),
               sub(".*", "L. infantum", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      data <- data[data$TimePoint == "0h",]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    infantum_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. infantum" = "#00A651",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    infantum_ClusterBarPlot <- ggplotly(infantum_ClusterBarPlot)
    
    print(infantum_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$infantum_Boxplot <- renderPlot({
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_infantum <- infantum_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_infantum <- infantum_MainDataTable2() 
    
    # Dummy plot
    
    infantum_Boxplot <-ggboxplot(data = plot_df_infantum, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_infantum)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- infantum_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_infantum) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_infantum <- jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        infantum_Boxplot <- ggboxplot(data = plot_df_infantum, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_infantum, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_infantum)+
          scale_color_brewer(palette="Dark2")+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        infantum_Boxplot <- list(infantum_Boxplot)
        
        boxplot_list <- c(boxplot_list, infantum_Boxplot)
        
      }
      
      infantum_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(infantum_Boxplot)
  })
  
  output$infantum_Dotplot <- renderPlot({
    
    ## Retrieve highlight data
    
    hl_df <- infantum_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_infantum <- infantum_MainDataTable1()
    
    ## Selection data
    
    jitter_data_infantum <- infantum_MainDataTable2()
    
    test_data_infantum <- infantum_TestDataTable2()
    
    # Dummy plot
    
    infantum_Dotplot <- ggplot(data = jitter_data_infantum, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_infantum) > 0) {
      
      y_pos_dot <- max(jitter_data_infantum$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
      
      # Retrieve line data  
        
      line_data_infantum <- jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
      
      # Format test data
      
      plot_test_data_infantum <- test_data_infantum[test_data_infantum$Majority.protein.IDs == HID[i],]
      plot_test_data_infantum <- plot_test_data_infantum[,c(4,5,3)]
      plot_test_data_infantum$p.adj <- format.pval(plot_test_data_infantum$p.adj, digits = 3)
        
      infantum_Dotplot <- ggdotplot(data = jitter_data_infantum[jitter_data_infantum$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                   position = position_jitter(0.2))+
        geom_line(data = line_data_infantum, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
        stat_pvalue_manual(plot_test_data_infantum, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
        scale_fill_brewer(palette="Dark2")+
        ylab("log2(LFQ)")+
        ggtitle(HID[i])+
        theme_minimal()+
        theme(plot.title = element_text(size=14,face = "bold"))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 14))+
        theme(legend.text = element_text(size = 12))+
        theme(legend.key.size = unit(1, 'cm'))
      
      infantum_Dotplot <- list(infantum_Dotplot)
      
      dotplot_list <- c(dotplot_list, infantum_Dotplot)
      
      }
      
      infantum_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(infantum_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$infantum_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$infantum_ui_datasetID
    exp_a <- input$infantum_TP
    exp_b <- input$infantum_referenceTP
    ID <- input$infantum_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$infantum_LFC)) 
    enrich_pvalue <- as.numeric(input$infantum_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == "All") {
      pg_quant <- full_infantum
    }else{
      pg_quant <- full_infantum %>% filter(Species == filter)
    }
    
    ## Change the column pater to user friendly reading (i. e., 0h instead of inft000).
    
    for (i in 1:length(hours)) {
      
      colnames(pg_quant) <- gsub(pattern = inf_sub[i], replacement = hours[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    volcanoplot <-
      ggplot(pg_quant, aes(x=difference_72h_0h, y=difference_72h_48h)) +
      theme_minimal()+
      ylab("-log10(p.value)") +
      xlab("log2(FoldChange)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    infantum_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
     
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log10(p.value)") +
        xlab("log2(FoldChange)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      infantum_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(infantum_Volcanoplot)
    
  })
  
  #### L. major ####
  
  #### Main table output  
  
  output$major_MainDataTable <- renderDataTable({
    major_MainDataTable1() %>%
      datatable(options = list(pageLength = 4, autoWidth = TRUE), 
                rownames = F) %>%
      formatRound('Value', digits = 2)
  })
  
  major_MainDataTable_proxy <- dataTableProxy('major_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$major_PCA <- renderPlot({
    
    # Load and format the data
    
    plot_df_major <- major_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_major))
    pca <- prcomp(t(na.omit(plot_df_major[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2 * -1
    scores$TimePoint <- factor(rep(c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"), each = 4), levels = c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"))
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    # Input plot
    
    major_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_point(size=3.5)+
      stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = TimePoint))+
      geom_text_repel(mapping = aes(label = Label),
                      size = 3,
                      fontface = 'bold',
                      color = 'black')+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_color_manual(values = colors_major)+
      scale_fill_manual(values = colors_major)+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))
    
    print(major_PCA)
    
  })
  
  output$major_HeatMap <- renderPlot({
    
    # Load and format the data
    
    plot_df_major <- major_FullDataTable1()
    lfq <- plot_df_major[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_major))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(hours)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = hours[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    # Input plot
    
    major_HeatMap <- gplots::heatmap.2(corr,
                                          trace="none",
                                          Colv=T,
                                          Rowv=T,
                                          dendrogram="column",
                                          density.info="density",
                                          srtCol=45,
                                          margins=c(10, 10), 
                                          key.xlab="Pearson's correlation coefficient",
                                          key.title="",
                                          keysize=1.5,
                                          col=hm_pal)
    
    print(major_HeatMap)
    
  })
  
  #### Clustering analysis output 
  
  output$major_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "major"
    
    pg_LFQ <- major_FullDataTable1()
    
    hl_df <- major_FullDataTable2()
    
    # User input variables
    
    xdim <- input$major_xDIM
    
    ydim <- input$major_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_major <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_major$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#EE2A7B"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_major$data[[1]][which(som_model_major$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    major_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_major$data[[1]][which(som_model_major$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
            
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      major_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(major_ClusterLinePlot)
    
  })
  
  output$major_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "major"
    
    pg_LFQ <- major_FullDataTable1()
    
    hl_df <- major_FullDataTable2()
    
    # User input variables
    
    xdim <- input$major_xDIM
    
    ydim <- input$major_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_major <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_major$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#EE2A7B"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_major$data[[1]][which(som_model_major$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LmjF.*", data$ID),
               sub(".*", "L. major", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      data <- data[data$TimePoint == "0h",]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    major_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. major" = "#EE2A7B",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    major_ClusterBarPlot <- ggplotly(major_ClusterBarPlot)
    
    print(major_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$major_Boxplot <- renderPlot({
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_major <- major_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_major <- major_MainDataTable2() 
    
    # Dummy plot
    
    major_Boxplot <-ggboxplot(data = plot_df_major, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_major)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- major_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_major) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_major <- jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        major_Boxplot <- ggboxplot(data = plot_df_major, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_major, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_major)+
          scale_color_brewer(palette="Dark2")+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        major_Boxplot <- list(major_Boxplot)
        
        boxplot_list <- c(boxplot_list, major_Boxplot)
        
      }
      
      major_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(major_Boxplot)
  })
  
  output$major_Dotplot <- renderPlot({
    
    ## Retrieve highlight data
    
    hl_df <- major_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_major <- major_MainDataTable1()
    
    ## Selection data
    
    jitter_data_major <- major_MainDataTable2()
    
    test_data_major <- major_TestDataTable2()
    
    # Dummy plot
    
    major_Dotplot <- ggplot(data = jitter_data_major, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_major) > 0) {
      
      y_pos_dot <- max(jitter_data_major$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_major <- jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        # Format test data
        
        plot_test_data_major <- test_data_major[test_data_major$Majority.protein.IDs == HID[i],]
        plot_test_data_major <- plot_test_data_major[,c(4,5,3)]
        plot_test_data_major$p.adj <- format.pval(plot_test_data_major$p.adj, digits = 3)
        
        major_Dotplot <- ggdotplot(data = jitter_data_major[jitter_data_major$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                      position = position_jitter(0.2))+
          geom_line(data = line_data_major, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
          stat_pvalue_manual(plot_test_data_major, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
          scale_fill_brewer(palette="Dark2")+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        major_Dotplot <- list(major_Dotplot)
        
        dotplot_list <- c(dotplot_list, major_Dotplot)
        
      }
      
      major_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(major_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$major_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$major_ui_datasetID
    exp_a <- input$major_TP
    exp_b <- input$major_referenceTP
    ID <- input$major_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$major_LFC)) 
    enrich_pvalue <- as.numeric(input$major_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == "All") {
      pg_quant <- full_major
    }else{
      pg_quant <- full_major %>% filter(Species == filter)
    }
    
    ## Change the column pater to user friendly reading (i. e., 0h instead of inft000).
    
    for (i in 1:length(hours)) {
      
      colnames(pg_quant) <- gsub(pattern = mjr_sub[i], replacement = hours[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    volcanoplot <-
      ggplot(pg_quant, aes(x=difference_72h_0h, y=difference_72h_48h)) +
      theme_minimal()+
      ylab("-log10(p.value)") +
      xlab("log2(FoldChange)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    major_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
      
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log10(p.value)") +
        xlab("log2(FoldChange)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      major_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(major_Volcanoplot)
    
  })
  
  
  #### L. mexicana ####
  
  #### Main table output  
  
  output$mexicana_MainDataTable <- renderDataTable({
    mexicana_MainDataTable1() %>%
      datatable(options = list(pageLength = 4, autoWidth = TRUE), 
                rownames = F) %>%
      formatRound('Value', digits = 2)
  })
  
  mexicana_MainDataTable_proxy <- dataTableProxy('mexicana_MainDataTable')
  
  #### Exploratory analysis output 
  
  output$mexicana_PCA <- renderPlot({
    
    # Load and format the data
    
    plot_df_mexicana <- mexicana_FullDataTable1()
    pca_columns <- grep(pattern = "^imputed.log2.LFQ.intensity.", colnames(plot_df_mexicana))
    pca <- prcomp(t(na.omit(plot_df_mexicana[,pca_columns])), scale.=TRUE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    scores <- data.frame(pca$x[,c("PC1","PC2")])
    scores$PC2 <- scores$PC2 * -1
    scores$TimePoint <- factor(rep(c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"), each = 4), levels = c("0h", "0.5h", "2h", "6h", "12h", "24h", "48h", "72h"))
    scores$Label <- gsub(".*_", replacement = "", x = gsub(pattern = "^imputed.log2.LFQ.intensity.",replacement = "",rownames(scores)))
    
    # Input plot
    
    mexicana_PCA <- ggplot(data=scores,aes(x=PC1, y=PC2, colour = TimePoint))+
      geom_point(size=3.5)+
      stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = TimePoint))+
      geom_text_repel(mapping = aes(label = Label),
                      size = 3,
                      fontface = 'bold',
                      color = 'black')+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      theme_minimal()+
      xlab(paste("PC1: ", round(percentVar[1] * 100,digits = 2), "%"))+
      ylab(paste("PC2: ", round(percentVar[2] * 100, digits = 2), "%"))+
      scale_color_manual(values = colors_mexicana)+
      scale_fill_manual(values = colors_mexicana)+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))
    
    print(mexicana_PCA)
    
  })
  
  output$mexicana_HeatMap <- renderPlot({
    
    # Load and format the data
    
    plot_df_mexicana <- mexicana_FullDataTable1()
    lfq <- plot_df_mexicana[, grep("^imputed.log2.LFQ.intensity.", names(plot_df_mexicana))]
    colnames(lfq) <- gsub("^imputed.log2.LFQ.intensity.", "", colnames(lfq))
    pattern <- unique(gsub(pattern = "_.*", replacement = "", colnames(lfq)))
    
    for (i in 1:length(hours)) {
      
      colnames(lfq) <- gsub(pattern = pattern[i], replacement = hours[i], colnames(lfq))
      
    }
    
    group <- as.factor(gsub("_[1-4]$", "", colnames(lfq)))
    lfq.avg <- t(apply(lfq, 1, function(x) tapply(x, group, mean,na.rm = T)))
    v <- apply(lfq, 1, sd)
    
    corr <- apply(lfq[rev(order(v)), ], 2, function(x) {
      apply(lfq[rev(order(v)), ], 2, function(y) {
        cor(x, y, use = "na.or.complete")
      })
    })
    
    # Input plot
    
    mexicana_HeatMap <- gplots::heatmap.2(corr,
                                          trace="none",
                                          Colv=T,
                                          Rowv=T,
                                          dendrogram="column",
                                          density.info="density",
                                          srtCol=45,
                                          margins=c(10, 10), 
                                          key.xlab="Pearson's correlation coefficient",
                                          key.title="",
                                          keysize=1.5,
                                          col=hm_pal)
    
    print(mexicana_HeatMap)
    
  })
  
  #### Clustering analysis output 
  
  output$mexicana_ClusterLinePlot <- renderPlot({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "mex"
    
    pg_LFQ <- mexicana_FullDataTable1()
    
    hl_df <- mexicana_FullDataTable2()
    
    # User input variables
    
    xdim <- input$mexicana_xDIM
    
    ydim <- input$mexicana_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_mex <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_mex$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    # This plot will always show, wether the user selects an input ID or not
    
    mycolors <- colorRampPalette(c("#000000", "#2B3990"))(length(nclust))
    
    plot_list <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_mex$data[[1]][which(som_model_mex$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
        geom_line(color="lightgrey")+
        stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
        scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
        ylab("z-score(LFQ)")+
        ggtitle(paste0("Cluster ",i))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(plot.title = element_text(size=22))+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      plot <- list(plot)
      
      plot_list <- c(plot_list, plot)
      
    }
    
    mexicana_ClusterLinePlot <- ggarrange(plotlist = plot_list)
    
    # Input plot
    
    # User dependant plot. The IDs the user select will show up in this plot, when present.
    
    if (nrow(hl_df) > 0) {
      
      HID <- unique(hl_df$Majority.protein.IDs)
      
      plot_list <- c()
      
      for (i in 1:length(nclust)) {
        
        data <- as.data.frame(som_model_mex$data[[1]][which(som_model_mex$unit.classif==i),])
        
        data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
        
        data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
        
        dataHL <- data[data$ID %in% HID,]
        
        if (nrow(dataHL) > 0) {
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          for (i in 1:length(HID)) {
            
            plot <- plot +
              geom_line(data=data[data$ID == HID[i],], aes(x=TimePoint, y=Value, group=ID,color = ID))+
              scale_color_brewer(name = "", palette = "Dark2")
            
            
          }
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }else{
          
          plot <- ggplot(data=data, aes(x=TimePoint, y=Value, group=ID)) +
            geom_line(color="lightgrey")+
            stat_summary(aes(y = Value,group=1,colour="Mean"), fun=mean, geom="line",group=1)+
            scale_color_manual(name = "", values = c("Mean" = mycolors[i]))+
            ylab("z-score(LFQ)")+
            ggtitle(paste0("Cluster ",i))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(plot.title = element_text(size=22))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
          plot <- list(plot)
          
          plot_list <- c(plot_list, plot)
          
        }
        
      }
      
      mexicana_ClusterLinePlot <- ggarrange(plotlist = plot_list)
      
    }
    
    print(mexicana_ClusterLinePlot)
    
  })
  
  output$mexicana_ClusterBarPlot <- renderPlotly({
    
    # Fixed variables
    
    set.seed(666)
    
    shortcut <- "mex"
    
    pg_LFQ <- mexicana_FullDataTable1()
    
    hl_df <- mexicana_FullDataTable2()
    
    # User input variables
    
    xdim <- input$mexicana_xDIM
    
    ydim <- input$mexicana_yDIM
    
    # xdim <- 2 
    # 
    # ydim <- 3
    
    # Script
    
    # For the SOM cluster analysis and  plots we have to standardize (compute de z-score) the data.
    
    # Select columns we need:
    
    mean_LFQ <- dplyr::select(pg_LFQ, c(Protein.IDs, Majority.protein.IDs, contains("mean2_naimp_")))
    
    pg_clust <- dplyr::select(mean_LFQ, 2,contains(shortcut))
    
    # Filter proteins by IQR
    
    pg_clust$max <- pg_clust %>% dplyr::select(contains('mean2_naimp_')) %>% apply(1, max)
    
    pg_clust <- pg_clust[which(pg_clust$max>quantile(pg_clust$max,0.25)),] %>% subset(select=-max)
    
    iqrs <- apply(pg_clust[,2:ncol(pg_clust)],1,IQR)
    
    pg_clust <-pg_clust[which(iqrs>quantile(iqrs, prob = c(0.10))),]
    
    # Data tidy
    
    colnames(pg_clust) <- c("Majority.Protein.IDs",paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "") )
    
    # Compute the z-score
    
    pg_clust[,2:ncol(pg_clust)] <- t(apply(pg_clust[,2:ncol(pg_clust)],1,function(x) effectsize::standardize(as.numeric(x),normalize=FALSE)))
    
    # Data tidy
    
    prep.pg_clust <- data.frame(as.matrix(pg_clust[,2:ncol(pg_clust)]))
    
    colnames(prep.pg_clust) <- paste(c("0","0.5","2","6","12","24",'48',"72"), 'h', sep = "")
    
    rownames(prep.pg_clust) <- pg_clust$Majority.Protein.IDs
    
    # Compute the SOM model
    
    som_grid <- somgrid(xdim = xdim, ydim = ydim, topo="hexagonal") # 3x3=9 clusters
    
    som_model_mex <- som(as.matrix(prep.pg_clust), 
                         grid=som_grid, 
                         rlen=1000, 
                         alpha=c(0.05,0.01))    
    
    nclust <- unique(som_model_mex$unit.classif)
    
    # Plot the clusters
    
    # Dummy plot:
    
    mycolors <- colorRampPalette(c("#000000", "#2B3990"))(length(nclust))
    
    plot_data <- c()
    
    for (i in 1:length(nclust)) {
      
      data <- as.data.frame(som_model_mex$data[[1]][which(som_model_mex$unit.classif==i),])
      
      data <- data %>% rownames_to_column('ID') %>% pivot_longer(-ID, names_to = "TimePoint", values_to = "Value")
      
      data$Species <- "Unknown"
      
      data$Species <-
        ifelse(grepl("ENSMUS.*", data$ID),
               sub(".*", "M. musculus", data$Species),
               sub(".*", "Unknown", data$Species))
      
      data$Species <-
        ifelse(grepl("LmxM.*", data$ID),
               sub(".*", "L. mexicana", data$Species),
               sub("", "", data$Species))
      
      data$TimePoint <- factor(data$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
      
      data <- data[data$TimePoint == "0h",]
      
      data <- data %>%
        group_by(Species) %>%
        summarize(Count = n())%>% 
        mutate(perc = Count/sum(Count))
      
      data$Cluster <- as.character(i)
      
      plot_data <- rbind(plot_data, data)
      
    }
    
    mexicana_ClusterBarPlot <- ggplot(plot_data, aes(fill=Species, y = Count , x= Cluster)) + 
      geom_bar(position="stack", stat="identity")+
      scale_fill_manual(values = c("L. mexicana" = "#2B3990",
                                   "M. musculus" = "#b4b4b4"))+
      theme_minimal()+
      ylab("Protein IDs")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.title = element_text(size = 16))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    mexicana_ClusterBarPlot <- ggplotly(mexicana_ClusterBarPlot)
    
    print(mexicana_ClusterBarPlot)
    
  })
  
  #### Pairwise analysis  time course output 
  
  output$mexicana_Boxplot <- renderPlot({
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_mexicana <- mexicana_MainDataTable1() 
    
    ## Selection data
    
    jitter_data_mexicana <- mexicana_MainDataTable2() 
    
    # Dummy plot
    
    mexicana_Boxplot <-ggboxplot(data = plot_df_mexicana, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
      scale_fill_manual(values = colors_mexicana)+
      scale_color_brewer(palette="Dark2")+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.text = element_text(size = 10))+
      theme(legend.key.size = unit(1, 'cm'))
    
    # Input plot
    
    hl_df <- mexicana_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    if (nrow(jitter_data_mexicana) > 0) {
      
      boxplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_mexicana <- jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        mexicana_Boxplot <- ggboxplot(data = plot_df_mexicana, x="TimePoint", y="Value", fill = "TimePoint",alpha = 0.3) +
          geom_jitter(data = jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],],aes(x=TimePoint, y=Value, color = Replica),width = 0.20)+
          geom_line(data = line_data_mexicana, aes(x=TimePoint, y=Value,group = 1))+
          scale_fill_manual(values = colors_mexicana)+
          scale_color_brewer(palette="Dark2")+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        mexicana_Boxplot <- list(mexicana_Boxplot)
        
        boxplot_list <- c(boxplot_list, mexicana_Boxplot)
        
      }
      
      mexicana_Boxplot <- ggarrange(plotlist = boxplot_list)
      
    }
    
    print(mexicana_Boxplot)
  })
  
  output$mexicana_Dotplot <- renderPlot({
    
    ## Retrieve highlight data
    
    hl_df <- mexicana_FullDataTable2()
    
    HID <- unique(hl_df$Majority.protein.IDs)
    
    # Load  and format the plot data
    
    ## Whole data
    
    plot_df_mexicana <- mexicana_MainDataTable1()
    
    ## Selection data
    
    jitter_data_mexicana <- mexicana_MainDataTable2()
    
    test_data_mexicana <- mexicana_TestDataTable2()
    
    # Dummy plot
    
    mexicana_Dotplot <- ggplot(data = jitter_data_mexicana, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    # Input plot
    
    if (nrow(jitter_data_mexicana) > 0) {
      
      y_pos_dot <- max(jitter_data_mexicana$Value)+0.2
      
      dotplot_list <- c()
      
      for (i in 1:length(HID)) {
        
        # Retrieve line data  
        
        line_data_mexicana <- jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],] %>%
          group_by(TimePoint) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE)) 
        
        # Format test data
        
        plot_test_data_mexicana <- test_data_mexicana[test_data_mexicana$Majority.protein.IDs == HID[i],]
        plot_test_data_mexicana <- plot_test_data_mexicana[,c(4,5,3)]
        plot_test_data_mexicana$p.adj <- format.pval(plot_test_data_mexicana$p.adj, digits = 3)
        
        mexicana_Dotplot <- ggdotplot(data = jitter_data_mexicana[jitter_data_mexicana$Majority.protein.IDs == HID[i],], x="TimePoint", y="Value", fill = "Replica",
                                      position = position_jitter(0.2))+
          geom_line(data = line_data_mexicana, aes(x=TimePoint, y=Value,group = 1),alpha = 0.8)+
          stat_pvalue_manual(plot_test_data_mexicana, y.position = y_pos_dot, step.increase = 0.1,label = "p.adj")+
          scale_fill_brewer(palette="Dark2")+
          ylab("log2(LFQ)")+
          ggtitle(HID[i])+
          theme_minimal()+
          theme(plot.title = element_text(size=14,face = "bold"))+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          theme(axis.text=element_text(size=18),
                axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 14))+
          theme(legend.text = element_text(size = 12))+
          theme(legend.key.size = unit(1, 'cm'))
        
        mexicana_Dotplot <- list(mexicana_Dotplot)
        
        dotplot_list <- c(dotplot_list, mexicana_Dotplot)
        
      }
      
      mexicana_Dotplot <- ggarrange(plotlist = dotplot_list)
      
    }
    
    print(mexicana_Dotplot)
  })
  
  #### Pairwise analysis  reference output 
  
  output$mexicana_Volcanoplot <- renderPlotly({
    
    #Input variables
    
    filter <- input$mexicana_ui_datasetID
    exp_a <- input$mexicana_TP
    exp_b <- input$mexicana_referenceTP
    ID <- input$mexicana_ui_selectedIDs
    enrich_s0 <- log2(2^as.numeric(input$mexicana_LFC)) 
    enrich_pvalue <- as.numeric(input$mexicana_Pvalue)
    
    # Fixed variables
    
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    val_count_a <- paste0('value_count_',exp_a) 
    val_count_b <- paste0('value_count_',exp_b)
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    
    min_quant_events <- 0 
    replicate_count <- 4
    enrich_c <- .05 
    volcano_threshold_linetype <- 2
    volcano_threshold_linesize <- 0.5
    difference_volcano_column <- 'mean2_naimp_meas_'
    volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
    
    # Dummy plot
    
    ## Filter the data frame for data set user input.
    
    if (filter == "All") {
      pg_quant <- full_mexicana
    }else{
      pg_quant <- full_mexicana %>% filter(Species == filter)
    }
    
    ## Change the column pater to user friendly reading (i. e., 0h instead of inft000).
    
    for (i in 1:length(hours)) {
      
      colnames(pg_quant) <- gsub(pattern = mex_sub[i], replacement = hours[i], colnames(pg_quant))
      
    }
    
    ## Dummy plot to show until user selection.
    
    volcanoplot <-
      ggplot(pg_quant, aes(x=difference_72h_0h, y=difference_72h_48h)) +
      theme_minimal()+
      ylab("-log10(p.value)") +
      xlab("log2(FoldChange)") +
      theme(legend.text = element_text(size = 12))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(plot.title = element_text(size=20))+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"))
    
    mexicana_Volcanoplot <-
      plotly::ggplotly(volcanoplot) %>%
      plotly::layout(showlegend = FALSE)
    
    # Input plot
    
    ## Upon proper user selection (if) produce the volcano plot.
    
    if (difference %in% colnames(pg_quant)) {
      
      ## Highlight user input, if any.
      
      pg_quant$highlight <- FALSE
      
      if (length(ID) != 0) {
        
        ids_searchstring <- paste(ID, collapse='|')
        
        matching_rows <- grep(ids_searchstring, as.character(pg_quant$Protein.IDs))
        
        pg_quant$highlight[matching_rows] <- TRUE
        
      }
      
      
      ## Filter count values. Standard set to 0, so they can check the dot all along the time course.
      ## In Mario's standard script is usually set to 2. But in this case, since it shown in the 
      ## box plot and in the dot plot, would be confusing not to show it here as well.
      ## Can be adjusted on the fixed variables section of this function.
      
      filter_min_quant_events <-
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      
      pg_quant <- pg_quant[filter_min_quant_events,]
      
      ## Select and highlight the enriched dots.
      
      my_enriched <- Enriched(pg_quant[c(difference, pvalue)],
                              c=enrich_c,
                              s0=enrich_s0,
                              pvalue=enrich_pvalue,
                              linetype=volcano_threshold_linetype,
                              size=volcano_threshold_linesize)
      pg_quant$enriched <- as.factor(my_enriched$enriched)
      
      
      my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
      my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
      
      ## Generate a user friendly plotly string for the hoover.
      
      plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Protein.IDs: %s<br />my_label: %s',
                                      pg_quant[[difference]], pg_quant[[pvalue]],
                                      sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                      pg_quant$my_label),
                              switch('Gene.names' %in% names(pg_quant),
                                     sprintf('<br />Gene.names: %s', pg_quant$Gene.names),
                                     NULL),
                              switch('Protein.names' %in% names(pg_quant),
                                     sprintf('<br />Protein.names: %s', pg_quant$Protein.names),
                                     NULL),
                              switch('description' %in% names(pg_quant),
                                     sprintf('<br />description: %s', pg_quant$description),
                                     NULL))
      
      pg_quant$my_text <- plotly_string
      
      ## Subset the data for the volcano plot.
      
      pg_quant_subset <-
        pg_quant %>%
        rowwise() %>%
        mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                   min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                                 'measured',
                                 ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                        'some imputed',
                                        'ref imputed'))) %>%
        select(Protein.IDs, my_label, matches(difference), matches(pvalue),
               my_text, highlight, enriched, imputation)
      
      ## Plot
      
      ## Select maximum x and y values of your volcano plot.
      
      max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]),
                         na.rm=TRUE)
      max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                           10^enrich_pvalue + 0.5))
      
      volcanoplot <-
        ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +
        geom_hline(yintercept=0, color='grey80', na.rm=TRUE) +
        geom_vline(xintercept=0, color='grey80', na.rm=TRUE) +
        plot(my_enriched) +
        geom_point(data=subset(pg_quant_subset, !highlight),
                   aes(alpha=my_points_nohighlight,
                       color=my_points_nohighlight,
                       fill=my_points_nohighlight,
                       size=my_points_nohighlight),
                   na.rm=TRUE) +
        geom_point(data=subset(pg_quant_subset, highlight),
                   aes(alpha=my_points_highlight,
                       color=my_points_highlight,
                       fill=my_points_highlight,
                       size=my_points_highlight),
                   na.rm=TRUE) +
        scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                           guide="none") +
        scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                           guide="none") +
        scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                          guide="none") +
        scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                          guide="none")+
        coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
        ylab("-log10(p.value)") +
        xlab("log2(FoldChange)") +
        annotate('text', label=exp_a, x=max_x_value / 2, y=-.4) +
        annotate('text', label=exp_b, x=-max_x_value / 2, y=-.4)+
        theme_minimal()+
        theme(legend.text = element_text(size = 14))+
        theme(legend.key.size = unit(1, 'cm'))+
        theme(legend.position="top")+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
      
      mexicana_Volcanoplot <-
        plotly::ggplotly(volcanoplot,tooltip='text') %>%
        plotly::layout(showlegend = FALSE)
      
    }
    
    print(mexicana_Volcanoplot)
    
  })
  
  #### Orthology ####
  
  output$Ortho_table <- renderDataTable({orth_df})
  
  ortho_MainDataTable <- reactive({
    ortho_main_table <- orth_df 
    ortho_main_table
  })
  
  ortho_MainDataTable_proxy <- dataTableProxy('Ortho_table')
  
  output$OrthoPlot <- renderPlot({
    
    OrthoPlot <- ggplot(data = df_infantum, aes(TimePoint, Value))+
      ylab("log2(LFQ)")+
      xlab("TimePoint")+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    
    id_input <- input$ui_orthoSelectedIDs
    
    # id_input <- c("OrthoGroup_7610","OrthoGroup_7832")
      
    id <- c()
    
    for (i in 1:length(id_input)) {
      
      loop_id <- orth_df[orth_df$OrthoGroup == id_input[i],]
      
      loop_leish <- loop_id$L_infantum.fasta
      
      if (length(loop_leish) > 0) {
        
        if (loop_leish == "*") {
          
          loop_leish <- loop_id$L_major.fasta
          
          }
        
        id <- c(id, loop_leish)
        
        }
    }
    
    id_pattern <- paste0(id, collapse = "|")
    
    df_f <- data.frame()
    
    if (length(id) > 0) {
      
      if (grepl(pattern = "LINF", x = id_pattern)) {
        
        df_infantum_loop <- df_infantum %>% filter(Majority.protein.IDs %in% id)
        
        if (nrow(df_infantum_loop) > 0) {
          
          df_f_inf <- orth_df[grep(pattern = id_pattern, orth_df$L_infantum.fasta),]
          
          df_f <- rbind(df_f, df_f_inf)
        }
      }
      
      if (grepl(pattern = "LmjF", x = id_pattern)) {
        
        df_major_loop <- df_major %>% filter(Majority.protein.IDs %in% id)
        
        if (nrow(df_major_loop) > 0) {
          
          df_f_mjr <- orth_df[grep(pattern = id_pattern, orth_df$L_major.fasta),]
          
          df_f <- rbind(df_f, df_f_mjr)
        }
        
      }
      
      if (grepl(pattern = "LmxM", x = id_pattern)) {
        
        df_mexicana_loop <- df_mexicana %>% filter(Majority.protein.IDs %in% id)
        
        if (nrow(df_mexicana_loop) > 0) {
          
          df_f_mex <- orth_df[grep(pattern = id_pattern, orth_df$L_mexicana.fasta),]
          
          df_f <- rbind(df_f, df_f_mex)
        }
        
      }
      
      
      
      if (nrow(df_f) > 0) {
        
        infantum_id <- unlist(strsplit(df_f$L_infantum.fasta,","))
        
        major_id <- unlist(strsplit(df_f$L_major.fasta,","))
        
        mexicana_id <- unlist(strsplit(df_f$L_mexicana.fasta,","))
        
        plot_df <- c()
        
        df_infantum <- df_infantum
        
          for (i in 1:length(infantum_id)) {
            
            if (infantum_id[i] != "*") {
            
            infantum_loop <- df_infantum[grep(pattern = infantum_id[i], df_infantum$Majority.protein.IDs),]
            
            plot_df <- rbind(plot_df,infantum_loop)
            
            }
            
          }
        
        df_major <- df_major
          
          for (i in 1:length(major_id)) {
            
            if (major_id[i] != "*") {
            
            major_loop <- df_major[grep(pattern = major_id[i], df_major$Majority.protein.IDs),]
            
            plot_df <- rbind(plot_df,major_loop)
            
            }
            
          }
          
      
        
        df_mexicana <- df_mexicana
        
          for (i in 1:length(mexicana_id)) {
            
            if (mexicana_id[i] != "*") {
              
            mexicana_loop <- df_mexicana[grep(pattern = mexicana_id[i], df_mexicana$Majority.protein.IDs),]
            
            plot_df <- rbind(plot_df,mexicana_loop)
            
            }
          }
          
        plot_df <- plot_df %>%
          group_by(TimePoint,Species, Majority.protein.IDs) %>%
          dplyr::summarize(Value = mean(Value, na.rm=TRUE))
        
        if (nrow(plot_df) > 0) {
          
          plot_df$TimePoint <- factor(plot_df$TimePoint, levels = c("0h","0.5h","2h","6h","12h","24h",'48h',"72h"))
          
          infantum_id_number <- length(unique(plot_df[plot_df$Species == "L. infantum",]$Majority.protein.IDs))
          
          major_id_number <- length(unique(plot_df[plot_df$Species == "L. major",]$Majority.protein.IDs))
          
          mexicana_id_number <- length(unique(plot_df[plot_df$Species == "L. mexicana",]$Majority.protein.IDs))
          
          colors_infantum <- c()
          
          if (infantum_id_number == 0) {
            
            colors_infantum <- c()
            
          }else{
            
            if (infantum_id_number == 1) {
              
              colors_infantum <- "#00A651"
              
            }else{
              colors_infantum <- colorRampPalette(c("#004722", "#00A651"))(infantum_id_number)
              
            }
            
          }
          
          colors_major <- c()
          
          if (major_id_number == 0) {
            
            colors_major <- c()
            
          }else{
            
            if (major_id_number == 1) {
              
              colors_major <- "#EE2A7B"
              
            }else{
              colors_major <- colorRampPalette(c("#651234", "#EE2A7B"))(major_id_number)
              
            }
            
          }
          
          colors_mexicana <- c()
          
          if (mexicana_id_number == 0) {
            
            colors_mexicana <- c()
            
          }else{
            
            if (mexicana_id_number == 1) {
              
              colors_mexicana <- "#2B3990"
              
            }else{
              colors_mexicana <- colorRampPalette(c("#12183D", "#2B3990"))(mexicana_id_number)
              
            }
            
          }
          
          colors <- c(colors_infantum, colors_major, colors_mexicana)
          
          OrthoPlot <- ggplot(data=plot_df, aes(x=TimePoint, y=Value, group=Majority.protein.IDs,color = Majority.protein.IDs)) +
            geom_line()+
            scale_color_manual("",values = colors)+
            ylab("log2(LFQ)")+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            theme(legend.text = element_text(size = 14))+
            theme(legend.key.size = unit(1, 'cm'))+
            theme(legend.position="top")+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
          
        }
        
      }
      
    }
    
    print(OrthoPlot)
    
  })
  
  #### Observe functions ####
  
  ### L. infantum 
  
  observe({
    updateSelectizeInput(
      session,
      'infantum_ui_selectedIDs',
      choices=infantum_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$infantum_MainDataTable_rows_selected)],
      selected=infantum_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$infantum_MainDataTable_rows_selected)]
    )
  })
  
  observeEvent(input$update_infantum_MainDataTable, {
    rows <-
      match(input$infantum_ui_selectedIDs,
            infantum_MainDataTable1()[['Majority.protein.IDs']])
    selectRows(infantum_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### L. major 
  
  observe({
    updateSelectizeInput(
      session,
      'major_ui_selectedIDs',
      choices=major_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$major_MainDataTable_rows_selected)],
      selected=major_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$major_MainDataTable_rows_selected)]
    )
  })
  
  observeEvent(input$update_major_MainDataTable, {
    rows <-
      match(input$major_ui_selectedIDs,
            major_MainDataTable1()[['Majority.protein.IDs']])
    selectRows(major_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### L. mexicana 
  
  observe({
    updateSelectizeInput(
      session,
      'mexicana_ui_selectedIDs',
      choices=mexicana_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$mexicana_MainDataTable_rows_selected)],
      selected=mexicana_MainDataTable1()[['Majority.protein.IDs']][as.numeric(input$mexicana_MainDataTable_rows_selected)]
    )
  })
  
  observeEvent(input$update_mexicana_MainDataTable, {
    rows <-
      match(input$mexicana_ui_selectedIDs,
            mexicana_MainDataTable1()[['Majority.protein.IDs']])
    selectRows(mexicana_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  ### Orthology
  
  observe({
    updateSelectizeInput(
      session,
      'ui_orthoSelectedIDs',
      choices=ortho_MainDataTable()[['OrthoGroup']][as.numeric(input$Ortho_table_rows_selected)],
      selected=ortho_MainDataTable()[['OrthoGroup']][as.numeric(input$Ortho_table_rows_selected)]
    )
  })
  
  observeEvent(input$update_ortho_MainDataTable, {
    rows <-
      match(input$ui_orthoSelectedIDs,
            ortho_MainDataTable()[['OrthoGroup']])
    selectRows(ortho_MainDataTable_proxy,
               as.numeric(rows))
  })
  
  
}
