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

ui <-  navbarPage("Leishmania Infectome",
            tabPanel("Welcome",
                     mainPanel(
                       h1("Abstract")
                     ),
                     mainPanel(
                       h4("This is a random text to show how the abstract for the project would look like"),
                       p("A 37-year-old Toru Watanabe has just arrived in Hamburg, West Germany. When he hears an orchestral cover of the Beatles' 
                       song Norwegian Wood, he is suddenly overwhelmed by feelings of loss and nostalgia. He thinks back to the 1960s, when so much 
                       happened that touched his life. Watanabe, his classmate Kizuki, and Kizuki's girlfriend Naoko are the best of friends. Kizuki
                       and Naoko are particularly close and feel as if they are soulmates, and Watanabe seems more than happy to be their enforcer. 
                       This idyllic existence is shattered by the unexpected suicide of Kizuki on his 17th birthday. Kizuki's death deeply touches 
                       both surviving friends; Watanabe feels the influence of death everywhere, while Naoko feels as if some integral part of her
                       has been permanently lost. The two of them spend more and more time together going for long walks on Sundays, although 
                       feelings for each other are never clarified in this interval. On the night of Naoko's 20th birthday, she feels especially 
                       vulnerable and they have sex, during which Watanabe realizes that she is a virgin. Afterward, Naoko leaves Watanabe a letter 
                       saying that she needs some time apart and is quitting college to go to a sanatorium. These events are set against a backdrop 
                       of civil unrest. The students at Watanabe's college go on strike and call for a revolution. Inexplicably, the students end 
                       their strike and act as if nothing had happened, which enrages Watanabe as a sign of hypocrisy. Watanabe is befriended by a 
                       fellow drama classmate, Midori Kobayashi. She is everything that Naoko is not? outgoing, vivacious, and supremely 
                       self-confident. Despite his love for Naoko, Watanabe finds himself attracted to Midori as well. Midori reciprocates his 
                       feelings, and their friendship grows during Naoko's absence. Watanabe and Midori share a special kind of relationship where 
                       both of them understand each other. Watanabe visits Naoko at her secluded mountain sanatorium near Kyoto. There he meets 
                       Reiko Ishida, an older patient there who has become Naoko's confidante. During this and subsequent visits, Reiko and Naoko 
                       reveal more about their past: Reiko talks about the cause of her downfall into mental illness and details the failure of her
                       marriage, while Naoko talks about the unexpected suicide of her older sister several years ago. When he returns to Tokyo,
                       Watanabe unintentionally alienates Midori through both his lack of consideration of her wants and needs, and his continuing
                       thoughts about Naoko. He writes a letter to Reiko, asking for her advice about his conflicted affections for both Naoko and
                       Midori. He does not want to hurt Naoko, but he does not want to lose Midori either. Reiko counsels him to seize this chance
                       for happiness and see how his relationship with Midori turns out. A later letter informs Watanabe that Naoko has killed 
                         herself. Watanabe, grieving and in a daze, wanders aimlessly around Japan, while Midori?with whom he hasn't kept in 
                         touch?wonders what has happened to him. After about a month of wandering, he returns to the Tokyo area and gets in contact
                         with Reiko, who leaves the sanatorium to come to visit. Reiko stays with Watanabe, and they have sex. It is through this 
                         experience and the intimate conversation that Watanabe and Reiko share that night, that he comes to realize that Midori is 
                         the most important person in his life. After he sees Reiko off, Watanabe calls Midori to declare his love for her. Midori
                         asks, Where are you now?, and the novel ends with Watanabe pondering that question. "),
                       h2("Experimental design"),
                       p("I think it will be a good idea to have a constant general overview on the experimentall setup or a summary of the database
                                        goal. Thus we can use this paragraph for a brief explanantion. Also the following two images can be changed to whatever
                                        we feel like."),
                       img(src = "ExperimentalSetup.png", height = 250, width = 500),
                       h2(""),
                       h2("References"),
                       p("Any references we deem necessary"),
                       h4("Bioinforamtic analysis tools:"),
                       p("Lechner, M., Findeiß, S., Steiner, L. et al. Proteinortho: Detection of (Co-)orthologs in large-scale analysis.
                         BMC Bioinformatics 12, 124 (2011). https://doi.org/10.1186/1471-2105-12-124"),
                       p("Lechner M, Hernandez-Rosales M, Doerr D, Wieseke N, Thévenin A, Stoye J, et al. (2014) Orthology Detection Combining 
                         Clustering and Synteny for Very Large Datasets. PLoS ONE 9(8): e105015. https://doi.org/10.1371/journal.pone.0105015"),
                       p("Plotly Technologies Inc. Collaborative data science. Plotly Technologies Inc. Montréal, QC. (2015) URL: https://plot.ly"),
                       p("R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical 
                       Computing, Vienna, Austria. URL https://www.R-project.org/."),
                       p("Wehrens R, Kruisselbrink J (2018). Flexible Self-Organizing Maps in kohonen 3.0.Journal of Statistical Software, 87(7), 
                         118. doi: 10.18637/jss.v087.i07"),
                       p("Wehrens R, Buydens LMC (2007). Self- and Super-Organizing Maps in R: The kohonen Package.
                         Journal of Statistical Software, 21(5), 1?19. doi: 10.18637/jss.v021.i05"),
                       p("Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. 
                         ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org. "),
                     )),
             tabPanel("L. infantum",
                      fluidRow(column(2,
                                      h2("Settings"),
                                      selectInput(inputId = "infantum_ui_datasetID", label = h4("Select a data set:"), choices = c("All",
                                                                                                                          unique(df_infantum$Species)))),
                               column(6,
                                      h2("Quantified proteins"),
                                      h2(""),
                                      dataTableOutput("infantum_MainDataTable")),
                               column(4,
                                      h2("Selected proteins:"),
                                      selectizeInput('infantum_ui_selectedIDs',
                                                     NULL,
                                                     choices=NULL,
                                                     multiple=TRUE),
                                      actionButton("update_infantum_MainDataTable",
                                                   "Update selection"))),
                                      # textInput(inputId = "ui_proteinID", label = h4("Protein ID input:"), value = "Enter text..."))),
                      navbarPage("Proteome Analysis",
                                 tabPanel("Exploratory Analysis",
                                          # h2("Sample correlation"),
                                          # p("To summarize the main characteristics of the leishmania infectome data set we perform an exploratory
                                          #   analysis, mainly focused on sample correlation. This way, we can assess for techincal sample 
                                          #   variance (replicates correlation) and experimental sample variance (groups, in this case,
                                          #   timepoints) in a graphical way. This analysis can be performed to the full data set or by 
                                          #   separating the data per species."),
                                          # selectInput(inputId = "infantum_datasetID", label = h4("Select a data set:"), choices = c("All",
                                          #                                                                                  unique(df_infantum$Species))),
                                          fluidRow(column(6,
                                                          h4("PCA plot"),
                                                          # p("Scatter plot depicting the output of a principal componetns analysis (PCA). Overall, 
                                                          #   replicates are expected to show relatively small distances between them, while sample 
                                                          #   groups often show higher distances between them, depending on the magnitude of 
                                                          #   intensity changes."),
                                                          plotOutput("infantum_PCA",height = "1000px")),
                                                   column(6,
                                                          h4("HeatMap plot"),
                                                          # p("Heat map depicting the pearson's correlation coefficient and the clustering between samples. Overall, 
                                                          #   replicates are expected to show high correlation between them, while sample 
                                                          #   groups might show lower correlations, depending on the magnitude of intensity changes."),
                                                          plotOutput("infantum_HeatMap",height = "1000px")))),
                                 tabPanel("Clustering analysis",
                                           # h2("Self organizing maps (SOMs)"),
                                           # p("Self Organizing Maps (SOMs) are a tool for visualizing patterns in high dimensional data by producing a 
                                           #   2 dimensional representation, which displays meaningful patterns in the higher dimensional structure. It is
                                           #   an unsupervised machine learning technique that preserves the topological structure of the data.The size of
                                           #   the map grid is user defined and is the result of the multiplication of X and Y dimensions. This will result
                                           #   in n clusters in which data is expect to behave in a similar pattern. Default X and Y dimension are set to
                                           #   2 and 3 respectively, but should be adjusted per data set. The first adjustment criteria is trying to maximize
                                           #   the number of clusters with an homogenous compositions (individuals close to the mean). The second criteria is
                                           #   not having clusters with a very low count of individuals. This analysis can be performed to the full data set or by 
                                           #   separating the data per species."),
                                           fluidRow(
                                                    # column(2,
                                                    #        selectInput(inputId = "infantum_datasetID2", label = h4("Select a data set:"), choices = c("All",
                                                    #                                                                        unique(df_infantum$Species)))),
                                                    # column(2,
                                                    #        textInput(inputId = "infantum_proteinID2", label = h4("Protein ID input:"), value = "Enter text...")),
                                                    column(4,
                                                           bsCollapsePanel("Advanced settings",
                                                                    column(6,
                                                                           sliderInput("infantum_xDIM", "Select X dimensions:",
                                                                                       min = 1, max = 3, 
                                                                                       value = 2, step = 1)),
                                                                    column(6,
                                                                           sliderInput("infantum_yDIM", "Select Y dimensions:",
                                                                                       min = 1, max = 3,
                                                                                       value = 3, step = 1))))),
                                           fluidRow(column(8,
                                                           h4("Lineplot: Resulting clusters and their mean"),
                                                           # p("Line plot representation of the resulting n clusters. For each cluster's protein ID,
                                                           #   the standarized LFQ intensity (z-score(LFQ), y-axis) per timepoint (x-axis) is 
                                                           #   depicted in lightgrey. The cluster mean is depicted in each species color gradient.
                                                           #   Finally, and upon Protein ID selection, the selected target is highlighted in orange."),
                                                           plotOutput("infantum_ClusterLinePlot",height = "1000px")),
                                                    column(4,
                                                           h4("Barplot: Protein count per cluster"),
                                                           # p("Barplot showing the protein in count and species distribution per cluster"),
                                                           plotlyOutput("infantum_ClusterBarPlot",height = "1000px")))),
                                 tabPanel("Pairwise Analysis - Time course",
                                          # h2("Time course comparison"),
                                          # p("The time course paiwise analysis results are depicted in this section. Thus, each timepoint is statistically
                                          #   compared to the one that precedes it.The analysis is performed to the full data set, but can
                                          #   be aslo visualized per species. The data distribution per timpoint is depicted in a boxplot shape.
                                          #   Upon protein ID input, its behaviour on the data distribution, as well as its statistical significance (line plot),
                                          #   is depicted. Finally, we can directly compare data between time course time points on a volcano plot."),
                                          # fluidRow(column(2,
                                          #                 selectInput(inputId = "infantum_speciesID", label = h4("Select data set:"), choices = c("All",
                                          #                                                                                                         unique(df_infantum$Species)))),
                                          #          column(2,
                                          #                 textInput(inputId = "infantum_proteinID", label = h4("Protein ID input:"), value = "Enter text...")),
                                          #          column(8)),
                                          fluidRow(column(6,
                                                          h4("Boxplot: Proteome distribution along the time course"),
                                                          # p("The protein IDs distribution is depicted as a boxplot. Thus, we have an indivdual boxplot per
                                                          #     each timepoint (x-axis) showing the proteins intensity, as log2(LFQ) (y-axis), distribution. Upon
                                                          #     protein ID input and when it is present on the data set, each ID replicate's intensity and mean line 
                                                          #     is depicted on the plot."),
                                                          plotOutput("infantum_Boxplot",height = "1000px")),
                                                   column(6,
                                                          h4("Dotplot: Selected ID statiscal significance"),
                                                          # p("Upon protein ID input and when present on the data set, a dotpot is generated. Thus, 
                                                          #   each replicate intensity, as log2(LFQ) (y-axis), per each timepoint (x-axis) is represented.
                                                          #   Additionally, a Welch's t-test is performed between adjacent timepoints whose results are plotted 
                                                          #   as adjusted (BH-FDR) p-values."),
                                                          plotOutput("infantum_Dotplot",height = "1000px")))),
                                 tabPanel("Pairwise Analysis - Reference",
                                          # h2("Against reference comparison"),
                                          h4("Volcano plot: Individual IDs fold change and p-value"),
                                          # p("Volcano plot showing the LFQ intensity difference (log2(FoldChange), x-axis) and its significance (-log10(p-value),
                                          #   y-axis) for each protein ID. Thus each ID is represented by a light grey dot, when not significant, or a by a black
                                          #   dot, when significant. Default significance threshold is set to p-value < 0.05 and abs(log2(FoldChange)) > 2, but
                                          #   can be ajusted to user needs. Finally, and upon ID selection, the dot is higlighted in orange, when not significant, 
                                          #   or in red, when significant."),
                                          fluidRow(column(2,
                                                         selectInput(inputId = "infantum_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_infantum$TimePoint))),
                                                   column(2,
                                                         selectInput(inputId = "infantum_TP", label = h4("Select a time point:"), choices = unique(df_infantum$TimePoint))),
                                                   column(4,
                                                          bsCollapsePanel("Advanced settings",
                                                                   column(6,
                                                                          selectInput(inputId = "infantum_Pvalue", label = h4("Select a p-value:"), choices = c(0.05,0.01))),
                                                                   column(6,
                                                                          sliderInput("infantum_LFC", "Select a log2(FoldChange):",
                                                                                      min = 0, max = 5,
                                                                                      value = 1, step = 0.1))))),
                                          plotlyOutput("infantum_Volcanoplot",height = "1000px")))),
                                                   # column(4,
                                                   #        h4("Volcano plot: Individual IDs fold change and p-value"),
                                                   #        p("blablabla volcanoplot"),
                                                   #        fluidRow(column(6,
                                                   #                        selectInput(inputId = "infantum_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_infantum$TimePoint))),
                                                   #                 column(6,
                                                   #                        selectInput(inputId = "infantum_TP", label = h4("Select a time point:"), choices = unique(df_infantum$TimePoint)))),
                                                   #        plotlyOutput("Plot_infantum_8",height = "1000px")))))),
            tabPanel("L. major",
                     fluidRow(column(2,
                                     h2("Settings"),
                                     selectInput(inputId = "major_ui_datasetID", label = h4("Select a data set:"), choices = c("All",
                                                                                                                         unique(df_major$Species)))),
                              column(6,
                                     h2("Quantified proteins"),
                                     h2(""),
                                     dataTableOutput("major_MainDataTable")),
                              column(4,
                                     h2("Selected proteins:"),
                                     selectizeInput('major_ui_selectedIDs',
                                                    NULL,
                                                    choices=NULL,
                                                    multiple=TRUE),
                                     actionButton("update_major_MainDataTable",
                                                  "Update selection"))),
                     # textInput(inputId = "ui_proteinID", label = h4("Protein ID input:"), value = "Enter text..."))),
                     navbarPage("Proteome Analysis",
                                tabPanel("Exploratory Analysis",
                                         # h2("Sample correlation"),
                                         # p("To summarize the main characteristics of the leishmania infectome data set we perform an exploratory
                                         #   analysis, mainly focused on sample correlation. This way, we can assess for techincal sample 
                                         #   variance (replicates correlation) and experimental sample variance (groups, in this case,
                                         #   timepoints) in a graphical way. This analysis can be performed to the full data set or by 
                                         #   separating the data per species."),
                                         # selectInput(inputId = "major_datasetID", label = h4("Select a data set:"), choices = c("All",
                                         #                                                                                  unique(df_major$Species))),
                                         fluidRow(column(6,
                                                         h4("PCA plot"),
                                                         # p("Scatter plot depicting the output of a principal componetns analysis (PCA). Overall, 
                                                         #   replicates are expected to show relatively small distances between them, while sample 
                                                         #   groups often show higher distances between them, depending on the magnitude of 
                                                         #   intensity changes."),
                                                         plotOutput("major_PCA",height = "1000px")),
                                                  column(6,
                                                         h4("HeatMap plot"),
                                                         # p("Heat map depicting the pearson's correlation coefficient and the clustering between samples. Overall, 
                                                         #   replicates are expected to show high correlation between them, while sample 
                                                         #   groups might show lower correlations, depending on the magnitude of intensity changes."),
                                                         plotOutput("major_HeatMap",height = "1000px")))),
                                tabPanel("Clustering analysis",
                                         # h2("Self organizing maps (SOMs)"),
                                         # p("Self Organizing Maps (SOMs) are a tool for visualizing patterns in high dimensional data by producing a 
                                         #   2 dimensional representation, which displays meaningful patterns in the higher dimensional structure. It is
                                         #   an unsupervised machine learning technique that preserves the topological structure of the data.The size of
                                         #   the map grid is user defined and is the result of the multiplication of X and Y dimensions. This will result
                                         #   in n clusters in which data is expect to behave in a similar pattern. Default X and Y dimension are set to
                                         #   2 and 3 respectively, but should be adjusted per data set. The first adjustment criteria is trying to maximize
                                         #   the number of clusters with an homogenous compositions (individuals close to the mean). The second criteria is
                                         #   not having clusters with a very low count of individuals. This analysis can be performed to the full data set or by 
                                         #   separating the data per species."),
                                         fluidRow(
                                           # column(2,
                                           #        selectInput(inputId = "major_datasetID2", label = h4("Select a data set:"), choices = c("All",
                                           #                                                                        unique(df_major$Species)))),
                                           # column(2,
                                           #        textInput(inputId = "major_proteinID2", label = h4("Protein ID input:"), value = "Enter text...")),
                                           column(4,
                                                  bsCollapsePanel("Advanced settings",
                                                                  column(6,
                                                                         sliderInput("major_xDIM", "Select X dimensions:",
                                                                                     min = 1, max = 3, 
                                                                                     value = 2, step = 1)),
                                                                  column(6,
                                                                         sliderInput("major_yDIM", "Select Y dimensions:",
                                                                                     min = 1, max = 3,
                                                                                     value = 3, step = 1))))),
                                         fluidRow(column(8,
                                                         h4("Lineplot: Resulting clusters and their mean"),
                                                         # p("Line plot representation of the resulting n clusters. For each cluster's protein ID,
                                                         #   the standarized LFQ intensity (z-score(LFQ), y-axis) per timepoint (x-axis) is 
                                                         #   depicted in lightgrey. The cluster mean is depicted in each species color gradient.
                                                         #   Finally, and upon Protein ID selection, the selected target is highlighted in orange."),
                                                         plotOutput("major_ClusterLinePlot",height = "1000px")),
                                                  column(4,
                                                         h4("Barplot: Protein count per cluster"),
                                                         # p("Barplot showing the protein in count and species distribution per cluster"),
                                                         plotlyOutput("major_ClusterBarPlot",height = "1000px")))),
                                tabPanel("Pairwise Analysis - Time course",
                                         # h2("Time course comparison"),
                                         # p("The time course paiwise analysis results are depicted in this section. Thus, each timepoint is statistically
                                         #   compared to the one that precedes it.The analysis is performed to the full data set, but can
                                         #   be aslo visualized per species. The data distribution per timpoint is depicted in a boxplot shape.
                                         #   Upon protein ID input, its behaviour on the data distribution, as well as its statistical significance (line plot),
                                         #   is depicted. Finally, we can directly compare data between time course time points on a volcano plot."),
                                         # fluidRow(column(2,
                                         #                 selectInput(inputId = "major_speciesID", label = h4("Select data set:"), choices = c("All",
                                         #                                                                                                         unique(df_major$Species)))),
                                         #          column(2,
                                         #                 textInput(inputId = "major_proteinID", label = h4("Protein ID input:"), value = "Enter text...")),
                                         #          column(8)),
                                         fluidRow(column(6,
                                                         h4("Boxplot: Proteome distribution along the time course"),
                                                         # p("The protein IDs distribution is depicted as a boxplot. Thus, we have an indivdual boxplot per
                                                         #     each timepoint (x-axis) showing the proteins intensity, as log2(LFQ) (y-axis), distribution. Upon
                                                         #     protein ID input and when it is present on the data set, each ID replicate's intensity and mean line 
                                                         #     is depicted on the plot."),
                                                         plotOutput("major_Boxplot",height = "1000px")),
                                                  column(6,
                                                         h4("Dotplot: Selected ID statiscal significance"),
                                                         # p("Upon protein ID input and when present on the data set, a dotpot is generated. Thus, 
                                                         #   each replicate intensity, as log2(LFQ) (y-axis), per each timepoint (x-axis) is represented.
                                                         #   Additionally, a Welch's t-test is performed between adjacent timepoints whose results are plotted 
                                                         #   as adjusted (BH-FDR) p-values."),
                                                         plotOutput("major_Dotplot",height = "1000px")))),
                                tabPanel("Pairwise Analysis - Reference",
                                         # h2("Against reference comparison"),
                                         h4("Volcano plot: Individual IDs fold change and p-value"),
                                         # p("Volcano plot showing the LFQ intensity difference (log2(FoldChange), x-axis) and its significance (-log10(p-value),
                                         #   y-axis) for each protein ID. Thus each ID is represented by a light grey dot, when not significant, or a by a black
                                         #   dot, when significant. Default significance threshold is set to p-value < 0.05 and abs(log2(FoldChange)) > 2, but
                                         #   can be ajusted to user needs. Finally, and upon ID selection, the dot is higlighted in orange, when not significant, 
                                         #   or in red, when significant."),
                                         fluidRow(column(2,
                                                         selectInput(inputId = "major_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_major$TimePoint))),
                                                  column(2,
                                                         selectInput(inputId = "major_TP", label = h4("Select a time point:"), choices = unique(df_major$TimePoint))),
                                                  column(4,
                                                         bsCollapsePanel("Advanced settings",
                                                                         column(6,
                                                                                selectInput(inputId = "major_Pvalue", label = h4("Select a p-value:"), choices = c(0.05,0.01))),
                                                                         column(6,
                                                                                sliderInput("major_LFC", "Select a log2(FoldChange):",
                                                                                            min = 0, max = 5,
                                                                                            value = 1, step = 0.1))))),
                                         plotlyOutput("major_Volcanoplot",height = "1000px")))),
            # column(4,
            #        h4("Volcano plot: Individual IDs fold change and p-value"),
            #        p("blablabla volcanoplot"),
            #        fluidRow(column(6,
            #                        selectInput(inputId = "major_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_major$TimePoint))),
            #                 column(6,
            #                        selectInput(inputId = "major_TP", label = h4("Select a time point:"), choices = unique(df_major$TimePoint)))),
            #        plotlyOutput("Plot_major_8",height = "1000px")))))),
            tabPanel("L. mexicana",
                     fluidRow(column(2,
                                     h2("Settings"),
                                     selectInput(inputId = "mexicana_ui_datasetID", label = h4("Select a data set:"), choices = c("All",
                                                                                                                               unique(df_mexicana$Species)))),
                              column(6,
                                     h2("Quantified proteins"),
                                     h2(""),
                                     dataTableOutput("mexicana_MainDataTable")),
                              column(4,
                                     h2("Selected proteins:"),
                                     selectizeInput('mexicana_ui_selectedIDs',
                                                    NULL,
                                                    choices=NULL,
                                                    multiple=TRUE),
                                     actionButton("update_mexicana_MainDataTable",
                                                  "Update selection"))),
                     # textInput(inputId = "ui_proteinID", label = h4("Protein ID input:"), value = "Enter text..."))),
                     navbarPage("Proteome Analysis",
                                tabPanel("Exploratory Analysis",
                                         # h2("Sample correlation"),
                                         # p("To summarize the main characteristics of the leishmania infectome data set we perform an exploratory
                                         #   analysis, mainly focused on sample correlation. This way, we can assess for techincal sample 
                                         #   variance (replicates correlation) and experimental sample variance (groups, in this case,
                                         #   timepoints) in a graphical way. This analysis can be performed to the full data set or by 
                                         #   separating the data per species."),
                                         # selectInput(inputId = "mexicana_datasetID", label = h4("Select a data set:"), choices = c("All",
                                         #                                                                                  unique(df_mexicana$Species))),
                                         fluidRow(column(6,
                                                         h4("PCA plot"),
                                                         # p("Scatter plot depicting the output of a principal componetns analysis (PCA). Overall, 
                                                         #   replicates are expected to show relatively small distances between them, while sample 
                                                         #   groups often show higher distances between them, depending on the magnitude of 
                                                         #   intensity changes."),
                                                         plotOutput("mexicana_PCA",height = "1000px")),
                                                  column(6,
                                                         h4("HeatMap plot"),
                                                         # p("Heat map depicting the pearson's correlation coefficient and the clustering between samples. Overall, 
                                                         #   replicates are expected to show high correlation between them, while sample 
                                                         #   groups might show lower correlations, depending on the magnitude of intensity changes."),
                                                         plotOutput("mexicana_HeatMap",height = "1000px")))),
                                tabPanel("Clustering analysis",
                                         # h2("Self organizing maps (SOMs)"),
                                         # p("Self Organizing Maps (SOMs) are a tool for visualizing patterns in high dimensional data by producing a 
                                         #   2 dimensional representation, which displays meaningful patterns in the higher dimensional structure. It is
                                         #   an unsupervised machine learning technique that preserves the topological structure of the data.The size of
                                         #   the map grid is user defined and is the result of the multiplication of X and Y dimensions. This will result
                                         #   in n clusters in which data is expect to behave in a similar pattern. Default X and Y dimension are set to
                                         #   2 and 3 respectively, but should be adjusted per data set. The first adjustment criteria is trying to maximize
                                         #   the number of clusters with an homogenous compositions (individuals close to the mean). The second criteria is
                                         #   not having clusters with a very low count of individuals. This analysis can be performed to the full data set or by 
                                         #   separating the data per species."),
                                         fluidRow(
                                           # column(2,
                                           #        selectInput(inputId = "mexicana_datasetID2", label = h4("Select a data set:"), choices = c("All",
                                           #                                                                        unique(df_mexicana$Species)))),
                                           # column(2,
                                           #        textInput(inputId = "mexicana_proteinID2", label = h4("Protein ID input:"), value = "Enter text...")),
                                           column(4,
                                                  bsCollapsePanel("Advanced settings",
                                                                  column(6,
                                                                         sliderInput("mexicana_xDIM", "Select X dimensions:",
                                                                                     min = 1, max = 3, 
                                                                                     value = 2, step = 1)),
                                                                  column(6,
                                                                         sliderInput("mexicana_yDIM", "Select Y dimensions:",
                                                                                     min = 1, max = 3,
                                                                                     value = 3, step = 1))))),
                                         fluidRow(column(8,
                                                         h4("Lineplot: Resulting clusters and their mean"),
                                                         # p("Line plot representation of the resulting n clusters. For each cluster's protein ID,
                                                         #   the standarized LFQ intensity (z-score(LFQ), y-axis) per timepoint (x-axis) is 
                                                         #   depicted in lightgrey. The cluster mean is depicted in each species color gradient.
                                                         #   Finally, and upon Protein ID selection, the selected target is highlighted in orange."),
                                                         plotOutput("mexicana_ClusterLinePlot",height = "1000px")),
                                                  column(4,
                                                         h4("Barplot: Protein count per cluster"),
                                                         # p("Barplot showing the protein in count and species distribution per cluster"),
                                                         plotlyOutput("mexicana_ClusterBarPlot",height = "1000px")))),
                                tabPanel("Pairwise Analysis - Time course",
                                         # h2("Time course comparison"),
                                         # p("The time course paiwise analysis results are depicted in this section. Thus, each timepoint is statistically
                                         #   compared to the one that precedes it.The analysis is performed to the full data set, but can
                                         #   be aslo visualized per species. The data distribution per timpoint is depicted in a boxplot shape.
                                         #   Upon protein ID input, its behaviour on the data distribution, as well as its statistical significance (line plot),
                                         #   is depicted. Finally, we can directly compare data between time course time points on a volcano plot."),
                                         # fluidRow(column(2,
                                         #                 selectInput(inputId = "mexicana_speciesID", label = h4("Select data set:"), choices = c("All",
                                         #                                                                                                         unique(df_mexicana$Species)))),
                                         #          column(2,
                                         #                 textInput(inputId = "mexicana_proteinID", label = h4("Protein ID input:"), value = "Enter text...")),
                                         #          column(8)),
                                         fluidRow(column(6,
                                                         h4("Boxplot: Proteome distribution along the time course"),
                                                         # p("The protein IDs distribution is depicted as a boxplot. Thus, we have an indivdual boxplot per
                                                         #     each timepoint (x-axis) showing the proteins intensity, as log2(LFQ) (y-axis), distribution. Upon
                                                         #     protein ID input and when it is present on the data set, each ID replicate's intensity and mean line 
                                                         #     is depicted on the plot."),
                                                         plotOutput("mexicana_Boxplot",height = "1000px")),
                                                  column(6,
                                                         h4("Dotplot: Selected ID statiscal significance"),
                                                         # p("Upon protein ID input and when present on the data set, a dotpot is generated. Thus, 
                                                         #   each replicate intensity, as log2(LFQ) (y-axis), per each timepoint (x-axis) is represented.
                                                         #   Additionally, a Welch's t-test is performed between adjacent timepoints whose results are plotted 
                                                         #   as adjusted (BH-FDR) p-values."),
                                                         plotOutput("mexicana_Dotplot",height = "1000px")))),
                                tabPanel("Pairwise Analysis - Reference",
                                         # h2("Against reference comparison"),
                                         h4("Volcano plot: Individual IDs fold change and p-value"),
                                         # p("Volcano plot showing the LFQ intensity difference (log2(FoldChange), x-axis) and its significance (-log10(p-value),
                                         #   y-axis) for each protein ID. Thus each ID is represented by a light grey dot, when not significant, or a by a black
                                         #   dot, when significant. Default significance threshold is set to p-value < 0.05 and abs(log2(FoldChange)) > 2, but
                                         #   can be ajusted to user needs. Finally, and upon ID selection, the dot is higlighted in orange, when not significant, 
                                         #   or in red, when significant."),
                                         fluidRow(column(2,
                                                         selectInput(inputId = "mexicana_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_mexicana$TimePoint))),
                                                  column(2,
                                                         selectInput(inputId = "mexicana_TP", label = h4("Select a time point:"), choices = unique(df_mexicana$TimePoint))),
                                                  column(4,
                                                         bsCollapsePanel("Advanced settings",
                                                                         column(6,
                                                                                selectInput(inputId = "mexicana_Pvalue", label = h4("Select a p-value:"), choices = c(0.05,0.01))),
                                                                         column(6,
                                                                                sliderInput("mexicana_LFC", "Select a log2(FoldChange):",
                                                                                            min = 0, max = 5,
                                                                                            value = 1, step = 0.1))))),
                                         plotlyOutput("mexicana_Volcanoplot",height = "1000px")))),
            # column(4,
            #        h4("Volcano plot: Individual IDs fold change and p-value"),
            #        p("blablabla volcanoplot"),
            #        fluidRow(column(6,
            #                        selectInput(inputId = "mexicana_referenceTP", label = h4("Select a reference time point:"), choices = unique(df_mexicana$TimePoint))),
            #                 column(6,
            #                        selectInput(inputId = "mexicana_TP", label = h4("Select a time point:"), choices = unique(df_mexicana$TimePoint)))),
            #        plotlyOutput("Plot_mexicana_8",height = "1000px")))))),
            tabPanel("Orthology",
                     navbarPage("Orthology analysis",
                                tabPanel("Proteinortho analysis",
                                         fluidRow(
                                           column(12,
                                                  h2("Orthology table:"),
                                                  dataTableOutput("Ortho_table"))),
                                         fluidRow(
                                           column(4,
                                                  h2("Selected orthology groups:"),
                                                  selectizeInput('ui_orthoSelectedIDs',
                                                                 NULL,
                                                                 choices=NULL,
                                                                 multiple=TRUE),
                                                  actionButton("update_ortho_MainDataTable",
                                                               "Update selection"))),
                                         h4("Line plot: Expression mean at each time point"),
                                         plotOutput("OrthoPlot",height = "1000px"))))
                     
            
  )
             
  