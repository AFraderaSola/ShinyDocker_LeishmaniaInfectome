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
######### Script ###########
############################

options(stringsAsFactors = FALSE)

source('myui.R', local = TRUE)
source('myserver.R')

shinyApp(
  ui = ui,
  server = server
)

