######################################
############# Libraries ##############
######################################

library(ggplot2)
library(RColorBrewer)

######################################
############## Script ################
######################################

pg_quant <- read.csv("./QuantifiedProteins_mexicana.csv")

pattern <- "Majority.protein|Species|imputed.log2.LFQ.intensity"

df <- pg_quant[,grep(pattern = pattern, colnames(pg_quant))]

plot_df <- df[,1:3]

plot_df$TimePoint <- gsub("imputed.log2.LFQ.intensity.|_.*","",colnames(plot_df)[3])

plot_df$Replica <- gsub("imputed.log2.LFQ.intensity.*_","",colnames(plot_df)[3])

plot_df$Replica <- as.character(plot_df$Replica)

colnames(plot_df)[3] <- "Value"

for (i in 4:ncol(df)) {
  
  loop <- df[,c(1:2,i)]
  
  loop$TimePoint <- gsub("imputed.log2.LFQ.intensity.|_.*","",colnames(loop)[3])
  
  loop$Replica <- gsub("imputed.log2.LFQ.intensity.*_","",colnames(loop)[3])
  
  plot_df$Replica <- as.character(plot_df$Replica)
  
  colnames(loop)[3] <- "Value"
  
  plot_df <- rbind(plot_df, loop)
  
}

write.csv(plot_df, "./Database_data_mexicana.csv",row.names = F)

pattern <- "Majority.protein|Species|^pvalue_"

df <- pg_quant[,grep(pattern = pattern, colnames(pg_quant))]

pattern <- "Majority.protein|Species|.*mext005_mext000|.*mext020_mext005|.*mext060_mext020|.*mext120_mext060|.*mext240_mext120|.*mext480_mext240|.*mext720_mext480"

df <- df[,grep(pattern = pattern, colnames(df))]

test_df <- df[,1:3]

test_df$group1 <- gsub("pvalue_mex.*_","",colnames(test_df)[3])

test_df$group2 <- gsub("pvalue_|_mex.*","",colnames(test_df)[3])

colnames(test_df)[3] <- "p.adj"

for (i in 4:ncol(df)) {
  
  loop <- df[,c(1:2,i)]
  
  loop$group1 <- gsub("pvalue_mex.*_","",colnames(loop)[3])
  
  loop$group2 <- gsub("pvalue_|_mex.*","",colnames(loop)[3])
  
  colnames(loop)[3] <- "p.adj"
  
  test_df <- rbind(test_df, loop)
  
}

write.csv(test_df, "./Test_data_mexicana.csv",row.names = F)
