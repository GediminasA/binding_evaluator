library(dplyr)
library(data.table)

df <- fread(cmd=paste("grep -A 3 RESULTS",snakemake@input$tsv),skip = 1,select=c("V1","V3")) %>%
    mutate(V1=toupper(paste("FREESASA",V1,sep="_")))
df <- data.table(t(df))
colnames(df) <- unlist(df[1,],use.names = F)
df[1,] <- NA 
df <- df %>%
    filter(row_number() %in% c( 2))
df$COMPLEXID <- snakemake@params$complexid
df$FRAMEID <- snakemake@params$frameid
df$PART <- snakemake@params$part
fwrite(file = snakemake@output$tsv,sep = "\t",x = df)

