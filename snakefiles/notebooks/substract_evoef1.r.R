options(warn = -2)
library(data.table)
library(dplyr)

mut = fread(cmd = paste(" grep '=' ", snakemake@input$mut , " | ", "grep 'Total' " ), header = F, col.names = c("Term","A","MUT")) %>%
    select("MUT")
wt = fread(cmd = paste(" grep '=' ", snakemake@input$wt , " | ", "grep 'Total' " ), header = F, col.names = c("Term","A","WT")) %>%
    select("WT")

df = bind_cols(mut,wt) %>%
    mutate(DDG = MUT - WT, MUTATIONS = snakemake@wildcards$mutations)

fwrite(x = df, file =   snakemake@output$mut)
