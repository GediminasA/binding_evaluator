library(data.table)
library(dplyr)

terms = fread(snakemake@input$table)
ddg = fread(snakemake@input$ddg, col.names = c("ddG"))
fwrite(
    x = bind_cols(select(terms, label),ddg),
    file = snakemake@output$ddg
    )

