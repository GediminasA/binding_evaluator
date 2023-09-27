library(dplyr)
library(data.table)

df <- lapply(X = snakemake@input, function(x){
    fread(x)
})
fwrite(x = bind_rows(df), file = snakemake@output[[1]])

