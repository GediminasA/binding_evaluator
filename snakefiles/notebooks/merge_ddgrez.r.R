library(dplyr)
library(data.table)


df <- lapply(X = snakemake@input, function(x){
    dfl <- fread(x) %>%
        mutate_at(vars(-("label")),as.double)
    return(dfl)
})
fwrite(x = bind_rows(df), file = snakemake@output[[1]])
