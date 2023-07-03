shh <- suppressPackageStartupMessages # It's a library, so shhh!
options(warn=-2)
shh(library(data.table))
shh(library(dplyr))

rezmut = bind_rows(lapply(X = snakemake@input$evoEF1_rez, fread))

datmut = fread(snakemake@input$muation_data)

rezmut <- rezmut %>%
    left_join(datmut, by = c("MUTATIONS"="MutationPDB"))
rezmut

fwrite(x = rezmut, file = snakemake@output$muation_data)
