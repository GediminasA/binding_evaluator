suppressWarnings({
library(data.table)
library(dplyr)
    })

out = do.call("rbind",lapply(snakemake@input$evoef_data, function(x){
    fread(x,
          sep=",",
          colClasses = "character"
         )
}))
fwrite(x = out, file = snakemake@output$evoef_data)

suppressWarnings({

out = do.call("rbind",lapply(snakemake@input$cleaned_haps, function(x){fread(x,sep=",",  colClasses = "character")}))
fwrite(x = out, file = snakemake@output$cleaned_haps)
})
