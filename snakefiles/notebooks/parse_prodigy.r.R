library(data.table)
library(dplyr)
library(stringr)



ln = readLines(snakemake@input$complex_prodigy)[[1]]
outdata = NA
if (grepl(pattern = "Predicted dissociation constant",x = ln)) {
    data = fread(cmd=paste("  sed 's/\\[/\\n/g' ",snakemake@input$complex_prodigy,
               " | grep -e 'No' -e 'Per' -e 'binding'  |
    sed 's/+//g' |
    sed 's/] //g' |
    sed 's/: /;/g' |
    sed 's/No./No/g' |
    sed 's/(kcal.mol-1)//g' |
    grep -v Reading "),sep=';',head=F,col.names = c("measure","value"))
    v1 = str_replace_all(string = data$measure,pattern = "-",replacement = "_")
    v2 = str_replace_all(string = v1,pattern = " ",replacement = "_")

    data$measure <- unlist(v2)
    outdata = t(data)
} else {
    head = "No_of_intermolecular_contacts	No_of_charged_charged_contacts	No_of_charged_polar_contacts	No_of_charged_apolar_contacts	No_of_polar_polar_contacts	No_of_apolar_polar_contacts	No_of_apolar_apolar_contacts	Percentage_of_apolar_NIS_residues	Percentage_of_charged_NIS_residues	Predicted_binding_affinity"
    head = strsplit(x=head,split="\t")[[1]]
    navals = rep(NA,length(head))
    datana <- data.table(measure=head,value=navals)
    outdata <- t(datana)
}


outdata = as.data.table(outdata)
colnames(outdata) <- unlist(outdata[1,],use.names = F)
outdata <- outdata[-c(1),]
outdata$COMPLEXID <- snakemake@params$id
fwrite(x = outdata,file = snakemake@output$complex_prodigy,row.names = F,col.names = T,sep="\t")
              
