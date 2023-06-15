library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

dfmap <- fread (snakemake@input$map_s_p)

PDB_TEMPLATE <- as.character(unlist(dfmap$PDB))
names(PDB_TEMPLATE) <- as.character(unlist(dfmap$TEMPLATE))
chain <- dfmap$CHAIN[[1]]

df <- fread(snakemake@input$data)
df$ID <- 1:nrow(df)

df2 <- df %>%
    mutate(Mutation = strsplit(str_trim(Mutation,side = "both"),split=",|;"))
df2 <- df2 %>%
    unnest(Mutation) %>%
    rowwise() %>%
    mutate(
        POSITION_TEMPLATE =substr(x = Mutation,start = 2,stop = nchar(Mutation)-1),
        WT =substr(x = Mutation,start = 1,stop = 1),
        SUB =substr(x = Mutation,start = nchar(Mutation),stop = nchar(Mutation)),
        ) %>%
    mutate(POSITION_PDB = PDB_TEMPLATE[POSITION_TEMPLATE]) 
df3 <- df2 %>%
    filter(!is.na(POSITION_PDB)) %>%
    #where deletions are - sbstityue to alanine
    mutate(SUB = ifelse(SUB=='-',"A",SUB)) %>%
    mutate(MutationPDB = paste(WT,POSITION_PDB,SUB,sep=""))
df4 <- df3 %>%
    group_by(ID) %>%
    summarise(
        Mutation=paste(Mutation, collapse=","),
        MutationPDB=paste(MutationPDB, collapse=","),
        Template = first(Template),
        PDB = first(PDB)
        
    ) %>%
    ungroup()
dfhap <- df2 %>%
    arrange(ID,as.numeric(POSITION_TEMPLATE),WT,SUB) %>%
     group_by(ID) %>%
     summarise(Haplotype = paste(Mutation,collapse=","))
df4 <- df4 %>%
    left_join(dfhap, by=c("ID")) %>%
    mutate(CHAIN=chain) %>%
    group_by(Mutation,MutationPDB,Template,PDB,CHAIN) %>%
    summarise(Haplotype=paste(Haplotype, collapse="|")) %>%
    ungroup()
fwrite(x = df4,file = snakemake@output[[1]])

# general info
dfhap <- df2 %>%
    arrange(ID,as.numeric(POSITION_TEMPLATE),WT,SUB) %>%
     group_by(ID) %>%
     summarise(Haplotype = paste(Mutation,collapse=","))
df_withhap <- df %>%
    left_join(dfhap, by="ID")
fwrite(x = df_withhap, file = snakemake@output[[2]])
