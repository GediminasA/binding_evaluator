options(warn=-1, message = FALSE)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

dfmap <- fread (snakemake@input$map_s_p)

PDB_TEMPLATE <- as.character(unlist(dfmap$PDB))
names(PDB_TEMPLATE) <- as.character(unlist(dfmap$TEMPLATE)) #map template to PDB map
chain <- dfmap$CHAIN[[1]]

df <- fread(snakemake@input$data)
df$ID <- 1:nrow(df)

df2 <- df %>%
    mutate(Mutations = sub(pattern = " ",replacement = ",", x = Mutations,fixed = T)) %>%
    mutate(Mutations = strsplit(str_trim(Mutations,side = "both"),split=",|;"))
df2 <- df2 %>%
    unnest(Mutations) %>%
    rowwise() %>%
    mutate(
        POSITION_TEMPLATE =substr(x = Mutations,start = 2,stop = nchar(Mutations)-1),
        WT =substr(x = Mutations,start = 1,stop = 1),
        SUB =substr(x = Mutations,start = nchar(Mutations),stop = nchar(Mutations)),
        ) %>%
    mutate(POSITION_PDB = PDB_TEMPLATE[POSITION_TEMPLATE]) 

df3 <- df2 %>%
    filter(!is.na(POSITION_PDB)) %>%
    #where deletions are - sbstityue to alanine
    mutate(SUB = ifelse(SUB=='-',"A",SUB)) %>%
    mutate(MutationPDB = paste(WT,POSITION_PDB,SUB,sep="")) %>%
    mutate(MutationPDBwChain = paste(WT,chain,POSITION_PDB,SUB,sep=""))

df4 <- df3 %>%
    group_by(ID) %>%
    summarise(
        Mutations=paste(Mutations, collapse=","),
        MutationsID = first(MutationsID),
        MutationPDB=paste(MutationPDB, collapse=","),
        MutationPDBwChain=paste(MutationPDBwChain, collapse=","),
        Template = first(Template),
        PDB = first(PDB)   
    ) %>%
    ungroup()

dfhap <- df2 %>%
    arrange(ID,as.numeric(POSITION_TEMPLATE),WT,SUB) %>%
     group_by(ID) %>%
     summarise(Haplotype = paste(Mutations,collapse=","))

df4 <- df4 %>%
    left_join(dfhap, by=c("ID")) %>%
    mutate(CHAIN=chain) %>%
    group_by(Mutations,MutationPDB,MutationPDBwChain,Template,PDB,CHAIN) %>%
    summarise(
        Haplotype=paste(Haplotype, collapse="|"),
        MutationsID=paste(MutationsID, collapse="|"),
    ) %>%
    ungroup()
fwrite(x = df4,file = snakemake@output[[1]])

# general info
dfhap <- df2 %>%
    arrange(ID,as.numeric(POSITION_TEMPLATE),WT,SUB) %>%
     group_by(ID) %>%
     summarise(Haplotype = paste(Mutations,collapse=","))
df_withhap <- df %>%
    left_join(dfhap, by="ID")
fwrite(x = df_withhap, file = snakemake@output[[2]])
df_withhap
