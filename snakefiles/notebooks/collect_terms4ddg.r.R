suppressWarnings(library(dplyr))
suppressWarnings(library(data.table))
suppressWarnings(library(tidyr))

df = data.frame(label=snakemake@wildcards$stem)

cad_df <- fread(snakemake@input$cad) %>%
    select(CADscore = V5)

cad_dS <- fread(snakemake@input$cad) %>%
    mutate(dS = V7 - V6) %>%
    select(dS)

provdf <- fread(cmd = paste("grep  -A10000000000000000000 SCORE",snakemake@input$provean), sep = "\t") %>%
    select(CS = SCORE) %>%
    group_by() %>%
    summarise(CS = sum(CS))

SA_com <- fread(snakemake@input$dssp_complex) %>%
    rename(SA_com = V1)
SA_part <- fread(snakemake@input$dssp_part) %>%
    rename(SA_part = V1)
ddG_EvoEF <- fread(snakemake@input$evoef) %>%
    rename(ddG_EvoEF = V1)

open_terms = fread(snakemake@input$openmm) %>%
    pivot_wider(names_from = V1,values_from = V2) #%>%
    #select(PotentialEnergy, HarmonicBondForce, PeriodicTorsionForce, CustomTorsionForce, CMAPTorsionForce, LJForce, LennardJones, CMMotionRemover, HarmonicAngleForce, LennardJones14, CustomGBForce, CoulombForce)



out_df <- bind_cols(df, cad_df, cad_dS, SA_com, SA_part, ddG_EvoEF, open_terms, provdf)


fwrite(file = snakemake@output[[1]], x = out_df)
