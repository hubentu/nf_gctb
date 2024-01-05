#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

sfile <- args[1]
info_f <- args[2]

library(tidyverse)

## info
l1 <- read_table(info_f, n_max=1, col_names = FALSE, show_col_types = F)
if(l1[1,1]!="Chrom"){
    info <- read_table(info_f, col_names = F)
    colnames(info)[1:6] <- c("Chrom", "ID", "GenPos", "PhysPos", "A1", "A2")
}else{
    info <- read_table(info_f, col_names = T)
}
info <- info |> rowwise() |>
    mutate(var_sort = paste(c(Chrom, PhysPos, sort(c(A1, A2))), collapse="_")) |>
    select(ID, var_sort)

print(tibble(info))

fun <- DataFrameCallback$new(function(x, pos){
    d1 <- x |>
        mutate(
            freq=case_when(                
                !is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) & (effect_allele_BBJ==effect_allele_BCAC) ~ (Freq_effect_BBJ*73225+Freq_effect_BCAC*27172)/(73225+27172),
                !is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) & (effect_allele_BBJ!=effect_allele_BCAC) ~ (Freq_effect_BBJ*73225+(1-Freq_effect_BCAC)*27172)/(73225+27172),
                !is.na(effect_allele_BBJ) & is.na(effect_allele_BCAC) ~ Freq_effect_BBJ,
                is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) ~ Freq_effect_BCAC)
        ) |>
        mutate(
            n=case_when(                
                !is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) & (effect_allele_BBJ==effect_allele_BCAC) ~ 106722,
                !is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) & (effect_allele_BBJ!=effect_allele_BCAC) ~ 106722,
                !is.na(effect_allele_BBJ) & is.na(effect_allele_BCAC) ~ 79550,
                is.na(effect_allele_BBJ) & !is.na(effect_allele_BCAC) ~ 27172)
        ) |>
        select(unique_SNP_id, effect_allele_meta, non_effect_allele_meta, freq, BETA_meta, SE_meta, P_meta, n) |>
    separate(unique_SNP_id, c("chr", "pos", "ref", "alt"), sep="_") |>
    rowwise() |> mutate(var_sort = paste(c(chr, pos, sort(c(effect_allele_meta, non_effect_allele_meta))), collapse="_"))

    lookup <- c(SNP="ID", A1="effect_allele_meta", A2="non_effect_allele_meta", freq="freq", b="BETA_meta", se="SE_meta", p="P_meta", n="n")

    d1a <- d1 |> inner_join(info, by = "var_sort") |>
        select(ID, effect_allele_meta, non_effect_allele_meta, freq, BETA_meta, SE_meta, P_meta, n) |>
        dplyr::rename(all_of(lookup))
    ## mutate(n=106722)
    return(d1a)
})

dat <- read_delim_chunked(sfile, fun, delim=" ")
write_tsv(dat, paste0(sub(".txt", "", sfile), ".ma"))
