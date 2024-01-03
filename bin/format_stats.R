#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

sfile <- args[1]
ldm <- args[2]

library(tidyverse)

## info
info_f <- list.files(ldm, ".*_chr.*.info$", full.names = TRUE)
info <- c()
for(i in info_f){
    l1 <- read_table(i, n_max=1, col_names = FALSE, show_col_types = F)
    if(l1[1,1]!="Chrom"){
        i1 <- read_table(i, col_names = F)
        colnames(i1)[1:6] <- c("Chrom", "ID", "GenPos", "PhysPos", "A1", "A2")
    }else{
        i1 <- read_table(i, col_names = T)
    }
    i1 <- i1 |> rowwise() |>
        mutate(var_sort = paste(c(Chrom, PhysPos, sort(c(A1, A2))), collapse="_")) |>
        select(ID, var_sort)
    info <- rbind(info, i1)
}

fun <- DataFrameCallback$new(function(x, pos){
    d1 <- x |>
        filter(effect_allele_iCOGs == effect_allele_Onco &
               effect_allele_Onco == effect_allele_meta) |>
        mutate(freq=(Freq_effect_iCOGs*37818+Freq_effect_Onco*58383)/(37818+58383)) |>
        select(var_name, effect_allele_meta, non_effect_allele_meta, freq, BETA_meta, SE_meta, P_meta) |>
    separate(var_name, c("chr", "pos", "ref", "alt"), sep="_") |>
    rowwise() |> mutate(var_sort = paste(c(chr, pos, sort(c(effect_allele_meta, non_effect_allele_meta))), collapse="_"))

    lookup <- c(SNP="ID", A1="effect_allele_meta", A2="non_effect_allele_meta", freq="freq", b="BETA_meta", se="SE_meta", p="P_meta")

    d1a <- d1 |> inner_join(info, by = "var_sort") |>
        select(ID, effect_allele_meta, non_effect_allele_meta, freq, BETA_meta, SE_meta, P_meta) |>
        dplyr::rename(all_of(lookup)) |>
        mutate(n=214675)
    return(d1a)
})

dat <- read_delim_chunked(sfile, fun, delim=" ")
write_tsv(dat, paste0(sub(".txt", "", sfile), ".ma"))



