library(tidyverse)

collapse_fx <- function(fx_list) {
  exclude <- c("putative", "protein", "containing", "domain", "to", "for", "a", "the",
               "during", "by", "may", "in", "and", "of", "between", "predicted", "predicted:",
               "hypothetical", "uncharacterized", "related", "similar", "unknown", "family", 
               "unnamed", "product", "NA", NA, "a", "b", "c", "d", "e", "", "has", "with",
               "an", "activity", "probable", "sp.", "possible", as.character(rep(0:100, each=1))
  )
  result <- data.frame(words = gsub(",","",unlist(strsplit(fx_list,"[[:space:]]")))) %>%
    mutate(words = as.character(tolower(words))) %>%
    filter(!words %in% exclude) %>%
    filter(!startsWith(words, "[") & !startsWith(words, "(") & !startsWith(words, "\\") & !startsWith(words, "/")
           & !startsWith(words, ";") & !startsWith(words, "<") & !startsWith(words, ">") & !startsWith(words, ".")
           & !startsWith(words, "-") & !startsWith(words, "—") & !startsWith(words, "-")) %>%
    filter(!endsWith(words, "]") & !endsWith(words, ")") & !endsWith(words, "\\") & !endsWith(words, "/")
           & !endsWith(words, ";") & !endsWith(words, "<") & !endsWith(words, ">") & !endsWith(words, ".")
           & !endsWith(words, "-") & !endsWith(words, "—") & !endsWith(words, "-")) %>%
    drop_na() %>%
    group_by(words) %>%
    summarize(count = n()) %>%
    arrange(desc(count)) 
  
  if (nrow(result) > 10) {
    return (toString(result$words[1:10]))
  }
  else {
    return (toString(result$words))
  }
}

process_ogs <- function(cutoff) {
  args <- commandArgs(trailingOnly = TRUE)
  ogs = read.delim(args[1]) %>%
    filter(Gene.Name != "Gene Name") %>%
    mutate(total_species = n_distinct(Gene.Name)) %>%
    group_by(OG5) %>%
    summarize(n_species = n_distinct(Gene.Name),
              total_species = mean(total_species),
              perc_species = round((n_species/total_species)*100)) %>%
    filter(perc_species >= cutoff) %>%
    mutate(OG5 = as.character(OG5))
  
  return (ogs)
}

read_fx <- function() {
  og_fx <- read.delim("/Users/katzlab_admin/Desktop/Caitlin/aa_deflines_OrthoMCL-5.txt", sep = "|", header=FALSE,
                      col.names = c("x1", "NP_id", "OG5", "fx")) %>%
    separate(col="OG5", into=c("prefix", "OG5"), sep="_", remove=TRUE) %>%
    dplyr::select(c("OG5", "fx")) 
  og_fx$OG5 <- trimws(og_fx$OG5, which = "both")
  
  return (og_fx)
}

filter_ogs <- function(ogs) {
  og_fx <- read_fx()
  og_fx_fil <- og_fx %>%
    inner_join(ogs, by="OG5") %>%
    mutate(fx = trimws(as.character(fx)), which="both") %>%
    group_by(OG5) %>%
    summarize(n_species = mean(n_species),
              total_species = mean(total_species),
              perc_species = mean(perc_species),
              fx = collapse_fx(fx)) %>%
    arrange(desc(perc_species)) %>%
    mutate(OG5 = paste0("OG5_", OG5))
  
  return (og_fx_fil)
}

main <- function(cutoff) {
  args <- commandArgs(trailingOnly = TRUE)
  
  ogs <- process_ogs(cutoff)
  
  if (nrow(ogs) > 0) {
    og_fx_fil <- filter_ogs(ogs)
    
    for (i in 1:nrow(og_fx_fil)) {
      print(paste0(as.character(og_fx_fil$OG5[i]), " hit by ", as.character(og_fx_fil$perc_species[i]), "% of taxa"))
    }
    write.csv(og_fx_fil, file=paste0(unlist(str_split(args[1], pattern=".tsv"))[1], "_filtered.tsv"), row.names=FALSE)
    written = TRUE
    return()
  }
  else {
    while (cutoff > 0) {
      print(cutoff)
      cutoff = cutoff - 5
      main(cutoff)
    }
    return("No qualifying OGs.")
  }
}

main(cutoff=75)


