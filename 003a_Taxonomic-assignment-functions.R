# ------------------------------------------------------------------------------
# ------------------- Taxonomic assignment functions -----------------------------
# ------------------------------------------------------------------------------

# AUTHOR: Paula Pappalardo - pappalardop@si.edu

# CITATION: There are three sets of functions that were used

# Functions developed by Paula Pappalardo for Pappalardo et al., (2025):
# Pappalardo, P., Hemmi, J. M., Machida, R. J., Leray, M., Collins, A. G., & Osborn, K. J. (2025). Taxon-specific BLAST percent identity thresholds for identification of unknown sequences using metabarcoding. Methods in Ecology and Evolution, 16, 2380–2394. https://doi.org/10.1111/2041-210X.70147 (Taxon-specific thresholds MS)
# Please cite the original publication if you reuse them

# Functions developed by Pappalardo & Palmer (2025) to compare BLAST hits 
# from different reference databases
# Pappalardo P. & E. Palmer. (2025). paulapappalardo/compare-reference-databases: First public release of comparing reference database functions (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.17478719

# Functions written by Paula Pappalardo for taxonomic standardization

# ---------Pappalardo et al (2025)
filterBlastResults <- function(results, n_to_keep, min_cov, min_length, max_evalue){
  # Keep only X filtered rows from blast results
  #
  # Args:
  # results: dataframe with blast results
  # n_to_keep: how many rows to keep after filtering - to speed up next steps
  # min_cov: minimum percent coverage 
  # min_length: minimum alignment length
  # max_evalue: maximum evalue to keep
  # How it works:
  # first, the function removes bad quality hits and blast to itself
  # then it slices part of the data to speed up computation later
  # TODO:  user needs to define min_cov, min_length, and max_evalue
  results_filtered <- results %>%
    filter(!(query_seqid == result_seqid)) %>%
    filter(percent_coverage >= min_cov, length >= min_length, evalue <= max_evalue) %>%
    group_by(query_seqid) %>% 
    slice_head(., n = n_to_keep)
  # return filtered dataset
  return(results_filtered)
}

tagIdenticalHits <- function(minidf){
  # Aim: identify which rows have identical similarity measures from BLAST result
  # Args
  #   minidf: subset of quality filtered BLAST results corresponding to one query_seq 
  # Notes: function meant to be used with group_by and group_map, grouping by query seq_id
  # 
  # Output: minidf with three extra columns:
  #   sim_values - string combining all the similarity measures
  #   tagged = TRUE in any of second or more rows will indicate if identical to first row (will be used by common ancestor function)
  #   best_check = "to check" for rows that had identical measures, and "is best" for the rest (will be used to split who needs common ancestor function)
  # set up columns we need
  working_df <- minidf %>% 
    mutate(sim_values = paste(evalue, bitscore, length, percent_coverage, nident, percent_identity, sep = "_"),
           tagged = NA) 
  # STEP 1: Loop by row to identify the rows that have similarity values identical to the first row
  # Note that this step will also tag T the hits with only one best hit
  for(i in 1:nrow(working_df)){
    #print(paste(working_df$query_seqid[1], "row", i, sep = " ")) # for debugging
    if (i == 1){
      answer <- T 
      working_df$tagged[1] <- answer 
    }
    else{
      answer <- working_df$sim_values[i] == working_df$sim_values[1]
      working_df$tagged[i] <- answer
    }
  }
  # STEP 2: add a column to separate the ones that had multiple hits with identical similarity measures
  tagged_df <- working_df %>% 
    mutate(best_check = ifelse(sum(tagged, na.rm = T) > 1, "to check", "is best")) 
  # remove things not in use
  rm(i, answer)
  # return final dataframe
  return(tagged_df)
}

# This function was edited from the getCommonAncestor_benchmarking
# in Pappalardo et al 2025
getCommonAncestor <- function(df_withtax){#last test on 2024-08-30
  # Aim: condense BLAST results with identical similarity measures
  # Args
  #   minidf: subset tagged BLAST results with higher taxonomy corresponding to one query_seq 
  # (NEEDS columns result_seqid, kingdom, phylum, class, order, family, genus, species, length, percent_coverage and percent_identity)
  # Dependencies: needs functions assign_unique, assign_unique_conservative, and assign_unique_ancestor_double
  # Output: reduced dataframe with only the higher taxonomy trusted for each unkwnon sequence
  #
  # First step: create a summarized dataset comparing taxonomic levels using custom made ancestor functions
  df_nodup <- df_withtax %>% 
    filter(tagged == TRUE) %>% # only combine taxonomy for the ones with identical matches
    group_by(query_seqid) %>% 
    summarize(result_seqid = assign_unique(result_seqid), # keep all match names with identical measure separated by dashes
              kingdom = assign_unique_conservative(kingdom), # if kingdom is not consistent, assign NA
              phylum = assign_unique_ancestor_double(vec_lower = phylum, vec_higher = kingdom),
              class = assign_unique_ancestor_double(vec_lower = class, vec_higher = phylum),
              order = assign_unique_ancestor_double(vec_lower = order, vec_higher = class),
              family = assign_unique_ancestor_double(vec_lower = family, vec_higher = order),
              genus = assign_unique_ancestor_double(vec_lower = genus, vec_higher = family),
              species = assign_unique_ancestor_double(vec_lower = species, vec_higher = genus),
              length = length[1],
              bitscore = bitscore[1],
              percent_coverage = percent_coverage[1],
              percent_identity = percent_identity[1],
              taxonomy_source = assign_unique(taxonomy_source)) 
  # Second step: if kingdom is NA, make everything NA (in case some very weird case is out there...unlikely)
  df_nodup_checked <- df_nodup %>% 
    mutate(phylum = ifelse(is.na(kingdom), NA_character_, phylum),
           class = ifelse(is.na(kingdom), NA_character_, class),
           order = ifelse(is.na(kingdom), NA_character_, order),
           family = ifelse(is.na(kingdom), NA_character_, family),
           genus = ifelse(is.na(kingdom), NA_character_, genus),
           species = ifelse(is.na(kingdom), NA_character_, species)
    ) %>% 
    # make kingdom Unknown if it was NA so that it count in tables and works for making sciname
    mutate(kingdom = ifelse(is.na(kingdom), "Unknown", kingdom))
  return(df_nodup_checked)
}

assign_unique <- function(vec){
  # Collapses unique values ignoring NAs
  # Args:
  # vec: a vector
  #TODO: Decide which outcome you want as "answer" in the else statement
  # if you want to see all the unique components, leave as is, separator is a dash
  # if you want to see NA when there is not agreement, comment out the second "answer" line within the else statement
  
  if(length(na.omit(vec)) == 0){ # check for a vector of ONLY NAs
    answer <- NA
  } else {
    answer <- paste(unique(na.omit(vec)), collapse = "-")
    #answer <- ifelse(grepl("-", answer), NA, answer)
  }
  return(answer)
}

assign_unique_conservative <- function(vec){
  # Function to assign taxonomic level IF similarity values are identical
  # Conservative, assumes that a NA compared with a true taxa level is NOT an agreement!
  # meant for higher taxonomic level (e.g., kingdom) where not having a name should not be allowed
  # last test on 2024-07-11
  # We have three scenarios: all NA, all tax levels (matching or not), mix of tax levels and NAs
  if(length(na.omit(vec)) == 0){ # check for a vector of ONLY NAs
    answer <- NA
  } else {if(sum(is.na(vec)) == 0){# check for vector of NO nas
    answer <- paste(unique(vec), collapse = "-") # different things separated by a dash
    answer <- ifelse(grepl("-", answer), NA, answer) # if there were different things make it NA
  } else{# there were some NAs implying conflict at high taxonomic level
    answer <- NA
  }
  }
  return(answer)
}

assign_unique_ancestor_double <- function(vec_lower, vec_higher){
  # Function to assign taxonomic level IF similarity values are identical
  # an NA in a species belonging to the same genus means keep species
  # an NA in genus of same family as the unique genus means keep genus
  # re tested on 2024-07-11
  # Args:
  # vec_lower: lower taxonomic rank we are checking (e.g., species)
  # vec_higher: higher taxonomic rank we are checking (e.g., genus if species was the lower one)
  # 
  # We have three scenarios: all NA, all tax levels (matching or not), mix of tax levels and NAs
  if(length(na.omit(vec_lower)) == 0){ # check for a vector of ONLY NAs
    answer <- NA
  } else {if(sum(is.na(vec_lower)) == 0){# check for vector of NO nas
    answer <- paste(unique(vec_lower), collapse = "-")
    answer <- ifelse(grepl("-", answer), NA, answer)
  } else{# there were some NAs - now we need to look up higher taxonomic level to make a decision
    # higher taxonomic level needs to be unique and not NAs
    if(length(na.omit(unique(vec_higher))) == 1 & length(na.omit(vec_higher)) != 0){
      answer <- paste(na.omit(unique(vec_lower)), collapse = "-")
      answer <- ifelse(grepl("-", answer), NA, answer)
    } else {answer <- NA}
  }
  }
  return(answer)
}

# ----------Pappalardo and Palmer (2025)

addScinameLevel <- function(mydf){
  # Finds the taxonomic resolution for a scientific name
  # Args:
  #  mydf: dataframe with tax columns for phylum, class, order, family, genus, and species 
  #  TODO: consider activating subfamily if you are dealing with BOLD output
  #  but currently BOLD taxonomy prepared by Paula will have this already
  # Ouput:
  #   same dataframe with additional column named sciname_level
  # This function is from XXX repo
  mydf_ed <- mydf %>% 
    mutate(sciname_level = case_when(sciname %in% na.omit(species) ~ "species",
                                     sciname %in% na.omit(genus) ~ "genus",
                                     sciname %in% na.omit(family) ~ "family",
                                     #sciname %in% na.omit(subfamily) ~ "subfamily", # if using
                                     sciname %in% na.omit(order) ~ "order",
                                     sciname %in% na.omit(class) ~ "class",
                                     sciname %in% na.omit(phylum) ~ "phylum",
                                     sciname %in% na.omit(kingdom) ~ "kingdom")) %>% 
    relocate(sciname_level, .after = sciname)
  # return modified dataframe
  return(mydf_ed)
}

labelFinalMatch <- function(mydf, label){
  # Condenses taxonomy to one tax string and adds reference database label
  # Args:
  #   mydf: blast final dataset (after quality filters, tags, and identical matches were solved)
  #   label: name to represent the origin database (e.g., "mlml" or "midori")
  # Output: the function will remove higher taxonomy columns and add a taxonomy string column
  # Usage: for dplyr use - labelFinalMatch(., label = "mlml")
  # 
  # combine all taxonomic levels into one string
  mydf_ed <- mydf %>% 
    mutate(tax_string = paste(kingdom, phylum, class, order, family, genus, species, sep = ";")) %>% 
    select(-kingdom, -phylum, -class, -order, -family, -genus, -species)
  # rename columns to indicate database source
  newnames <- paste(names(mydf_ed), label, sep = "_")
  names(mydf_ed) <- newnames
  return(mydf_ed)
}

compareBOLDandMIDORI <- function(mydf_bold, mydf_midori, identity_th){
  # Compare matches to bold and midori for same sequence
  # Args
  # mydf_bold: dataframe with blast results from BOLD
  # mydf_midori: dataframe with blast results from MIDORI
  # identity_th: which identity threshold you want to assume "good matches"
  # both dataframes need columns for similarity metrics and taxonomic levels 
  # prepare data to be compared
  # midori
  mydf_midori_ed <- mydf_midori %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "midori") %>% 
    dplyr::rename(query_seqid = query_seqid_midori)
  # bold
  mydf_bold_ed <- mydf_bold %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "bold") %>% 
    dplyr::rename(query_seqid = query_seqid_bold)
  
  # merge and compare matches
  mydf_ed <- mydf_midori_ed %>% 
    full_join(mydf_bold_ed, by = "query_seqid") %>%
    # first check for agreement on sciname and sciname level
    mutate(sciname_agree = ifelse(sciname_bold == sciname_midori, TRUE, FALSE),
           sciname_level_agree = ifelse(sciname_level_bold == sciname_level_midori, TRUE, FALSE)) %>% 
    # then identify the best identity for later
    mutate(best_db_identity  = case_when(is.na(percent_identity_bold) & !is.na(percent_identity_midori) ~ "keep midori - bold NA",
                                         is.na(percent_identity_midori) & !is.na(percent_identity_bold) ~ "keep bold - midori NA",
                                         percent_identity_midori > percent_identity_bold ~ "keep midori",
                                         percent_identity_midori <= percent_identity_bold ~ "keep bold",
                                         percent_identity_midori == percent_identity_bold ~ "identical",
                                         T ~ "CHECK")) %>% 
    # now try to come up with the outcomes based on the info we have:
    mutate(match_outcome = case_when(sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_midori >= identity_th & percent_identity_bold >= identity_th ~ "tricky and relevant - separate function",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_midori <= identity_th & percent_identity_bold <= identity_th ~ "tricky and irrelevant - separate function",
                                     percent_identity_bold > identity_th & percent_identity_midori <= identity_th ~ "keep bold - bold higher",
                                     percent_identity_midori > identity_th & percent_identity_bold <= identity_th ~ "keep midori - midori higher",
                                     is.na(percent_identity_bold) & !is.na(percent_identity_midori) ~ "keep midori - bold NA",
                                     is.na(percent_identity_midori) & !is.na(percent_identity_bold)  ~ "keep bold - midori NA",
                                     is.na(tax_string_midori) & is.na(tax_string_bold) ~ "no match from either",
                                     sciname_agree == TRUE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == TRUE & sciname_level_agree == FALSE ~ best_db_identity,# few cases, e.g. Oomycota
                                     T ~ "CHECK"))
  # create ordered factor of the taxonomic level
  higher_tax <- factor(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                       levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = ifelse(which(higher_tax == sciname_level_bold) > which(higher_tax == sciname_level_midori), "keep bold", "keep midori")) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2)
  
  return(mydf_ed_final)
  
}

pickFinalTax_BoldVsMidori <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: classified dataframe after running compareDatabases(), compareTaxMatches(), and classifyAssignments()
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and will add column database_source to track where the match came from
  # Usage: pickFinalTax(.) 
  #
  mydf_ed <- mydf %>% 
    dplyr::mutate(tax_string = case_when(grepl("keep midori", match_outcome) ~ tax_string_midori,
                                         grepl("keep bold", match_outcome) ~ tax_string_bold,
                                         T ~ "CHECK"),
                  sciname = case_when(grepl("keep midori", match_outcome) ~ sciname_midori,
                                      grepl("keep bold", match_outcome) ~ sciname_bold,
                                      T ~ "CHECK"),
                  sciname_level = case_when(grepl("keep midori", match_outcome) ~ sciname_level_midori,
                                            grepl("keep bold", match_outcome) ~ sciname_level_bold,
                                            T ~ "CHECK"),
                  taxonomy_source = case_when(grepl("keep midori", match_outcome) ~ taxonomy_source_midori,
                                              grepl("keep bold", match_outcome) ~ taxonomy_source_bold,
                                              T ~ "CHECK"),
                  database_source = case_when(grepl("keep midori", match_outcome) ~ "MIDORI2",
                                              grepl("keep bold", match_outcome)~ "BOLD",
                                              T ~ "CHECK"),
                  match_name = case_when(grepl("keep midori", match_outcome) ~ result_seqid_midori,
                                         grepl("keep bold", match_outcome)~ result_seqid_bold,
                                         T ~ "CHECK"),
                  bitscore = case_when(grepl("keep midori", match_outcome) ~ bitscore_midori,
                                       grepl("keep bold", match_outcome) ~ bitscore_bold,
                                       T ~ NA),
                  percent_coverage = case_when(grepl("keep midori", match_outcome) ~ percent_coverage_midori,
                                               grepl("keep bold", match_outcome) ~ percent_coverage_bold,
                                               T ~ NA),
                  percent_identity = case_when(grepl("keep midori", match_outcome) ~ percent_identity_midori,
                                               grepl("keep bold", match_outcome) ~ percent_identity_bold,
                                               T ~ NA)) %>%
    dplyr::select(query_seqid, database_source, match_name, bitscore, percent_coverage, percent_identity, sciname, sciname_level, tax_string, taxonomy_source) 
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
}

compareMLMLandLavrador<- function(mydf_mlml, mydf_lavrador, identity_th){
  # Compare matches between MLML and Lavrador for each unknown sequence
  # Args
  # mydf_mlml: dataframe with blast results from BOLD
  # mydf_labrador: dataframe with blast results from MIDORI
  # identity_th: which identity threshold you want to assume "good matches"
  # both dataframes need columns for similarity metrics and taxonomic levels 
  # prepare data to be compared
  # midori
  mydf_mlml_ed <- mydf_mlml %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "mlml") %>% 
    dplyr::rename(query_seqid = query_seqid_mlml)
  # bold
  mydf_lavrador_ed <- mydf_lavrador %>% 
    addScinameLevel() %>% 
    labelFinalMatch(., label = "lavrador") %>% 
    dplyr::rename(query_seqid = query_seqid_lavrador)
  
  # merge and compare matches
  mydf_ed <- mydf_mlml_ed %>% 
    full_join(mydf_lavrador_ed, by = "query_seqid") %>%
    # first check for agreement on sciname and sciname level
    mutate(sciname_agree = ifelse(sciname_mlml == sciname_lavrador, TRUE, FALSE),
           sciname_level_agree = ifelse(sciname_level_mlml == sciname_level_lavrador, TRUE, FALSE)) %>% 
    # then identify the best identity for later
    mutate(best_db_identity  = case_when(is.na(percent_identity_mlml) & !is.na(percent_identity_lavrador) ~ "keep lavrador - mlml NA",
                                         is.na(percent_identity_lavrador) & !is.na(percent_identity_mlml) ~ "keep mlml - lavrador NA",
                                         percent_identity_lavrador > percent_identity_mlml ~ "keep lavrador",
                                         percent_identity_lavrador <= percent_identity_mlml ~ "keep mlml",
                                         percent_identity_lavrador == percent_identity_mlml ~ "identical",
                                         T ~ "CHECK")) %>% 
    # now try to come up with the outcomes based on the info we have:
    mutate(match_outcome = case_when(sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_lavrador >= identity_th & percent_identity_mlml >= identity_th ~ "tricky and relevant - separate function",
                                     sciname_agree == FALSE & sciname_level_agree == FALSE &
                                       percent_identity_lavrador <= identity_th & percent_identity_mlml <= identity_th ~ "tricky and irrelevant - separate function",
                                     percent_identity_mlml > identity_th & percent_identity_lavrador <= identity_th ~ "keep mlml - mlml higher",
                                     percent_identity_lavrador > identity_th & percent_identity_mlml <= identity_th ~ "keep lavrador - lavrador higher",
                                     is.na(percent_identity_mlml) & !is.na(percent_identity_lavrador) ~ "keep lavrador - mlml NA",
                                     is.na(percent_identity_lavrador) & !is.na(percent_identity_mlml)  ~ "keep mlml - lavrador NA",
                                     is.na(tax_string_lavrador) & is.na(tax_string_mlml) ~ "no match from either",
                                     sciname_agree == TRUE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == FALSE & sciname_level_agree == TRUE ~ best_db_identity,
                                     sciname_agree == TRUE & sciname_level_agree == FALSE ~ best_db_identity,# few cases, e.g. Oomycota
                                     T ~ "CHECK"))
  # create ordered factor of the taxonomic level
  higher_tax <- factor(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                       levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), ordered = T)
  
  # separate "Tricky" cases where we want to compare taxonomic resolution achieved by each database
  mydf_ed_1 <- mydf_ed %>% 
    filter(grepl("tricky", match_outcome)) %>% 
    rowwise() %>% 
    mutate(match_outcome = ifelse(which(higher_tax == sciname_level_mlml) > which(higher_tax == sciname_level_lavrador), "keep mlml", "keep lavrador")) %>% 
    ungroup()
  
  mydf_ed_2 <- mydf_ed %>% 
    filter(!grepl("tricky", match_outcome))
  
  mydf_ed_final <- bind_rows(mydf_ed_1, mydf_ed_2)
  
  return(mydf_ed_final)
  
}

pickFinalTax_MLMLvsLavrador <- function(mydf){
  # Pick best final id and keep final percent identity
  # Args
  #   mydf: classified dataframe after running compareDatabases(), compareTaxMatches(), and classifyAssignments()
  # Output
  #   The function will keep only the higher taxonomy of the best hit, and will add column database_source to track where the match came from
  # Usage: pickFinalTax(.) 
  #
  mydf_ed <- mydf %>% 
    mutate(tax_string = case_when(grepl("keep mlml", match_outcome) ~ tax_string_mlml,
                                  grepl("keep lavrador", match_outcome) ~ tax_string_lavrador,
                                  T ~ "CHECK"),
           sciname = case_when(grepl("keep mlml", match_outcome) ~ sciname_mlml,
                               grepl("keep lavrador", match_outcome) ~ sciname_lavrador,
                               T ~ "CHECK"),
           sciname_level = case_when(grepl("keep mlml", match_outcome) ~ sciname_level_mlml,
                                     grepl("keep lavrador", match_outcome) ~ sciname_level_lavrador,
                                     T ~ "CHECK"),
           taxonomy_source = case_when(grepl("keep mlml", match_outcome) ~ taxonomy_source_mlml,
                                       grepl("keep lavrador", match_outcome) ~ taxonomy_source_lavrador,
                                       T ~ "CHECK"),
           database_source = case_when(grepl("keep mlml", match_outcome) ~ "mlml",
                                       grepl("keep lavrador", match_outcome)~ "lavrador",
                                       T ~ "CHECK"),
           match_name = case_when(grepl("keep mlml", match_outcome) ~ result_seqid_mlml,
                                  grepl("keep lavrador", match_outcome)~ result_seqid_lavrador,
                                  T ~ "CHECK"),
           bitscore = case_when(grepl("keep mlml", match_outcome) ~ bitscore_mlml,
                                grepl("keep lavrador", match_outcome) ~ bitscore_lavrador,
                                T ~ NA),
           percent_coverage = case_when(grepl("keep mlml", match_outcome) ~ percent_coverage_mlml,
                                        grepl("keep lavrador", match_outcome) ~ percent_coverage_lavrador,
                                        T ~ NA),
           percent_identity = case_when(grepl("keep mlml", match_outcome) ~ percent_identity_mlml,
                                        grepl("keep lavrador", match_outcome) ~ percent_identity_lavrador,
                                        T ~ NA)) %>%
    select(query_seqid, database_source, match_name, bitscore, percent_coverage, percent_identity, sciname, sciname_level, tax_string, taxonomy_source) 
  # return modified dataframe that includes only the final chosen hit from each database
  return(mydf_ed)
}

#--------------Functions for taxonomy standardization
# developed by Paula Pappalardo for this and other projects

matchToGBIF <- function(species_list){
  library(rgbif)
  # match species list to GBIF taxonomic backbone
  # Args
  # species_list = vector with list of species to search
  # Create GBIF names dictionary
  list_matched <- map_dfr(species_list, name_backbone) %>%
    select(verbatim_name, canonicalName,
           species, rank,
           status, confidence,
           matchType, kingdom,
           phylum, class,
           order, family,
           genus, synonym)
  return(list_matched)
}


# GBIF taxonomy using name_backbone function - it sometimes fails to return a match
# Worms taxonomy checked using classification(taxa_name, db = "worms") from taxize
# Have in mind that same names are ambiguous, the correct classification for our 
# purposes was checked in the full NCBI blast results to find the correct match 
# (e.g., the correct phylum for ambiguous taxa)
gbif_updates <- tribble(
  ~verbatim_name, ~kingdom, ~phylum, ~class, ~order, ~family, ~genus, ~species, ~taxonomy_source,
  # the "manual fix - gbif" usually are cases where the higher taxonomy is available for a species in the same genus
  # As
  "Alaria crassifolia", "Chromista", "Ochrophyta", "Phaeophyceae", "Laminariales", "Alariaceae", "Alaria", "Alaria crassifolia","manual fix - gbif",
  "Alcyonidium",  "Animalia", "Bryozoa", "Gymnolaemata", "Ctenostomatida", "Alcyonidiidae", "Alcyonidium", NA_character_,"manual fix - worms",
  "Aegina", "Animalia", "Cnidaria", "Hydrozoa", "Narcomedusae", "Aeginidae", "Aegina", NA_character_, "manual fix - worms",
  "Achlya", "Chromista", "Oomycota", "Peronosporea", "Saprolegniales", "Saprolegniaceae", "Achlya", NA_character_,"manual fix - gbif",
  "Alcyonidium sp.", "Animalia", "Bryozoa", "Gymnolaemata", "Ctenostomatida", "Alcyonidiidae", "Alcyonidium", NA_character_,"manual fix - worms",
  "Amathia", "Animalia", "Bryozoa", "Gymnolaemata", "Ctenostomatida", "Vesiculariidae", "Amathia", NA_character_,"manual fix - worms",
  "Amathia sp. bowerbankia", "Animalia", "Bryozoa", "Gymnolaemata", "Ctenostomatida", "Vesiculariidae", "Amathia", NA_character_,"manual fix - worms",
  "Amathia sp. Bowerbankia", "Animalia", "Bryozoa", "Gymnolaemata", "Ctenostomatida", "Vesiculariidae", "Amathia", NA_character_,"manual fix - worms",
  "Ampedus cordatus", "Animalia","Arthropoda", "Insecta", "Coleoptera", "Elateridae", "Ampedus", "Ampedus cordatus", "manual fix - gbif",
  "Aprion virescens", "Animalia", "Chordata", "Teleostei", "Eupercaria incertae sedis", "Lutjanidae", "Aprion", "Aprion virescens","manual fix - worms",
  "Aurelia", "Animalia", "Cnidaria", "Scyphozoa", "Semaeostomeae", "Ulmaridae", "Aurelia", NA_character_, "manual fix - gbif",
  # Bs
  "Balanus", "Animalia","Arthropoda", "Maxillopoda", "Sessilia", "Balanidae", "Balanus", NA_character_, "manual fix - gbif",
  "Baicaloclepsis grubei", "Animalia", "Annelida", "Clitellata", "Rhynchobdellida", "Glossiphoniidae", "Baicaloclepsis", "Baicaloclepsis grubei", "manual fix -ncbi",
  "Barentsia sp.", "Animalia", "Entoprocta", NA_character_, NA_character_, "Barentsiidae", "Barentsia", NA_character_, "manual fix - worms",
  "Beroe cucumis", "Animalia", "Ctenophora", "Nuda", "Beroida", "Beroidae", "Beroe", "Beroe cucumis", "manual fix - worms",
  # Cs
  "Calanus", "Animalia", "Arthropoda", "Copepoda", "Calanoida", "Calanidae", "Calanus", NA_character_, "manual fix - gbif",
  "Cafileria marina", "Chromista", NA_character_, "Bigyra", "Bicosoecida", NA_character_, "Cafileria", "Cafileria marina", "manual fix = gbif/ncbi",
  "Calephelis virginiensis", "Animalia", "Arthropoda", "Insecta", "Lepidoptera", "Riodinidae", "Calephelis", "Calephelis virginiensis", "manual fix",
  "Caprella", "Animalia", "Arthropoda", "Malacostraca", "Amphipoda", "Caprellidae", "Caprella", NA_character_, "manual fix - gbif",
  "Caprella sp.", "Animalia", "Arthropoda", "Malacostraca", "Amphipoda", "Caprellidae", "Caprella", NA_character_, "manual fix - gbif",
  "Cecidomyiidae",  "Animalia", "Arthropoda", "Insecta", "Diptera", "Cecidomyiidae", NA_character_, NA_character_, "manual fix - gbif",
  "Cephalothrix", "Animalia", "Nemertea", "Paleonemertea", "Archinemertea", "Cephalotrichidae", "Cephalothrix", NA_character_, "manual fix - gbif",
  "Chaetognatha", "Animalia", "Chaetognatha", NA_character_, NA_character_,  NA_character_, NA_character_, NA_character_, "manual fix - worms",
  "Chlorella desiccata", "Plantae", "Chlorophyta", "Trebouxiophyceae", "Chlorellales", "Chlorellaceae", "Chlorella", "Chlorella desiccata","manual fix - gbif",
  "Chlorella", "Plantae", "Chlorophyta", "Trebouxiophyceae", "Chlorellales", "Chlorellaceae", "Chlorella",  NA_character_,"manual fix - gbif",
  "Chone", "Animalia", "Annelida", "Polychaeta", "Sabellida", "Sabellidae", "Chone", NA_character_, "manual fix - worms",
  "Clione", "Animalia", "Mollusca", "Gastropoda", "Pteropoda", "Clionidae", "Clione", NA_character_, "manual fix - gbif",
  "Chlorophyta", "Plantae", NA_character_, NA_character_, NA_character_,  NA_character_, NA_character_, NA_character_, "manual fix - worms",
  "Clytia",  "Animalia", "Cnidaria", "Hydrozoa", "Leptothecata", "Campanulariidae", "Clytia", NA_character_, "manual fix - gbif",
  "Cochliopodium crystalli", "Protozoa", "Amoebozoa", "Discosea","Himatismenida", "Cochliopodiidae", "Cochliopodium", "Cochliopodium crystalli", "manual fix - gbif",
  "Cochliopodium marrii", "Protozoa", "Amoebozoa", "Discosea","Himatismenida", "Cochliopodiidae", "Cochliopodium", "Cochliopodium marrii", "manual fix - gbif",
  "Coccomyxa", "Plantae", "Chlorophyta", NA_character_, "Trebouxiophyceae", NA_character_, "Coccomyxa", NA_character_, "manual fix - gbif",
  "Coralline crust", "Plantae", "Rhodophyta", "Florideophyceae", "Corallinales", "Corallinaceae", "Coralline", NA_character_, "manual fix - gbif",
  "Corella sp.", "Animalia", "Chordata", "Ascidiacea", "Phlebobranchia", "Corellidae", "Corella", NA_character_, "manual fix - worms",
  "Corophiinae", "Animalia", "Arthropoda", "Malacostraca", "Amphipoda", "Corophiidae", NA_character_, NA_character_,"manual fix - worms",
  "Coryne", "Animalia", "Cnidaria", "Hydrozoa", "Anthoathecata", "Corynidae", "Coryne", NA_character_, "manual fix - gbif",
  "Crisia Occidentalis", "Animalia", "Bryozoa", "Stenolaemata", "Cyclostomatida", "Crisiidae", "Crisia", "Crisia occidentalis", "manual fix - worms",
  "Cryptomonas gyropyrenoidosa", "Chromista", "Cryptophyta", "Cryptophyceae", "Cryptomonadales", "Cryptomonadaceae", "Cryptomonas", "Cryptomonas gyropyrenoidosa","manual fix - gbif",
  "Cunea russae", "Protozoa", "Amoebozoa", "Discosea", "Dactylopodida", "Paramoebidae", "Cunea", "Cunea russae", "manual fix - gbif/ncbi",
  "Cunea thuwala", "Protozoa", "Amoebozoa", "Discosea", "Dactylopodida", "Paramoebidae", "Cunea", "Cunea thuwala", "manual fix - gbif/ncbi",
  "Cunea profundata", "Protozoa", "Amoebozoa", "Discosea", "Dactylopodida", "Paramoebidae", "Cunea", "Cunea profundata", "manual fix - gbif/ncbi",
  "Cunea", "Protozoa", "Amoebozoa", "Discosea", "Dactylopodida", "Paramoebidae", "Cunea", NA_character_, "manual fix - gbif/ncbi",
  "Cyanea", "Animalia", "Cnidaria", "Scyphozoa", "Semaeostomeae", "Cyaneidae", "Cyanea", NA_character_, "manual fix - worms",
  # Ds
  "Desmarestia", "Chromista", "Ochrophyta", "Phaeophyceae", "Desmarestiales", "Desmarestiaceae", "Desmarestia", NA_character_, "manual fix - gbif",
  "Dictyosiphon", "Chromista", "Animalia", "Phaeophyceae", "Ectocarpales", "Chordariaceae", "Dictyosiphon", NA_character_, "manual fix - worms",
  "Dictyota", "Chromista", "Ochrophyta", "Phaeophyceae", "Dictyotales", "Dictyotaceae", "Dictyota", NA_character_, "manual fix - gbif",
  "Diplosoma", "Animalia", "Chordata", "Ascidiacea", "Aplousobranchia", "Didemnidae", "Diplosoma", NA_character_, "manual fix - gbif",
  "Diptera", "Animalia","Arthropoda", "Insecta", "Diptera", NA_character_,  NA_character_, NA_character_, "manual fix - gbif",
  "Dorvillea schistomeringos sp.", "Animalia", "Annelida", "Polychaeta", "Eunicida", "Dorvilleidae", "Dorvillea", "Dorvillea schistomeringos", "manual fix - gbif",
  "Dorvillea schistomeringos longicornis", "Animalia", "Annelida", "Polychaeta", "Eunicida", "Dorvilleidae", "Dorvillea", "Dorvillea schistomeringos", "manual fix - gbif",
  "Doto", "Animalia", "Mollusca", "Gastropoda", "Nudibranchia", "Dotidae", "Doto", NA_character_, "manual fix - gbif",
  "Doto form", "Animalia", "Mollusca", "Gastropoda", "Nudibranchia", "Dotidae", "Doto", NA_character_, "manual fix - gbif",
  # Es
  "Elachista", "Chromista", "Ochrophyta", "Phaeophyceae", "Ectocarpales", "Chordariaceae","Elachista", NA_character_, "manual fix - gbif",
  "Elamenopsis kempi", "Animalia", "Arthropoda", "Malacostraca","Decapoda", "Hymenosomatidae", "Elamenopsis", "Elamenopsis kempi", "manual fix",
  "Electra", "Animalia", "Bryozoa", "Gymnolaemata", "Cheilostomatida", "Electridae", "Electra", NA_character_, "manual fix - gbif",
  "Eteone", "Animalia", "Annelida", "Polychaeta", "Phyllodocida","Phyllodocidae", "Eteone", NA_character_, "manual fix - worms",
  "Eusyllinae", "Animalia", "Annelida", "Polychaeta", "Phyllodocida", "Syllidae", NA_character_, NA_character_, "manual fix - worms",# it's a subfamily
  "Exogoninae", "Animalia", "Annelida", "Polychaeta", "Phyllodocida", "Syllidae", NA_character_, NA_character_, "manual fix - worms",# it's a subfamily
  # Fs
  "Fissipedicella orientalis", "Chromista", "Animalia", "Phaeophyceae", "Ralfsiales", "Ralfsiaceae", "Fissipedicella", "Fissipedicella orientalis", "manual fix - worms",
  "Flabellina", "Animalia", "Mollusca", "Gastropoda", "Nudibranchia", "Flabellinidae", "Flabellina", NA_character_, "manual fix - worms",
  # Gs
  "Gammarus salinus", "Animalia",  "Arthropoda", "Malacostraca", "Amphipoda", "Gammaridae", "Gammarus", "Gammarus salinus", "manual fix - worms",
  "Gnorisphaeroma oregonesis", "Animalia", "Arthropoda", "Malacostraca", "Isopoda", "Sphaeromatidae", "Gnorimosphaeroma", "Gnorimosphaeroma oregonensis", "manual fix - gbif",
  "Grantia", "Animalia", "Porifera", "Calcarea", "Leucosolenida", "Grantiidae", "Grantia", NA_character_, "manual fix - gbif",
  "Globisporangium tenuihyphum", "Chromista", "Oomycota", "Peronosporea", "Peronosporales", "Pythiaceae", "Globisporangium", "Globisporangium tenuihyphum","manual fix - gbif",
  # Hs
  "Haptoglossa", "Chromista", "Oomycota", "Peronosporea", "Haptoglossales", "Haptoglossaceae", "Haptoglossa",NA_character_, "manual fix - gbif",
  "Hanseniella", "Animalia","Arthropoda", "Symphyla",NA_character_, "Scutigerellidae", "Hanseniella",NA_character_, "manual fix - gbif",
  "Heterococcus caespitosus", "Chromista","Gyrista","Xanthophyceae","Tribonematales","Heteropediaceae","Heterococcus","Heterococcus caespitosus", "manual fix - ncbi", 
  "Hexanauplia", "Animalia","Arthropoda", "Hexanauplia", NA_character_,NA_character_, NA_character_, NA_character_, "manual fix - worms",
  "Histiona aroides", "Protozoa", "Loukozoa", "Jakobea", "Jakobida", "Histionidae", "Histiona", "Histiona aroides", "manual fix - gbif",
  "Hyalidae", "Animalia","Arthropoda", "Malacostraca", "Amphipoda", "Hyalidae", NA_character_, NA_character_, "manual fix - worms",
  # Is
  "Indohya lynbeazleyae", "Animalia",  "Arthropoda", "Arachnida", "Pseudoscorpiones", "Hyidae", "Indohya", "Indohya lynbeazleyae", "manual fix - gbif",
  "Isopoda", "Animalia","Arthropoda", "Malacostraca", "Isopoda", NA_character_, NA_character_, NA_character_, "manual fix - worms",
  # Js
  "Jakoba bahamiensis", "Protozoa", "Loukozoa", "Jakobea", "Jakobida", "Jakobidae", "Jakoba", "Jakoba bahamiensis", "manual fix - gbif",
  "Jania verrucosa", "Plantae", "Rhodophyta", "Florideophyceae", "Corallinales", "Corallinaceae", "Jania", NA_character_, "manual fix - worms", # uncertain status
  "Januini", "Animalia", "Annelida", "Sabellida", "Serpulidae", NA_character_, NA_character_,NA_character_, "manual fix - worms",
  # Ks
  "Kamptozoa", "Animalia", "Entoprocta", NA_character_, "Solitaria", "Barentsiidae", "Urnatella", "Urnatella gracilis", "manual fix - worms",
  "Korotnevella heteracantha", "Protozoa", "Amebozoa", "Lobosa", "Amoebida", "Paramoebidae", "Korotnevella", "Korotnevella heteracantha","manual fix - gbif",
  # Ls
  "Leiosolenus curtus", "Animalia",  "Mollusca", "Bivalvia", "Mytilida", "Mytilidae", "Leisolenus", "Leiosolenus curtus", "manual fix - gbif",
  "Lemonias zygia", "Animalia", "Arthropoda", "Insecta", "Lepidoptera", "Riodinidae", "Lemonias", "Lemonias zygia", "manual fix - inaturalist",
  "Leucocytozoon grallariae", "Chromista","Apicomplexa","Aconoidasida","Haemosporida","Leucocytozoidae","Leucocytozoon","Leucocytozoon grallariae", "manual fix - gbif/ncbi",
  "Leucocytozoon neotropicalis", "Chromista","Apicomplexa","Aconoidasida","Haemosporida","Leucocytozoidae","Leucocytozoon","Leucocytozoon neotropicalis", "manual fix - gbif/ncbi",
  "Limpet", "Animalia",  "Mollusca", "Gastropoda", NA_character_, NA_character_,  NA_character_,  NA_character_,  "manual fix",
  # Ms
  "malaiseGraci01 Malaise4560", "Animalia", "Arthropoda", "Insecta", "Lepidoptera", "Gracillariidae", NA_character_, NA_character_, "manual fix - bold",
  "Macrobranchium", "Animalia",  "Arthropoda", "Malacostraca", "Decapoda", "Palaemonidae",  "Macrobranchium", NA_character_, "manual fix",
  "Marophrys", "Protozoa",NA_character_,"Centroplasthelida","Acanthocystida",NA_character_,"Marophrys",NA_character_, "manual fix - gbif/ncbi",
  "Mazzaella canaliculata", "Plantae", "Rhodophyta", "Florideophyceae", "Gigartinales", "gigartinaceae", "Mazzaella", "Mazzaella canaliculata","manual fix - gbif",
  "Monosiga brevicollis", "Protozoa", "Choanozoa", "Choanoflagellatea", "Choanoflagellida", "Codosiginidae", "Monosiga","Monosiga brevicollis", "manual fix - gbif",
  "Morulina",  "Animalia",  "Arthropoda", "Collembola", "Poduromorpha", "Neanuridae", "Morulina", NA_character_, "manual fix - NCBI",
  "Microcosmus sp.", "Animalia", "Chordata", "Ascidiacea", "Stolidobranchia", "Pyuridae", "Microcosmus", NA_character_, "manual fix - worms",
  "Microcosmus ", "Animalia", "Chordata", "Ascidiacea", "Stolidobranchia", "Pyuridae", "Microcosmus", NA_character_, "manual fix - worms",
  "Microcosmus", "Animalia", "Chordata", "Ascidiacea", "Stolidobranchia", "Pyuridae", "Microcosmus", NA_character_, "manual fix - worms",
  "Microheliella maris", "Protozoa", NA_character_,NA_character_,NA_character_,NA_character_,"Microheliella","Microheliella maris","manual fix - publication",# https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2022.1030570/full
  "Mucochytrium quahogii", "Chromista","Bigyra","Bigyra","Thraustochytrida","Thraustochytriaceae","Mucochytrium","Mucochytrium quahogii", "manual fix - ncbi",
  "Myxicola Infundibulum", "Animalia", "Annelida", "Polychaeta", "Sabellida", "Sabellidae", "Myxicola", "Myxicola Infundibulum", "manual fix - gbif",
  # Ns
  "Navicula", "Chromista", "Ochrophyta", "Bacillariophyceae", "Naviculales", "Naviculaceae", "Navicula",  NA_character_, "manual fix - gbif",
  "Nereis", "Animalia", "Annelida", "Polychaeta", "Phyllodocida", "Nereididae", "Nereis", NA_character_, "manual fix - gbif",
  "Nitzschia anatoliensis", "Chromista", "Ochrophyta", "Bacillariophyceae","Bacillariales", "Bacillariaceae", "Nitzschia",  NA_character_, "manual fix - ncbi", 
  "Nitzschia", "Chromista", "Ochrophyta", "Bacillariophyceae","Bacillariales", "Bacillariaceae", "Nitzschia",  NA_character_, "manual fix - gbif",
  # Os
  "Obelia", "Animalia", "Cnidaria", "Hydrozoa", "Leptothecata", "Campanulariidae", "Obelia", NA_character_, "manual fix - worms",
  "Ophelia", "Animalia", "Annelida", "Polychaeta", NA_character_, "Opheliidae", "Ophelia",NA_character_, "manual fix - gbif",
  "Ovalopodium rosalinum", "Protozoa", "Amoebozoa", "Discosea", "Himatismenida", "Cochliopodiidae", "Ovalopodium", "Ovalopodium rosalinum", "manual fix - gbif/ncbi",
  "Oligochaeta", "Animalia", "Annelida", "Clitellata", NA_character_,  NA_character_, NA_character_, NA_character_,"manual fix - worms",
  # Ps
  "Patinella", "Animalia", "Bryozoa", "Stenolaemata", "Cyclostomatida", "Lichenoporidae", "Patinella" , NA_character_, "manual fix - gbif/worms",
  "Pectinaria" , "Animalia", "Annelida", "Polychaeta", NA_character_, "Pectinariidae", "Pectinaria" ,NA_character_, "manual fix - gbif",
  "Pelagia", "Animalia", "Cnidaria", "Scyphozoa", "Semaeostomeae", "Pelagiidae", "Pelagia", NA_character_, "manual fix - gbif",
  "Pholoe", "Animalia", "Annelida", "Polychaeta", "Phyllodocida","Sigalionidae", "Pholoe", NA_character_, "manual fix - gbif/worms",
  "Phoridae",  "Animalia", "Arthropoda","Insecta",  "Diptera", "Phoridae",  NA_character_, NA_character_, "manual fix - gbif",
  "Phyllodoce", "Animalia", "Annelida", "Polychaeta", "Phyllodocida","Phyllodocidae", "Phyllodoce", NA_character_, "manual fix - gbif",
  "Physa", "Animalia", "Mollusca", "Gastropoda", NA_character_, "Physidae", "Physa", NA_character_, "manual fix - ncbi",
  "Pinnularia", "Chromista", "Ochrophyta", "Bacillariophyceae", "Naviculales", "Pinnulariaceae", "Pinnularia", NA_character_, "manual fix - ncbi",
  "Pileolaria sp.", "Animalia", "Annelida", "Polychaeta", "Sabellida", "Serpulidae", "Pileolaria", NA_character_, "manual fix - worms",
  "Planopodium haveli",  "Protozoa", "Amoebozoa", "Discosea", "Himatismenida", NA_character_, "Planopodium",  "Planopodium haveli", "manual fix - paper", #https://pubmed.ncbi.nlm.nih.gov/33709902/
  "Pleonosporium australicum",  "Plantae", "Rhodophyta", "Florideophyceae", "Ceramiales", "Wrangeliaceae", "Pleonosporium", "Pleonosporium australicum","manual fix - gbif",
  "Plumularioidea", "Animalia", "Cnidaria", "Hydrozoa", "Leptothecata",  NA_character_, NA_character_, NA_character_, "manual fix - gbif/worms",
  "Polychaeta", "Animalia", "Annelida", "Polychaeta", NA_character_,  NA_character_, NA_character_, NA_character_, "manual fix - gbif/worms",
  "Polydora", "Animalia", "Annelida", "Polychaeta", NA_character_, "Spionidae", "Polydora", NA_character_, "manual fix - gbif",
  "Polysiphonia", "Plantae", "Rhodophyta", "Florideophyceae", "Ceramiales", "Rhodomelaceae", "Polysiphonia",NA_character_, "manual fix - gbif",
  "Procephalotrix spiralis", "Animalia", "Nemertea", "Paleonemertea", "Archinemertea", "Cephalotrichidae", "Cephalotrix", "Cephalotrix spiralis", "manual fix - worms",
  # Qs
  # Rs
  "Radix", "Animalia", "Mollusca", "Gastropoda", NA_character_,"Lymnaeidae" ,"Radix", NA_character_,"manual fix - gbif",
  "Ramipedicella longicellularis", "Chromista", "Animalia", "Phaeophyceae", "Ralfsiales", "Ralfsiaceae","Ramipedicella", NA_character_,"manual fix - gbif",
  "Reclinomonas americana", "Protozoa", "Loukozoa", "Jakobea", "Jakobida", "Histionidae", "Reclinomonas", "Reclinomonas americana","manual fix - gbif",
  # Ss
  "Sacoglossa", "Animalia", "Mollusca", "Gastropoda", NA_character_, NA_character_, NA_character_, NA_character_, "manual fix - worms",
  "Sagitta", "Animalia", "Chaetognatha", "Sagittoidea", "Aphragmophora", "Sagittidae", "Sagitta", NA_character_,"manual fix - gbif",
  "Salvatoria brevipharyngea", "Animalia", "Annelida", "Polychaeta", "Phyllodocida","Syllidae", "Salvatoria", "Salvatoria brevipharyngea", "manual fix - worms",
  "Satsuma", "Animalia", "Mollusca", "Gastropoda", "Stylommatophora", "Camaenidae","Satsuma", NA_character_, "manual fix - gbif",
  "Schizostauron trachyderma", "Chromista", "Ochrophyta", "Bacillariophyceae", "Naviculales", "Stauroneidaceae", "Schizostauron", "Schizostauron trachyderma","manual fix - gbif",
  "Sciaridae", "Animalia", "Arthropoda", "Insecta", "Diptera", "Sciaridae", NA_character_, NA_character_, "manual fix",
  "Scytodes", "Animalia","Arthropoda","Arachnida","Araneae","Scytodidae","Scytodes", NA_character_,"manual fix - ncbi",
  "Spio", "Animalia", "Annelida", "Polychaeta", NA_character_, "Spionidae", "Spio", NA_character_,"manual fix - gbif",
  "Stratonice succinea", "Animalia", "Annelida", "Polychaeta", "Phyllodocida", "Nereididae", "Alitta", "Alitta succinea", "manual fix - worms",
  "Symmachia probetor", "Animalia","Arthropoda","Insecta","Lepidoptera","Riodinidae","Symmachia","Symmachia probetor", "manual fix - ncbi",
  "Symmachia tricolor", "Animalia","Arthropoda","Insecta","Lepidoptera","Riodinidae","Symmachia","Symmachia tricolor", "manual fix - ncbi",
  # Ts
  "Terebellidae", "Animalia", "Annelida", "Polychaeta", "Terebellida", "Terebellidae", NA_character_, NA_character_, "manual fix - worms",
  "Tetracladium", "Fungi", "Ascomycota", "Leotiomycetes", "Heleotiales", "Helotiaceae", "Tetracladium", NA_character_,"manual fix - gbif",
  "Tilopteridalean", "Chromista", "Ochrophyta", "Phaeophyceae", NA_character_, NA_character_, NA_character_, NA_character_, "manual fix - ncbi",
  "Trichonotus setiger", "Animalia", "Chordata", "Teleostei", "Gobiiformes", "Trichonotidae", "Trichonotus", "Trichonotus setiger","manual fix - worms",
  "Tubulanus", "Animalia", "Nemertea", "Palaeonemertea", "Tubuliformes", "Tubulanidae", "Tubulanus", NA_character_, "manual fix - worms",
  "tubulipora sp.", "Animalia", "Bryozoa", "Stenolaemata", "Cyclostomatida", "Tubuliporidae", "Tubulipora", NA_character_,"manual fix - worms",
  "Tunicata", "Animalia", "Chordata", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "manual fix - worms",
  # Us
  "Ulva", "Plantae", "Chlorophyta", "Ulvophyceae", "Ulvales", "Ulvaceae", "Ulva", NA_character_, "manual fix - gbif",
  # Vs
  "Verhoeffiella", "Animalia", "Arthropoda", "Collembola", "Entomobryomorpha", "Entomobryidae", "Verhoeffiella", NA_character_, "manual fix - gbif",
  "Vermiviatum covidum", "Animalia", "Platyhelminthes", NA_character_, "Tricladida", "Geoplanidae", "Vermiviatum","Vermiviatum covidum", "manual fix - gbif",
  "Vexillifera", "Protozoa", "Amoebozoa", "Discosea", "Dactylopodida", "Vexilliferidae", "Vexillifera", NA_character_, "manual fix - worms",
  "Vischeria", "Chromista", "Ochrophyta", "Eutigmatophyceae", "Eustigmatales", "Chlorobotryaceae", "Vischeria", NA_character_, "manual fix - gbif and ncbi",
  # Ws
  # Xs
  # Ys
  # Zs
)