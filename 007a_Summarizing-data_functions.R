# ------------------------------------------------------------------------------
# -------------- Summarizing data and plotting functions -----------------------
# ------------------------------------------------------------------------------

# AUTHOR: Paula Pappalardo - pappalardop@si.edu

# Custom taxa for all good metazoans
custom_taxa_all <- c("Annelida" =  "#71284E",
                     "Bryozoa" = "#8E4C76",
                     "Chaetognatha" =  "#297F9E",
                     "Chordata:Teleostei" = "#68ADC7",
                     "Cnidaria" = "#B0DAE7",
                     "Echinodermata" = "#C199C0",
                     "Mollusca:Bivalvia" = "#238B45",
                     "Mollusca:Gastropoda" = "#74C476",
                     "Mollusca:Polyplacophora" = "#A1D99B",
                     "Nemertea" = "#D8C1D7",
                     "Phoronida" = "#8E8E8E",
                     "Porifera" = "#A86F9A",
                     "Platyhelminthes" = "#FEEDA0",
                     "Ctenophora" = "#5A5A5A",
                     "Arthropoda:Branchiopoda" = "#5C2D13",
                     "Arthropoda:Copepoda:Cyclopoida" = "#FEC44F",
                     "Arthropoda:Copepoda:Calanoida" = "#FF7F00",
                     "Arthropoda:Copepoda:Harpacticoida" =  "#FE9929",
                     "Arthropoda:Copepoda:unresolved" = "#FDAE6B",
                     "Arthropoda:Malacostraca" = "#92634C",
                     "Arthropoda:Thecostraca" = "#AE8A78",
                     "Arthropoda:Ostracoda" = "#C8AD9F",
                     "Arthropoda:unresolved" = "#76462D",
                     "Rotifera" = "#C7E9C0")

# niceplot format
niceplot <- theme_bw(base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "bottom")

addCustomTaxa <- function(mydf){
  # add custom taxa for customized color palette
  mydf_ed <- mydf %>% 
    mutate(custom_taxa = case_when(phylum == "Annelida" ~ "Annelida",
                                   phylum == "Bryozoa" ~ "Bryozoa",
                                   phylum == "Chaetognatha" ~ "Chaetognatha",
                                   phylum == "Chordata" ~ "Chordata:Teleostei",
                                   phylum == "Cnidaria" ~ "Cnidaria",
                                   phylum == "Echinodermata" ~ "Echinodermata",
                                   phylum == "Mollusca" & class == "Bivalvia" ~ "Mollusca:Bivalvia",
                                   phylum == "Mollusca" & class == "Gastropoda" ~ "Mollusca:Gastropoda",
                                   phylum == "Mollusca" & class == "Polyplacophora" ~ "Mollusca:Polyplacophora",
                                   phylum == "Nemertea" ~ "Nemertea",
                                   phylum == "Phoronida" ~ "Phoronida",
                                   phylum == "Porifera" ~ "Porifera",
                                   phylum == "Platyhelminthes" ~ "Platyhelminthes",
                                   phylum == "Ctenophora" ~ "Ctenophora",
                                   phylum == "Arthropoda" & class == "Branchiopoda" ~ "Arthropoda:Branchiopoda",
                                   phylum == "Arthropoda" & class == "Copepoda" & order == "Cyclopoida" ~ "Arthropoda:Copepoda:Cyclopoida",
                                   phylum == "Arthropoda" & class == "Copepoda" & order == "Calanoida" ~ "Arthropoda:Copepoda:Calanoida",
                                   phylum == "Arthropoda" & class == "Copepoda" & order == "Harpacticoida" ~ "Arthropoda:Copepoda:Harpacticoida",
                                   phylum == "Arthropoda" & class == "Copepoda" & is.na(order) ~"Arthropoda:Copepoda:unresolved",
                                   phylum == "Arthropoda" & class == "Malacostraca" ~ "Arthropoda:Malacostraca",
                                   phylum == "Arthropoda" & class == "Thecostraca" ~ "Arthropoda:Thecostraca",
                                   phylum == "Arthropoda" & class == "Ostracoda" ~ "Arthropoda:Ostracoda",
                                   phylum == "Arthropoda" & is.na(class) ~ "Arthropoda:unresolved",
                                   phylum == "Rotifera" ~ "Rotifera"))
}        