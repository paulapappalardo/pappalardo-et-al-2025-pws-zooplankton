# ------------------------------------------------------------------------------
# -------------- Phyloseq data analysis and plotting functions -----------------------
# ------------------------------------------------------------------------------

# AUTHOR: Paula Pappalardo - pappalardop@si.edu

#----------------------------------FOR PLOTS-----------------------------------
library(ggplot2)
library(extrafont)

# The color blind friendly palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9",
               "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7")

# Load fonts for PDF output
loadfonts(device = "pdf")

#----------------------------------Functions-----------------------------------

# fix factor levels to be in chronological (or desired) order for better plots
fixFactorsOrder <- function(myphyloseq){
  sample_data(myphyloseq)$month = factor(sample_data(myphyloseq)$month,
                                         levels = c("April", "May", "June", "July", "August", "September"),
                                         labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))
  #sample_data(myphyloseq)$Week = factor(sample_data(myphyloseq)$Week, levels = c("Week1","Week2","Week3","Week4","Week5", "Week6", "Week7", "Week8", "Week9", "Week10", "Week11", "Week12", "Week13", "Week14", "Week15", "Week16", "Week17", "Week18", "Week19", "Week20", "Week21"))
  # sample_data(myphyloseq)$spawning = factor(sample_data(myphyloseq)$spawning,
  #                                                 levels = c("Onset", "Peak", "Decrease", "Past"))
  sample_data(myphyloseq)$year = factor(sample_data(myphyloseq)$year,
                                        levels = c("2017", "2018", "2019", "2021"))
  return(myphyloseq)
}

addPlanktonGroupToTaxTable <- function(myphyloseq){
  # Adds plankton group to simple tax table
  #
  library(taxize)
  # extract tax table as dataframe
  updatedtaxtable <- tax_table(myphyloseq) %>%
    as.data.frame() %>% 
    dplyr::select(kingdom, phylum, class, order, family, genus, species, sciname) %>% 
    addPlanktonGroup() %>% 
    dplyr::relocate(plankton_group, .after = kingdom)
  # convert taxtable to phyloseq format
  mat <- updatedtaxtable %>% as.matrix(.)
  rownames(mat) <- rownames(updatedtaxtable)
  taxa_updated <- tax_table(mat)
  # update the phyloseq object tax table
  tax_table(myphyloseq) <- taxa_updated
  # return updated phyloseq
  myphyloseq
}

# The "addPlanktonGroup" function was downloaded from Pappalardo (2025)
# https://doi.org/10.5281/zenodo.17467576

addPlanktonGroup <- function(mydf){
  # Adds a column with the Plankton group classification 
  #
  # Args:
  #   dataframe including columns with higher taxonomy, must include
  #   columns for phylum, class, order, family, genus, and species
  # Returns:
  #   The same dataframe with an additional column named plankton_group.
  #   Plankton groups were Meroplankton and Holoplankton, but I also included
  #   Benthic category, Holo/Mero when uncertain, and highlighted taxa that are parasites.
  #
  # Packages required:
  library(taxize)
  library(dplyr)
  #
  # GBIF taxonomic backbone does not have class Actinopteri,
  # we can get orders using taxize in case we need to separate fishes
  # same with suborder Hyperiidea to separate holoplanktonic amphipods
  
  # fishes
  actinopteri <- downstream("Actinopteri", db ='worms', downto = "Order")
  acti_orders <- actinopteri[[1]]$name
  rm(actinopteri)
  # hyperiids
  hyperiids <- downstream("Hyperiidea", db ='worms', downto = "Family")
  hype_families <- hyperiids[[1]]$name
  # Create the new column "plankton_group" based on taxonomic information
  mydf_ed <- mydf %>%
    arrange(phylum, class, order) %>%
    mutate(plankton_group = NA) %>% 
    mutate(plankton_group = case_when(
      #
      # ----Phylum Annelida----
      #
      # Polychaetes are tricky, they include a mix of strategies 
      #class == "Polychaeta"  ~ "Holo/Mero", # below we have specific assignments for smaller taxonomic groups
      class == "Polychaeta" & genus %in% c("Marphysa", "Cirriformia", "Diopatra", "Hemipodia", "Pholoides", "Pholoe", "Glycinde", "Protodorvillea", "Sthenelais", "Dorvillea", "Protodrilus", "Polygordius", "Glycera", "Ophryotrocha", "Pectinaria", "Ophelia") ~ "Meroplankton",
      class == "Polychaeta" & genus %in% c("Dinophilus") ~ "Holoplankton",
      class == "Polychaeta" & species %in% c("Armandia brevis", "Schistomeringos longicornis", "Capitella teleta", "Scoloplos armiger", "Hydroides elegans", "Cistenides granulata", "Terebellides horikoshii", "Ficopomatus enigmaticus") ~ "Meroplankton",
      class == "Polychaeta" & species %in% c("Naineris dendritica", "Dimorphilus gyrociliatus", "Ctenodrilus serratus", "Terebellides stroemii") ~ "Benthic",
      # 6 Glycera species in InvertTraits with planktonic larvae
      # Fam Cirritulidae are benthic: https://en.wikipedia.org/wiki/Cirratulidae
      # Genus Cirriformia has 3 members with planktonic larva in InvertTraits
      # Hydroides elegans planktotrophic larvae: https://invasions.si.edu/nemesis/species_summary/68295
      # 5 Pholoe in invertTraits with planktonic larva
      # Terebellides horikoshii is best guess based on their benthic larvae and reports of lecithotrophic larva: https://www.mdpi.com/1424-2818/13/2/60 
      # 3 Ophelia in InvertTraits with planktonic non feeding larva, benthic adult
      # 2 Pectinaria in InvertTraits with planktonic larva, benthic adult
      # Cistenides granulata confirmed trocophora by Shanks book of Pacific invertebrate larvae
      # in InvertTraits, most Ophryotrocha seems to have planktonic-nonfeeding larva (except O. vivipara) so I added as meroplankton
      # Scoloplos armiger is reported as exhibiting poecilogony in InvertTraits, adding it as meroplankton
      # Terebellides stroemii is aplanktonic in InvertTraits, so was classified as Benthics
      # Dimorphilus gyrociliatus is aplanktonic larva based on InvertTraits and meiofauna based on https://www.martinduranlab.com/dimorphilus
      # Schistomeringos longicornis has "benthic" in SeaLifeBase and free spawing larva according to  https://scholarsbank.uoregon.edu/xmlui/bitstream/handle/1794/6123/6.pdf?sequence=16&isAllowed=y
      # Capitella teleta is planktonic-notfeeding larvae in InvertTraits, and benthic according https://en.wikipedia.org/wiki/Capitella_teleta
      # Ctenodrilus serratus Westheide et al 2003 says it is meiofauna without planktonic larva
      # Pholoides have planktonic lecithotrophic larva according to: https://www.invertebase.org/portal/taxa/index.php?taxon=162460&clid=16
      # Glycinde have planktonic larva according to: https://scholarsbank.uoregon.edu/xmlui/bitstream/handle/1794/6123/6.pdf?sequence=16&isAllowed=y
      # InverTraits also has two Glycinde species as planktonic larva
      # It seems that Dinophilus are free living found in tidal pools, and have direct development (info agrees with 1 species from InvertTraits):
      # https://www.journals.uchicago.edu/doi/epdf/10.2307/1535836
      # https://en.wikipedia.org/wiki/Dinophilus_taeniatus
      # Protodorvillea seems to have representatives with pelagic larva: https://www.marlin.ac.uk/habitats/detail/1115/protodorvillea_kefersteini_and_other_polychaetes_in_impoverished_circalittoral_mixed_gravelly_sand
      # Sthenelais has at least two other species with planktonic development according https://link.springer.com/article/10.1007/BF02366201
      # Another Dorvillea from California has planktonic-nonfeeding larva (Blake 1975 for D. rudolphi)
      # Plate and Husemann 1994_Identification guide to the planktonic polychaete larvae had a lot of info to identify meroplanktonic species
      phylum == "Annelida" & class == "Polychaeta" & family == "Amphinomidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Thalassematidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Alciopidae" ~ "Holoplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Chrysopetalidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Lopadorrhynchidae" ~ "Holoplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Tomopteridae" ~ "Holoplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Typhloscolecidae" ~ "Holoplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Oweniidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Capitellidae" ~ "Meroplankton", # adults in muddy sediments, 9 species in InvertTraits with planktonic larva
      phylum == "Annelida" & class == "Polychaeta" & family == "Sabellidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Sabellariidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Poecilochaetidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Nephtyidae" ~ "Meroplankton", # InvertTraits has all 9 species with info as having planktonic larvae
      phylum == "Annelida" & class == "Polychaeta" & family == "Dorvilleidae" ~ "Meroplankton", # InvertTraits has 10 of 11 species with info as having planktonic larvae
      phylum == "Annelida" & class == "Polychaeta" & family == "Spionidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Ampharetidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Terebellidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Chaetopteridae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Magelonidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Polynoidae" ~ "Meroplankton",
      phylum == "Annelida" & class == "Polychaeta" & family == "Hesionidae" ~ "Meroplankton",# from Paula's InvertTraits database
      phylum == "Annelida" & class == "Polychaeta" & genus == "Leitoscoloplos" ~ "Meroplankton",# from Paula's InvertTraits database
      phylum == "Annelida" & class == "Polychaeta" & family == "Phyllodocidae" ~ "Meroplankton", #https://en.wikipedia.org/wiki/Phyllodocidae
      phylum == "Annelida" & class == "Polychaeta" & family == "Nereididae" ~ "Meroplankton", #https://en.wikipedia.org/wiki/Nereididae
      phylum == "Annelida" & class == "Polychaeta" & family == "Syllidae" ~ "Meroplankton", #https://en.wikipedia.org/wiki/Syllidae
      phylum == "Annelida" & class == "Polychaeta" & family == "Urechidae" ~ "Meroplankton", #https://en.wikipedia.org/wiki/Urechis_caupo
      phylum == "Annelida" & class == "Clitellata" & family == "Naididae" ~ "Benthic", # book series "Freshwater Invertebrates of the United States." third edition of Volume II, "Oligochaeta" (2012), by Brinkhurst and Marchese
      #
      # ----Phylum Mollusca----
      #
      phylum == "Mollusca" & class== "Gastropoda" & species == "Rictaxis punctocaelatus" ~ "Meroplankton",
      phylum == "Mollusca" & class== "Gastropoda" & species == "Physella acuta" ~ "benthic - freshwater",
      phylum == "Mollusca" & class== "Gastropoda" & family == "Pyramidellidae" ~ "Meroplankton - parasite of invertebrates",
      phylum == "Mollusca" & class== "Gastropoda" & (family %in% c("Plakobranchidae", "Hermaeidae", "Limapontiidae", "Lottiidae"))  ~ "Meroplankton",
      phylum == "Mollusca" & class== "Gastropoda" & (family  %in% c("Atlantidae", "Pterotracheidae", "Carinariidae")) ~ "Holoplankton",
      phylum == "Mollusca" & class== "Gastropoda" & order!= "Pteropoda" ~ "Meroplankton",
      # 3 Tonicellidae in InvertTraits with planktonic-nonfeeding larvae
      phylum == "Mollusca" & class== "Polyplacophora" & family %in% c("Mopaliidae", "Tonicellidae") ~ "Meroplankton",
      phylum == "Mollusca" & order== "Pteropoda" ~ "Holoplankton",
      phylum == "Mollusca" & class== "Bivalvia" ~ "Meroplankton",
      phylum == "Mollusca" & class== "Cephalopoda" ~ "Meroplankton",
      phylum == "Mollusca" & class== "Scaphopoda" ~ "Meroplankton",
      #
      # ----Phylum Cnidaria----
      #
      phylum == "Cnidaria" & class == "Hydrozoa" & family == "Geryoniidae" ~ "Holoplankton",
      phylum == "Cnidaria" & class == "Hydrozoa" & order == "Trachymedusae" ~ "Holoplankton",
      phylum == "Cnidaria" & class == "Hydrozoa" & order == "Narcomedusae" ~ "Holoplankton",
      phylum == "Cnidaria" & class == "Hydrozoa" & order == "Siphonophorae" ~ "Holoplankton",
      phylum == "Cnidaria" & class == "Anthozoa" ~ "Meroplankton",
      phylum == "Cnidaria" & class == "Hydrozoa" ~ "Meroplankton",
      phylum == "Cnidaria" & class == "Scyphozoa" ~ "Meroplankton",
      phylum == "Cnidaria" & order == "Stauromedusae" ~ "Meroplankton", # looks like benthic adults, crawling larva, free spawing eggs https://en.wikipedia.org/wiki/Stauromedusae
      #
      # ----Phylum Chordata----
      #
      # TODO: think about your samples and modify the fish classification if needed
      phylum == "Chordata" & order %in% acti_orders ~ "Holoplankton", # but if this is from plankton samples where it could be eggs/larvae it could be ichtyoplankton 
      phylum == "Chordata" & class== "Appendicularia" ~ "Holoplankton",
      phylum == "Chordata" & class== "Thaliacea" ~ "Holoplankton",
      phylum == "Chordata" & class== "Ascidiacea" ~ "Meroplankton",
      phylum == "Chordata" & class== "Leptocardii" ~ "Benthic",
      #
      # ----Phylum Arthropoda----
      #
      phylum == "Arthropoda" & class== "Maxillopoda" & order == "Pedunculata" & species == "Pollicipes polymerus" ~ "Meroplankton",
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Decapoda" &
        (family == "Penaeidae" | family == "Palaemonidae" | family == "Hippolytidae" |
           family == "Pandalidae" | family == "Alpheidae") ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Decapoda" &
        !(family == "Penaeidae" | family == "Palaemonidae" | family == "Hippolytidae" |
            family == "Pandalidae" | family == "Alpheidae") ~ "Meroplankton",
      phylum == "Arthropoda" & class== "Maxillopoda" & order == "Sessilia"~ "Meroplankton",
      phylum == "Arthropoda" & class== "Thecostraca" & order == "Sessilia"~ "Meroplankton",
      phylum == "Arthropoda" & class== "Thecostraca" & order == "Balanomorpha"~ "Meroplankton",
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Euphausiacea" ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Tanaidacea" ~ "Benthic",
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Cumacea" ~ "Benthic",
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Isopoda" ~ "Benthic", 
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Mysida" ~ "Benthic", 
      phylum == "Arthropoda" & class== "Malacostraca" &  order == "Stomatopoda" ~ "Meroplankton",
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Amphipoda" & family == "Corophiidae" ~ "Benthic", # https://www.mapress.com/zootaxa/2009/f/zt02260p379.pdf
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Amphipoda" & family == "Ischyroceridae" ~ "Benthic", # https://biodiversity.org.au/afd/taxa/Ischyroceridae
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Amphipoda" & family %in% hype_families ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Amphipoda"  ~ "Benthic",
      phylum == "Arthropoda" & class== "Malacostraca" & order == "Lophogastrida" ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Copepoda" & species %in% c("Bomolochus cuneatus", "Mytilicola orientalis",
                                                                   "Holobomolochus spinulus", "Chondracanthus gracilis",
                                                                   "Ergasilus turgidus", "Eudactylina similis",
                                                                   "Pseudocharopinus dentatus", "Haemobaphes diceraus",
                                                                   "Nectobrachia indivisa") ~ "Holoplankton - fish parasites",
      phylum == "Arthropoda" & class== "Copepoda" & family == "Monstrillidae" ~ "Holoplankton - invert parasites",
      phylum == "Arthropoda" & class== "Copepoda" & family == "Caligidae" ~ "Holoplankton - fish parasites",
      phylum == "Arthropoda" & class== "Copepoda" & order %in% c("Calanoida", "Cyclopoida") ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Copepoda" & order == "Harpacticoida" & genus == "Alteutha" ~ "Holoplankton", # this genus has pelagic adults
      phylum == "Arthropoda" & class== "Copepoda" & order == "Harpacticoida" & genus == "Microsetella" ~ "Holoplankton", # this genus has pelagic adults
      phylum == "Arthropoda" & class== "Copepoda" & order == "Harpacticoida" & genus == "Macrosetella"~ "Holoplankton", # this genus has pelagic adults
      phylum == "Arthropoda" & class== "Copepoda" & order == "Harpacticoida" & genus == "Euterpina" ~ "Holoplankton", # this genus has pelagic adults
      phylum == "Arthropoda" & class== "Copepoda" & order == "Harpacticoida" ~ "Benthic", # most of adults are benthic but they have naupliar stages (alghough the nauplii may be benthic...)
      phylum == "Arthropoda" & class== "Ostracoda" ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Branchiopoda" ~ "Holoplankton",
      phylum == "Arthropoda" &  class== "Hexanauplia"   ~ "Holoplankton",
      phylum == "Arthropoda" & class== "Pycnogonida" ~ "Meroplankton",
      #
      # ----Phylum Nemertea----
      #
      phylum == "Nemertea" & species == "Poseidonemertes collaris" ~ "Meroplankton", #https://www.researchgate.net/profile/Svetlana-Maslakova/publication/273513549_From_trochophore_to_pilidium_and_back_again_-_a_larva's_journey/links/57154fb208ae16479d8ac88a/From-trochophore-to-pilidium-and-back-again-a-larvas-journey.pdf
      phylum == "Nemertea" & species %in% c("Paranemertes californica", "Gurjanovella littoralis",
                                            "Emplectonema viride", "Carcinonemertes epialti",
                                            "Procephalotrix spiralis", "Cephalothrix spiralis",
                                            "Amphiporus formidabilis") ~ "Meroplankton", # Mendes 2021 thesis: https://www.teses.usp.br/teses/disponiveis/41/41131/tde-09062022-173135/publico/Cecili_Mendes.pdf#page=25
      phylum == "Nemertea" & genus == "Quasitetrastemma" ~ "Meroplankton", #https://link.springer.com/content/pdf/10.1134/S1063074008040081.pdf 
      phylum == "Nemertea" & genus == "Zygonemertes" ~ "Meroplankton", # https://www.researchgate.net/profile/Svetlana-Maslakova/publication/273513549_From_trochophore_to_pilidium_and_back_again_-_a_larva's_journey/links/57154fb208ae16479d8ac88a/From-trochophore-to-pilidium-and-back-again-a-larvas-journey.pdf
      phylum == "Nemertea" & order == "Heteronemertea" ~ "Meroplankton",  #Shanks's book on Pacific invertebrates
      phylum == "Nemertea" ~ "Holo/Mero",
      #Quasitetrastemma https://link.springer.com/article/10.1134/S1063074008040081
      #
      # ----Other phyla----
      #
      phylum == "Sipuncula" | class == "Sipuncula" ~ "Meroplankton", # depending on taxonomy Sipuncula is as class (e.g., NCBI)
      phylum == "Bryozoa" ~ "Meroplankton",
      phylum == "Ctenophora" ~ "Holoplankton",
      phylum == "Brachiopoda" ~ "Meroplankton",
      phylum == "Chaetognatha" ~ "Holoplankton",
      phylum == "Hemichordata" ~ "Meroplankton",
      phylum == "Echinodermata" ~ "Meroplankton",
      phylum == "Platyhelminthes" ~ "Meroplankton",
      phylum == "Porifera" ~ "Meroplankton",
      phylum == "Entoprocta" ~ "Meroplankton",
      phylum == "Nematoda" ~ "Meroplankton",
      phylum == "Rotifera" ~ "Holoplankton",
      phylum == "Phoronida" ~ "Meroplankton",
      phylum == "Kinorhyncha" ~ "Meroplankton", #https://en.wikipedia.org/wiki/Kinorhyncha
      phylum == "Unidentified" ~ "N/A",
      #
      # ----Add a flag for unclassified taxa
      #
      # Unassigned taxa get a "CHECK" flag, the goal is improving the function over time
      T ~ "CHECK"))
  return(mydf_ed)
}

