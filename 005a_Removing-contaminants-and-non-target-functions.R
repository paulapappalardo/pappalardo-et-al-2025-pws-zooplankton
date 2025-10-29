# ------------------------------------------------------------------------------
# ------------Removing contaminants and non-target functions -------------------
# ------------------------------------------------------------------------------

# AUTHOR: Paula Pappalardo - pappalardop@si.edu

# handy Paula addition from phyloseqCompanion similar function
taxa.data.frame <- function(ps) {
  return(as(phyloseq::tax_table(ps), "data.frame"))
}

simplifyTaxTable <- function(myphyloseq){
  # Keeps only taxonomic levels on tax table
  #
  # Args:
  # phyloseq object
  # Returns:
  # same phyloseq object but with taxa table including only taxonomic columns
  #
  # extract tax table as dataframe
  updatedtaxtable <- taxa.data.frame(myphyloseq) %>%
    select(kingdom, phylum, class, order, family, genus, species, sciname) 
  # convert taxtable to phyloseq format
  mat <- updatedtaxtable %>% as.matrix(.)
  rownames(mat) <- rownames(updatedtaxtable)
  taxa_updated <- tax_table(mat)
  # update the phyloseq object tax table
  tax_table(myphyloseq) <- taxa_updated
  # return updated phyloseq
  myphyloseq
}
