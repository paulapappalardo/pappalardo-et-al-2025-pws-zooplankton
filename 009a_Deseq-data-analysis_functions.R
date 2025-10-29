# ------------------------------------------------------------------------------
# -------------- Phyloseq data analysis and plotting functions -----------------------
# ------------------------------------------------------------------------------

# AUTHOR: Paula Pappalardo - pappalardop@si.edu
# handy Paula addition from phyloseqCompanion similar function
taxa.data.frame <- function(ps) {
  return(as.data.frame(phyloseq::tax_table(ps), "data.frame"))
}

# I made a function to quickly see the significant results and the taxa
getSignficantOnly <- function(res, phyloseq){
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  if(isEmpty(sigtab)){
    print("there are no taxa with differential abundance")
  } else {
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq)[rownames(sigtab), ], "matrix"))
    return(sigtab)
  }
}


