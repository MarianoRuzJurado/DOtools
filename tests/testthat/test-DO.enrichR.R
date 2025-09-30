library(testthat)
library(DOtools)
library(dplyr)

test_that("DO.enrichR basic functionality with mocked enrichr", {
  # Mock DGE data
  df_DGE <- data.frame(
    gene = paste0("Gene", 1:4),
    avg_log2FC_SC_wilcox = c(1, 0.5, -1, -0.7),
    p_val_SC_wilcox = c(0.01, 0.02, 0.01, 0.03),
    stringsAsFactors = FALSE
  )

  # Dummy enrichr function that always returns a data frame
  dummy_enrichr <- function(genes, databases) {
    res <- lapply(databases, function(db) {
      data.frame(
        Term = paste0("GO_term_", 1:2),
        P.value = c(0.01, 0.02),
        stringsAsFactors = FALSE
      )
    })
    names(res) <- databases
    return(res)
  }

  # Mock enrichr during the test
  with_mock(
    `enrichR::enrichr` = dummy_enrichr,
    `enrichR::setEnrichrSite` = function(site) NULL,
    {
      result <- DO.enrichR(
        df_DGE = df_DGE,
        gene_column = "gene",
        pval_column = "p_val_SC_wilcox",
        log2fc_column = "avg_log2FC_SC_wilcox",
        pval_cutoff = 0.05,
        log2fc_cutoff = 0.25,
        path = NULL,
        filename = "",
        species = "Human",
        go_catgs = c("GO_Biological_Process_2023")
      )

      expect_true(is.data.frame(result))
      expect_true(all(c("Database", "State", "Term", "P.value") %in% colnames(result)))
      expect_true(all(result$State %in% c("enriched", "depleted")))
      expect_true(all(result$Database == "GO_Biological_Process_2023"))
    }
  )
})
