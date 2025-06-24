.onLoad <- function(libname, pkgname) {
  bioc_pkgs <- c("limma", "ComplexHeatmap", "edgeR")
  
  missing_pkgs <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "))
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    BiocManager::install(missing_pkgs, ask = FALSE)
  }
}
