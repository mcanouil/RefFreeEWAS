HNSCC <- load("data-raw/HNSCC.Rdata")
save(list = HNSCC, file = "data/HNSCC.rda", compress = "bzip2", version = 2, ascii = FALSE)

RefFreeEWAS <- load("data-raw/RefFreeEWAS.Rdata")
save(list = RefFreeEWAS, file = "data/RefFreeEWAS.rda", compress = "bzip2", version = 2, ascii = FALSE)

writeLines(
  text = c(
    paste("HNSCC:", paste(HNSCC, collapse = " ")),
    paste("RefFreeEWAS:", paste(RefFreeEWAS, collapse = " "))
  ),
  con = "data/datalist"
)
