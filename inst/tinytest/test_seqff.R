library(tinytest)
library(RWisecondorX)

support <- RWisecondorX:::seqff_load_support_data()
counts <- data.frame(
  binName = support$bininfo$binName,
  counts = as.integer((support$bininfo$GC * 1000) + seq_len(nrow(support$bininfo)) %% 17L + 50L),
  stringsAsFactors = FALSE
)

ff <- seqff_predict(counts, input_type = "counts")

expect_equal(names(ff), c("SeqFF", "Enet", "WRSC"))
expect_true(all(is.finite(ff)))
