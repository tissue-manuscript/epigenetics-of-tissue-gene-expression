library(RnBeads)

parallel.setup(20)

rns <- load.rnb.set(path="rnbSet_unnormalized.zip")
rnb.set <- rns
report.dir <- "test_output"

rnb.options(
  identifiers.column = "sampleID", 
  import.bed.style = "bismarkCov", 
  assembly = "rn6",
  filtering.sex.chromosomes.removal = TRUE,
  disk.dump.big.matrices = TRUE,
  disk.dump.bigff = TRUE,
  logging.disk = TRUE,
  enforce.memory.management = TRUE
)


rnb.set.unfiltered <- rnb.set

result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports=report.dir)

rnb.set <- result$rnb.set

rnb.run.exploratory(rnb.set, report.dir)

rnb.run.differential(rnb.set, report.dir)
