library(tinytest)
library(RWisecondorX)

donor_bam <- system.file("extdata", "nipter_conformance_fixture.bam",
                         package = "RWisecondorX")
if (!nzchar(donor_bam)) {
  exit_file("Bundled NIPTeR conformance fixture is missing; run `make fixtures`.")
}

tmp_dir <- tempfile("bam_sim_cohort_", tmpdir = getwd())
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(tmp_dir)) {
  stop("Failed to create temporary BAM cohort simulation directory.",
       call. = FALSE)
}

out_dir <- file.path(tmp_dir, "cohort")
manifest <- simulate_trisomy_cohort(
  bams = c(donor_bam, donor_bam),
  donor_ids = c("d1", "d2"),
  out_dir = out_dir,
  trisomy_chromosomes = "21",
  fetal_fraction = 0.20,
  mosaic = c(FALSE, TRUE),
  seed = 100L,
  index = TRUE,
  include_donors = TRUE
)

manifest_path <- attr(manifest, "manifest_path")
expect_true(file.exists(manifest_path),
            info = "simulate_trisomy_cohort writes a manifest file")

manifest_disk <- utils::read.table(manifest_path, header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE)
expect_equal(nrow(manifest), 6L,
             info = "manifest contains donor negatives and four positives")
expect_equal(nrow(manifest_disk), 6L,
             info = "written manifest row count matches returned manifest")

donor_rows <- manifest[manifest$cohort_role == "donor_negative", , drop = FALSE]
sim_rows <- manifest[manifest$cohort_role == "simulated_positive", , drop = FALSE]

expect_equal(nrow(donor_rows), 2L,
             info = "manifest includes both donor negatives")
expect_equal(nrow(sim_rows), 4L,
             info = "manifest includes one simulation per donor and mosaic setting")
expect_true(all(donor_rows$output_bam == donor_rows$source_bam),
            info = "donor negatives point back to the source BAM")

expect_true(all(file.exists(sim_rows$output_bam)),
            info = "all simulated BAMs are written")
expect_true(all(file.exists(sim_rows$index_path)),
            info = "all simulated BAMs are indexed")

for (donor_id in c("d1", "d2")) {
  donor_sim <- sim_rows[sim_rows$donor_id == donor_id, , drop = FALSE]
  expect_equal(nrow(donor_sim), 2L,
               info = sprintf("donor %s has two simulated positives", donor_id))
  non_mosaic <- donor_sim[!donor_sim$mosaic, , drop = FALSE]
  mosaic <- donor_sim[donor_sim$mosaic, , drop = FALSE]
  expect_true(mosaic$output_records > non_mosaic$output_records,
              info = sprintf("mosaic simulation for %s keeps more records",
                             donor_id))
}

unlink(tmp_dir, recursive = TRUE)
