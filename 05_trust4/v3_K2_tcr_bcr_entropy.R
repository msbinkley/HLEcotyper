# =============================================================================
# v3_K2_tcr_bcr_entropy.R — TCR and BCR Shannon entropy per sample from TRUST4 report files, run on the Mac
# (author 2026-11-04: "download the bcr tcr data locally and process locally").
#
# Same calculation as the original analysis (computational/trust4/scripts/2024-01-30_BCR.Rmd, which calls getClonality and
# getClonalityTCR in trust4_metric_functions.R):
#   * read the report; keep rows with count > 0
#   * TCR rows: V, J or C starts with TR. BCR heavy rows: V, J or C starts with IGH.
#   * TCR entropy = -sum(p x log2 p), with p = frequency / sum(frequency) over the TCR rows
#   * BCR entropy = the same over heavy-chain clones, after grouping related CDR3s into clusters (below)
# Reading the report: header = FALSE, because the header line starts with "#" and R skips it as a comment. The original
#   Rmd also has a chunk with header = TRUE, which used the first data row as the header and lost it: 33 of the 129 samples
#   of the old table (checked 2026-11-04; all 129 are reproduced exactly once this is accounted for).
# One change: the original aligned each cluster with msa(ClustalW), then removed the gaps again. The alignment does not
#   affect the entropy: only the ungapped sequences are used afterwards. msa is not installed here, so this script skips
#   it. Running this script on the original reports checks that nothing changed (see the Run lines).
#   MergeSeqs, CreateMotifList and MergeMotifs are copied from trust4_metric_functions.R; MergeSeqs no longer handles gaps.
#
# Usage: Rscript v3_K2_tcr_bcr_entropy.R <folder with report files, searched recursively> <output tsv>
# Run:   new reports: data/input/tcr_bcr/2026-11-04_trust4_fastp_reports (fastp-trimmed FASTQs, TRUST4 v1.1.10)
#        old reports: copied from analysis/ajayss/tcellreceptor_analysis/cHL (raw FASTQs, TRUST4 v1.0.8)
# Output columns: sample, TCR_entropy, BCR_entropy, n_TCR_rows, n_BCR_heavy_rows
# No random numbers are used.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

# Step 1. Arguments and report files
arguments <- commandArgs(trailingOnly = TRUE)
report_folder <- arguments[1]
output_file <- arguments[2]
report_files <- list.files(report_folder, pattern = "report\\.tsv$", recursive = TRUE, full.names = TRUE)
cat("report files:", length(report_files), "\n")

# Step 2. Sample name from a file name. New: TRUST_cHL01_report.tsv. Old: cHL01.trust4.out_report.tsv or
#         HL001-T1.trust4.out_report.tsv; the "-T1" part is dropped so old and new samples have the same name.
sample_from_file <- function(file_path) {
  file_name <- basename(file_path)
  file_name <- sub("^TRUST_", "", file_name)
  file_name <- sub("_report\\.tsv$", "", file_name)
  file_name <- sub("\\.trust4\\.out$", "", file_name)
  file_name <- sub("-.*$", "", file_name)
  return(file_name)
}

# Step 3. Helper functions copied from trust4_metric_functions.R
MergeSeqs <- function(seqs, dd) {
  tmp.vv <- order(nchar(seqs))
  seqs <- seqs[tmp.vv]
  dd <- dd[tmp.vv, ]
  newseqs <- seqs
  seq.labs <- rep(0, length(newseqs))
  nn <- length(seq.labs)
  while (TRUE) {
    for (ii in 1:nn) {
      if (length(grep(newseqs[ii], newseqs[-ii])) > 0) {
        seq.labs[ii] <- 1
      }
    }
    if (sum(seq.labs) == 0) {
      break
    }
    vv0 <- which(seq.labs == 0)
    newseqs <- newseqs[vv0]
    seq.labs <- seq.labs[vv0]
    nn <- length(newseqs)
    seqs <- seqs[vv0]
    dd <- dd[vv0, ]
  }
  return(list(SS = seqs, DD = dd))
}

CreateMotifList <- function(mm) {
  tmp.mat <- matrix(unlist(strsplit(rep(mm, 8, ""), "")), nrow = 8, byrow = TRUE)
  diag(tmp.mat) <- "."
  mm.list <- apply(tmp.mat, 1, paste, collapse = "")
  return(mm.list)
}

MergeMotifs <- function(motif.list) {
  # merge motifs by allowing one mismatch
  unique.motifs <- c()
  for (mm in motif.list) {
    mm.list <- CreateMotifList(mm)
    sign <- 0
    for (tmp.mm in mm.list) {
      if (length(grep(tmp.mm, unique.motifs)) > 0) {
        sign <- 1
        break
      }
    }
    if (sign == 0) {
      unique.motifs <- c(unique.motifs, mm)
    }
  }
  return(unique.motifs)
}

# Step 4. Shannon entropy (log2) of a set of clone frequencies, normalised to sum to 1
shannon_entropy <- function(frequencies) {
  proportions <- frequencies / sum(frequencies)
  return(-sum(proportions * log2(proportions)))
}

# Step 5. BCR entropy of one sample, the logic of getClonality(start = 3, end = 10). Positions 3-10 of each complete CDR3
#         give an 8-amino-acid motif; motifs one mismatch apart are merged; CDR3s matching a motif form a cluster, kept if
#         more than 3 sequences remain after removing sequences contained in a longer one. A kept cluster counts as one clone.
MOTIF_START <- 3
MOTIF_END <- 10

bcr_heavy_entropy <- function(heavy_rows) {
  unique_rows <- heavy_rows[!duplicated(heavy_rows$CDR3nt), ]
  complete_cdr3 <- unique_rows$CDR3aa[unique_rows$is_complete == "Y"]
  motif_counts <- table(substr(complete_cdr3, MOTIF_START, MOTIF_END))
  motif_counts <- sort(motif_counts, decreasing = TRUE)
  motif_counts <- motif_counts[which(nchar(names(motif_counts)) == (MOTIF_END - MOTIF_START + 1))]
  motifs <- MergeMotifs(names(motif_counts))

  clustered_nt <- c()
  cluster_frequencies <- c()
  for (motif in motifs) {
    # rows whose CDR3 matches this motif with up to one mismatch
    matching_rows <- c()
    for (motif_pattern in CreateMotifList(motif)) {
      matching_rows <- c(matching_rows, grep(motif_pattern, unique_rows$CDR3aa))
    }
    matching_rows <- unique(matching_rows)
    if (length(matching_rows) < 2) {
      next
    }
    cluster_nt <- unique(as.character(unique_rows[matching_rows, "CDR3nt"]))
    cluster_rows <- unique_rows[matching_rows, ]
    rownames(cluster_rows) <- cluster_rows$CDR3nt
    cluster_rows <- cluster_rows[cluster_nt, ]
    if (length(cluster_nt) < 3) {
      next
    }
    merged <- MergeSeqs(cluster_nt, cluster_rows)
    if (is.null(dim(merged$DD))) {
      next
    }
    if (nrow(merged$DD) <= 3) {
      next
    }
    clustered_nt <- c(clustered_nt, cluster_nt)
    cluster_frequencies <- c(cluster_frequencies, sum(merged$DD$frequency))
  }

  # clones = one per kept cluster + every row not in a kept cluster
  other_frequencies <- heavy_rows$frequency[!(heavy_rows$CDR3nt %in% clustered_nt)]
  return(shannon_entropy(c(cluster_frequencies, other_frequencies)))
}

# Step 6. Process every report
result_rows <- list()
for (report_file in report_files) {
  sample_name <- sample_from_file(report_file)
  if (length(readLines(report_file)) <= 1) {
    cat("empty report, skipped:", sample_name, "\n")
    next
  }
  report <- read.table(file = report_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(report) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "C", "cid", "cid_full_length")
  report <- report[report$count > 0, ]
  report$is_complete <- ifelse(report$CDR3aa != "partial" & report$CDR3aa != "out_of_frame" &
                                 !grepl("^_", report$CDR3aa) & !grepl("^\\?", report$CDR3aa), "Y", "N")

  tcr_rows <- report[grepl("^TR", report$V) | grepl("^TR", report$J) | grepl("^TR", report$C), ]
  heavy_rows <- report[grepl("^IGH", report$V) | grepl("^IGH", report$J) | grepl("^IGH", report$C), ]

  result_rows[[length(result_rows) + 1]] <- data.frame(
    sample = sample_name,
    TCR_entropy = shannon_entropy(tcr_rows$frequency),
    BCR_entropy = bcr_heavy_entropy(heavy_rows),
    n_TCR_rows = nrow(tcr_rows),
    n_BCR_heavy_rows = nrow(heavy_rows))
  cat(length(result_rows), sample_name, "\n")
}

# Step 7. Write the table
results <- do.call(rbind, result_rows)
if (anyDuplicated(results$sample) > 0) {
  stop("duplicated sample names: ", paste(results$sample[duplicated(results$sample)], collapse = ", "))
}
write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written:", output_file, "\n")
