suppressPackageStartupMessages({
  library(Biostrings)
})

# Working directory: C:/Users/baysa/Documents/Weedgenes
setwd("C:/Users/baysa/Documents/Weedgenes")

cds_file <- "JG_gene_highquality_set_aed03_cleaned.CDS"
stopifnot(file.exists(cds_file))

cds_all <- readDNAStringSet(cds_file)

get_kmers_single <- function(dna_seq, k = 21L) {
  n <- length(dna_seq)
  if (n < k) return(character(0))
  as.character(DNAStringSet(Views(dna_seq, start = seq_len(n - k + 1L), width = k)))
}

k <- 21L
valid_cds <- cds_all[width(cds_all) >= k]
stopifnot(length(valid_cds) > 0L)

# Count in how many CDS entries each unique 21-mer appears
per_seq_unique_kmers <- lapply(valid_cds, function(s) unique(get_kmers_single(s, k = k)))
kmer_doc_freq <- table(unlist(per_seq_unique_kmers, use.names = FALSE))

# Score each CDS sequence (subtract 1 to exclude self-occurrence)
cds_scores_list <- lapply(seq_along(valid_cds), function(i) {
  seq_i <- valid_cds[[i]]
  all_km <- get_kmers_single(seq_i, k = k)
  total_km <- length(all_km)
  if (total_km == 0L) return(NULL)

  other_counts <- as.integer(kmer_doc_freq[all_km]) - 1L
  shared_km <- sum(other_counts > 0L)
  unique_km <- total_km - shared_km

  data.frame(
    gene_id = names(valid_cds)[i],
    cds_length = length(seq_i),
    total_21mers = total_km,
    unique_21mers = unique_km,
    shared_21mers = shared_km,
    unique_fraction = round(unique_km / total_km, 4),
    stringsAsFactors = FALSE
  )
})

cds_scores <- do.call(rbind, cds_scores_list)
clean_pool <- subset(cds_scores, cds_length >= 250L)
clean_pool <- clean_pool[order(-clean_pool$unique_fraction,
                              clean_pool$shared_21mers,
                              -clean_pool$cds_length), ]
rownames(clean_pool) <- NULL

write.csv(clean_pool, "cds_21mer_uniqueness_ranking.csv", row.names = FALSE)
cat("Saved cds_21mer_uniqueness_ranking.csv with", nrow(clean_pool), "rows.\n")
print(head(clean_pool, 10))
                               
suppressPackageStartupMessages({
  library(Biostrings)
})

setwd("C:/Users/baysa/Documents/Weedgenes")

cds_file <- "JG_gene_highquality_set_aed03_cleaned.CDS"
stopifnot(file.exists(cds_file))
cds_all <- readDNAStringSet(cds_file)

# Use top 5 cleaner genes from cds_21mer_uniqueness_ranking.csv if available;
# otherwise fall back to specified target IDs.
if (file.exists("cds_21mer_uniqueness_ranking.csv")) {
  ranking_df <- read.csv("cds_21mer_uniqueness_ranking.csv", stringsAsFactors = FALSE)
  target_ids <- head(ranking_df$gene_id, 5)
} else {
  target_ids <- c(
    "maker-h1tg000001l-exonerate_protein2genome-gene-1.0-mRNA-1",
    "maker-h1tg000001l-exonerate_protein2genome-gene-2.0-mRNA-1"
  )
}

target_idx <- which(names(cds_all) %in% target_ids)
stopifnot(length(target_idx) > 0L)

target_cds <- cds_all[target_idx]
background_cds <- cds_all[-target_idx]

extract_bg_kmers <- function(dna_set, k = 21L) {
  valid_seqs <- dna_set[width(dna_set) >= k]
  if (length(valid_seqs) == 0L) return(character(0))
  kmers <- unlist(lapply(valid_seqs, function(s) {
    as.character(DNAStringSet(Views(s, start = seq_len(length(s) - k + 1L), width = k)))
  }), use.names = FALSE)
  unique(kmers)
}

bg_21mers <- extract_bg_kmers(background_cds, k = 21L)

score_windows <- function(seq_dna, gene_id, bg_kmers,
                          win_size = 250L, step_size = 25L,
                          k = 21L, gc_min = 0.35, gc_max = 0.55) {
  n <- length(seq_dna)
  if (n < win_size) return(NULL)

  starts <- seq(1L, n - win_size + 1L, by = step_size)
  rows <- lapply(starts, function(st) {
    en <- st + win_size - 1L
    sub_seq <- subseq(seq_dna, st, en)
    if (letterFrequency(sub_seq, "N", as.prob = FALSE) > 0L) return(NULL)

    gc <- sum(letterFrequency(sub_seq, c("G", "C"), as.prob = TRUE))
    w_kmers <- as.character(DNAStringSet(Views(sub_seq, start = seq_len(win_size - k + 1L), width = k)))
    hits <- sum(w_kmers %in% bg_kmers)

    data.frame(
      gene_id = gene_id,
      start = st,
      end = en,
      length = win_size,
      gc_content = round(gc, 4),
      gc_in_range = (gc >= gc_min) & (gc <= gc_max),
      bg_21mer_hits = hits,
      sequence = as.character(sub_seq),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

all_candidates <- do.call(rbind, mapply(
  FUN = score_windows,
  seq_dna = target_cds,
  gene_id = names(target_cds),
  MoreArgs = list(bg_kmers = bg_21mers),
  SIMPLIFY = FALSE
))

stopifnot(!is.null(all_candidates), nrow(all_candidates) > 0L)

# Select up to 2 windows per gene: strict zero-hit + target GC if present, else lowest hits + GC near 0.45
top_candidates <- do.call(rbind, lapply(split(all_candidates, all_candidates$gene_id), function(df) {
  strict_df <- subset(df, bg_21mer_hits == 0L & gc_in_range == TRUE)
  pool <- if (nrow(strict_df) > 0L) strict_df else df
  pool <- pool[order(pool$bg_21mer_hits, !pool$gc_in_range, abs(pool$gc_content - 0.45), pool$start), ]
  head(pool, 2)
}))
rownames(top_candidates) <- NULL

write.csv(top_candidates, "candidate_dsrna_summary.csv", row.names = FALSE)

cand_seqs <- DNAStringSet(top_candidates$sequence)
names(cand_seqs) <- paste0(
  top_candidates$gene_id, "_",
  top_candidates$start, "-", top_candidates$end,
  "_GC", top_candidates$gc_content,
  "_hits", top_candidates$bg_21mer_hits
)
writeXStringSet(cand_seqs, "candidate_dsrna_sequences.fasta")

# Export top 1 candidate per gene
top1_df <- do.call(rbind, lapply(split(top_candidates, top_candidates$gene_id), function(df) head(df, 1)))
rownames(top1_df) <- NULL
top1_seqs <- DNAStringSet(top1_df$sequence)
names(top1_seqs) <- paste0(
  top1_df$gene_id, "_",
  top1_df$start, "-", top1_df$end,
  "_GC", top1_df$gc_content,
  "_hits", top1_df$bg_21mer_hits
)
writeXStringSet(top1_seqs, "selected_top1_per_gene.fasta")

cat("Exported candidate_dsrna_summary.csv, candidate_dsrna_sequences.fasta, and selected_top1_per_gene.fasta\n")
print(top_candidates[, c("gene_id", "start", "end", "gc_content", "bg_21mer_hits")])
suppressPackageStartupMessages({
  library(Biostrings)
})

setwd("C:/Users/baysa/Documents/Weedgenes")

cds_file <- "JG_gene_highquality_set_aed03_cleaned.CDS"
stopifnot(file.exists(cds_file))
cds_all <- readDNAStringSet(cds_file)

# Use top 5 cleaner genes from cds_21mer_uniqueness_ranking.csv if available;
# otherwise fall back to specified target IDs.
if (file.exists("cds_21mer_uniqueness_ranking.csv")) {
  ranking_df <- read.csv("cds_21mer_uniqueness_ranking.csv", stringsAsFactors = FALSE)
  target_ids <- head(ranking_df$gene_id, 5)
} else {
  target_ids <- c(
    "maker-h1tg000001l-exonerate_protein2genome-gene-1.0-mRNA-1",
    "maker-h1tg000001l-exonerate_protein2genome-gene-2.0-mRNA-1"
  )
}

target_idx <- which(names(cds_all) %in% target_ids)
stopifnot(length(target_idx) > 0L)

target_cds <- cds_all[target_idx]
background_cds <- cds_all[-target_idx]

extract_bg_kmers <- function(dna_set, k = 21L) {
  valid_seqs <- dna_set[width(dna_set) >= k]
  if (length(valid_seqs) == 0L) return(character(0))
  kmers <- unlist(lapply(valid_seqs, function(s) {
    as.character(DNAStringSet(Views(s, start = seq_len(length(s) - k + 1L), width = k)))
  }), use.names = FALSE)
  unique(kmers)
}

bg_21mers <- extract_bg_kmers(background_cds, k = 21L)

score_windows <- function(seq_dna, gene_id, bg_kmers,
                          win_size = 250L, step_size = 25L,
                          k = 21L, gc_min = 0.35, gc_max = 0.55) {
  n <- length(seq_dna)
  if (n < win_size) return(NULL)

  starts <- seq(1L, n - win_size + 1L, by = step_size)
  rows <- lapply(starts, function(st) {
    en <- st + win_size - 1L
    sub_seq <- subseq(seq_dna, st, en)
    if (letterFrequency(sub_seq, "N", as.prob = FALSE) > 0L) return(NULL)

    gc <- sum(letterFrequency(sub_seq, c("G", "C"), as.prob = TRUE))
    w_kmers <- as.character(DNAStringSet(Views(sub_seq, start = seq_len(win_size - k + 1L), width = k)))
    hits <- sum(w_kmers %in% bg_kmers)

    data.frame(
      gene_id = gene_id,
      start = st,
      end = en,
      length = win_size,
      gc_content = round(gc, 4),
      gc_in_range = (gc >= gc_min) & (gc <= gc_max),
      bg_21mer_hits = hits,
      sequence = as.character(sub_seq),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

all_candidates <- do.call(rbind, mapply(
  FUN = score_windows,
  seq_dna = target_cds,
  gene_id = names(target_cds),
  MoreArgs = list(bg_kmers = bg_21mers),
  SIMPLIFY = FALSE
))

stopifnot(!is.null(all_candidates), nrow(all_candidates) > 0L)

# Select up to 2 windows per gene: strict zero-hit + target GC if present, else lowest hits + GC near 0.45
top_candidates <- do.call(rbind, lapply(split(all_candidates, all_candidates$gene_id), function(df) {
  strict_df <- subset(df, bg_21mer_hits == 0L & gc_in_range == TRUE)
  pool <- if (nrow(strict_df) > 0L) strict_df else df
  pool <- pool[order(pool$bg_21mer_hits, !pool$gc_in_range, abs(pool$gc_content - 0.45), pool$start), ]
  head(pool, 2)
}))
rownames(top_candidates) <- NULL

write.csv(top_candidates, "candidate_dsrna_summary.csv", row.names = FALSE)

cand_seqs <- DNAStringSet(top_candidates$sequence)
names(cand_seqs) <- paste0(
  top_candidates$gene_id, "_",
  top_candidates$start, "-", top_candidates$end,
  "_GC", top_candidates$gc_content,
  "_hits", top_candidates$bg_21mer_hits
)
writeXStringSet(cand_seqs, "candidate_dsrna_sequences.fasta")

# Export top 1 candidate per gene
top1_df <- do.call(rbind, lapply(split(top_candidates, top_candidates$gene_id), function(df) head(df, 1)))
rownames(top1_df) <- NULL
top1_seqs <- DNAStringSet(top1_df$sequence)
names(top1_seqs) <- paste0(
  top1_df$gene_id, "_",
  top1_df$start, "-", top1_df$end,
  "_GC", top1_df$gc_content,
  "_hits", top1_df$bg_21mer_hits
)
writeXStringSet(top1_seqs, "selected_top1_per_gene.fasta")

cat("Exported candidate_dsrna_summary.csv, candidate_dsrna_sequences.fasta, and selected_top1_per_gene.fasta\n")
print(top_candidates[, c("gene_id", "start", "end", "gc_content", "bg_21mer_hits")])

suppressPackageStartupMessages({
  library(Biostrings)
})

setwd("C:/Users/baysa/Documents/Weedgenes")

cds_file <- "JG_gene_highquality_set_aed03_cleaned.CDS"
summary_file <- "candidate_dsrna_summary.csv"
stopifnot(file.exists(cds_file), file.exists(summary_file))

cds_all <- readDNAStringSet(cds_file)
cand_df <- read.csv(summary_file, stringsAsFactors = FALSE)

# Use top 1 window per gene from candidate_dsrna_summary.csv
targets <- do.call(rbind, lapply(split(cand_df, cand_df$gene_id), function(df) {
  df <- df[order(df$bg_21mer_hits, abs(df$gc_content - 0.45), df$start), ]
  df[1, c("gene_id", "start", "end"), drop = FALSE]
}))
rownames(targets) <- NULL

design_flanking_primers <- function(gene_id, ds_start, ds_end, cds_set,
                                    p_len = 20L, amp_min = 70L, amp_max = 180L,
                                    min_gap = 10L, gc_range = c(0.35, 0.70),
                                    tm_range = c(52, 68)) {
  if (!(gene_id %in% names(cds_set))) return(NULL)

  # cds_set[[gene_id]] is a DNAString -> use length(), not width()
  gene_seq <- cds_set[[gene_id]]
  cds_len <- length(gene_seq)
  pairs <- list()

  scan_side <- function(reg_seq, reg_offset, label, side) {
    n <- length(reg_seq)
    if (n < (p_len * 2L + 10L)) return()

    starts <- seq_len(n - p_len + 1L)
    cands_list <- lapply(starts, function(st) {
      en <- st + p_len - 1L
      sub_f <- subseq(reg_seq, st, en)
      if (letterFrequency(sub_f, "N", as.prob = FALSE) > 0L) return(NULL)

      p_f <- as.character(sub_f)
      p_r <- as.character(reverseComplement(sub_f))
      gc <- sum(letterFrequency(sub_f, c("G", "C"), as.prob = TRUE))
      tm <- round(4 * sum(letterFrequency(sub_f, c("G", "C"))) +
                    2 * sum(letterFrequency(sub_f, c("A", "T"))), 1)

      if (gc >= gc_range[1] && gc <= gc_range[2] && tm >= tm_range[1] && tm <= tm_range[2]) {
        data.frame(
          abs_start = reg_offset + st - 1L,
          abs_end = reg_offset + en - 1L,
          tm = tm,
          gc = round(gc, 3),
          fwd_seq = p_f,
          rev_seq = p_r,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    })

    cdf <- do.call(rbind, cands_list)
    if (is.null(cdf) || nrow(cdf) < 2L) return()

    for (i in seq_len(nrow(cdf))) {
      for (j in seq_len(nrow(cdf))) {
        if (i == j) next
        amp_len <- cdf$abs_end[j] - cdf$abs_start[i] + 1L
        outside_ok <- if (side == "up") {
          cdf$abs_end[i] < cdf$abs_start[j] && cdf$abs_end[j] < ds_start
        } else {
          cdf$abs_end[i] < cdf$abs_start[j] && cdf$abs_start[i] > ds_end
        }
        if (outside_ok && amp_len >= amp_min && amp_len <= amp_max) {
          pairs[[length(pairs) + 1L]] <<- data.frame(
            gene_id = gene_id,
            location = label,
            ds_window = paste0(ds_start, "-", ds_end),
            fwd_pos = paste0(cdf$abs_start[i], "-", cdf$abs_end[i]),
            rev_pos = paste0(cdf$abs_start[j], "-", cdf$abs_end[j]),
            amplicon_bp = amp_len,
            tm_diff = round(abs(cdf$tm[i] - cdf$tm[j]), 2),
            fwd_seq = cdf$fwd_seq[i],
            rev_seq = cdf$rev_seq[j],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  up_end <- ds_start - min_gap - 1L
  dn_start <- ds_end + min_gap + 1L

  if (up_end >= 35L) {
    scan_side(subseq(gene_seq, 1L, up_end), 1L, "upstream_of_dsRNA", "up")
  }
  if (dn_start <= cds_len - 35L) {
    scan_side(subseq(gene_seq, dn_start, cds_len), dn_start, "downstream_of_dsRNA", "dn")
  }

  if (length(pairs) == 0L) return(NULL)
  out <- do.call(rbind, pairs)
  out <- out[order(out$tm_diff, abs(out$amplicon_bp - 110)), ]
  head(out, 3)
}

qpcr_list <- lapply(seq_len(nrow(targets)), function(r) {
  design_flanking_primers(
    gene_id = targets$gene_id[r],
    ds_start = targets$start[r],
    ds_end = targets$end[r],
    cds_set = cds_all
  )
})

qpcr_candidates <- do.call(rbind, qpcr_list)

if (!is.null(qpcr_candidates) && nrow(qpcr_candidates) > 0L) {
  rownames(qpcr_candidates) <- NULL
  write.csv(qpcr_candidates, "rt_qpcr_primer_candidates.csv", row.names = FALSE)
  cat("Saved rt_qpcr_primer_candidates.csv with", nrow(qpcr_candidates), "primer pairs.\n")
  print(qpcr_candidates)
} else {
  cat("No flanking primer pairs met current thresholds.\n")
}

