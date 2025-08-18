# Purpose: Extract the effect variable ("x" or "y") from TCEP pair description files,
#          handling non-UTF-8 encodings and heterogeneous arrow notations.
#
# Usage:
#   source("tcep_extract_effect_fixed_utf8.R")
#
# Output:
#   A CSV "ground_truth_effects.csv" written to the same directory as the pair files,
#   with two columns: pair, effect (values: "x", "y", or NA).

pairs_dir <- "xxxxxx/pairs"
out_csv   <- file.path(pairs_dir, "ground_truth_effects.csv")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------- Robust whitespace trimming (encoding-safe) --------------------
trimws_safe <- function(x) {
    if (!length(x)) return(x)
    x <- as.character(x)
    x <- gsub("^[[:space:]]+", "", x, perl = TRUE)
    x <- gsub("[[:space:]]+$", "", x, perl = TRUE)
    x
}

# -------------------- Robust line reader with multi-encoding fallback ----------------
safe_read_lines <- function(path) {
    # Try common encodings including UTF-16 variants; convert everything to UTF-8.
    enc_try <- c("UTF-8", "CP932", "Shift_JIS", "Windows-1252", "latin1",
                 "GB18030", "GBK", "UTF-16LE", "UTF-16BE", "")
    for (enc in enc_try) {
        res <- try({
            con <- file(path, open = "r", encoding = enc)
            on.exit(try(close(con), silent = TRUE), add = TRUE)
            ln <- readLines(con, warn = FALSE)
            ln_u8 <- iconv(ln, from = enc %||% "", to = "UTF-8", sub = "byte")
            ln_u8[is.na(ln_u8)] <- ""
            ln_u8
        }, silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
    }
    # Raw bytes fallback
    raw <- readBin(path, what = "raw", n = file.size(path))
    txt <- rawToChar(raw, multiple = TRUE)
    iconv(txt, from = "", to = "UTF-8", sub = "byte")
}

# -------------------- Parse one line for causal arrow --------------------
# Supports:
#   - Right arrows: "x -> y", "x --> y", "x - - > y", "x → y", "x ⟶ y", "x ⇒ y"  -> effect = right side
#   - Left arrows : "x <- y", "x < - - y", "x ← y", "x ⟵ y", "x ⇐ y"             -> effect = left side
parse_effect_line <- function(ln) {
    s <- tolower(ln)
    
    # Left arrows (effect = left side)
    re_left_ascii   <- "\\b([xy])\\b\\s*<\\s*(?:-\\s*)+\\s*\\b([xy])\\b"
    re_left_unicode <- "\\b([xy])\\b\\s*[←⟵⇐]\\s*\\b([xy])\\b"
    
    m <- regexec(re_left_ascii, s, perl = TRUE)
    r <- regmatches(s, m)[[1]]
    if (length(r) == 3) return(r[2])
    
    m <- regexec(re_left_unicode, s, perl = TRUE)
    r <- regmatches(s, m)[[1]]
    if (length(r) == 3) return(r[2])
    
    # Right arrows (effect = right side)
    re_right_ascii   <- "\\b([xy])\\b\\s*(?:-\\s*)+>\\s*\\b([xy])\\b"
    re_right_unicode <- "\\b([xy])\\b\\s*[→⟶⇒]\\s*\\b([xy])\\b"
    
    m <- regexec(re_right_ascii, s, perl = TRUE)
    r <- regmatches(s, m)[[1]]
    if (length(r) == 3) return(r[3])
    
    m <- regexec(re_right_unicode, s, perl = TRUE)
    r <- regmatches(s, m)[[1]]
    if (length(r) == 3) return(r[3])
    
    return(NA_character_)
}

# -------------------- Extract effect from an entire description file -----------------
extract_effect_from_lines <- function(lines) {
    if (length(lines) == 0) return(NA_character_)
    
    L    <- trimws_safe(lines)
    L    <- L[nzchar(L)]
    L_lc <- tolower(L)
    
    # 1) Prefer the "ground truth" block (case-insensitive; colon optional).
    gt_hits <- grep("\\bground\\s*truth\\b\\s*:?", L_lc, perl = TRUE)
    if (length(gt_hits) >= 1) {
        for (idx in gt_hits) {
            # Same line first (e.g., "Ground truth: x -> y")
            eff <- parse_effect_line(L[idx])
            if (!is.na(eff)) return(eff)
            
            # Then the following lines (up to 15 or until a blank line/paragraph break).
            start <- idx + 1
            end   <- min(length(L), start + 14)
            block <- c()
            for (i in start:end) {
                if (!nzchar(L[i])) break
                block <- c(block, L[i])
            }
            for (ln in block) {
                eff <- parse_effect_line(ln)
                if (!is.na(eff)) return(eff)
            }
        }
    }
    
    # 2) Fallback: scan the whole file, keep the last arrow occurrence to avoid noise.
    last_eff <- NA_character_
    for (i in seq_along(L)) {
        eff <- parse_effect_line(L[i])
        if (!is.na(eff)) last_eff <- eff
    }
    if (!is.na(last_eff)) return(last_eff)
    
    # 3) Semantic fallback: handle verbal statements like "x causes y",
    #    using labels extracted from "x = ..." and "y = ..." lines.
    get_label <- function(pattern) {
        m <- regexpr(pattern, L_lc, perl = TRUE)
        if (any(m > 0)) {
            i <- which(m > 0)[1]
            txt <- sub(pattern, "\\1", L_lc[i], perl = TRUE)
            txt <- gsub("[^a-z0-9\\s_/\\-]+", " ", txt) # simplify punctuation
            return(trimws_safe(txt))
        }
        NA_character_
    }
    x_lab <- get_label("^\\s*x\\s*=\\s*(.+)$")
    y_lab <- get_label("^\\s*y\\s*=\\s*(.+)$")
    
    if (!is.na(x_lab) && !is.na(y_lab)) {
        q <- function(s) gsub("([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1", s, perl = TRUE)
        flex <- function(s) gsub("\\s+", "\\\\s+", q(s))  # allow whitespace variation
        
        re_xy <- paste0("\\b", flex(x_lab), "\\b.*\\b(cause|causes|caused|determine|determines|influence|influences)\\b.*\\b", flex(y_lab), "\\b")
        re_yx <- paste0("\\b", flex(y_lab), "\\b.*\\b(cause|causes|caused|determine|determines|influence|influences)\\b.*\\b", flex(x_lab), "\\b")
        
        if (any(grepl(re_xy, L_lc, perl = TRUE))) return("y") # x causes y -> effect = y
        if (any(grepl(re_yx, L_lc, perl = TRUE))) return("x") # y causes x -> effect = x
    }
    
    # 4) No parse
    return(NA_character_)
}

# -------------------- Batch over all *_des.txt files --------------------
des_files <- list.files(pairs_dir, pattern = "^pair\\d{4}_des\\.txt$", full.names = TRUE)
if (length(des_files) == 0) stop("No *_des.txt files found in directory: ", pairs_dir)

pair_idx <- as.integer(sub("^pair(\\d{4})_des\\.txt$", "\\1", basename(des_files)))
ord <- order(pair_idx)
des_files <- des_files[ord]
pair_idx  <- pair_idx[ord]

effects  <- character(length(des_files))
warn_cnt <- 0L

for (k in seq_along(des_files)) {
    f <- des_files[k]
    lines <- safe_read_lines(f)
    eff <- extract_effect_from_lines(lines)
    if (is.na(eff)) {
        warn_cnt <- warn_cnt + 1L
        message(sprintf("[WARN] Could not parse %s -> set to NA", basename(f)))
    }
    effects[k] <- eff %||% NA_character_
}

# -------------------- Manual curation for known problematic pairs --------------------
# Academic note for record-keeping:
# After running the script, pairs 72, 86, 88, 97, 98, 106, and 108 cannot be automatically resolved.
# The underlying cause is that the corresponding description files contain either character-encoding
# anomalies or non-standard, heterogeneous causal annotations (e.g., spaced arrows "x - - > y",
# left arrows "x <- y", purely verbal statements such as "x causes y", or no explicit ground-truth
# arrow line, often interleaved with UPDATE notes), which prevents rule-based parsing from reliably
# extracting the effect variable. Accordingly, we set the effect variable to y for pairs 72, 86, 88,
# 97, and 98, and to x for pairs 106 and 108.

manual_effects <- c(`72`="y", `86`="y", `88`="y", `97`="y", `98`="y", `106`="x", `108`="x")

for (p in names(manual_effects)) {
    pos <- which(pair_idx == as.integer(p))
    if (length(pos)) {
        needs_fill <- is.na(effects[pos]) | !nzchar(effects[pos])
        if (any(needs_fill)) {
            effects[pos[needs_fill]] <- manual_effects[[p]]
            message(sprintf("[INFO] Manually set effect for pair%04d -> %s",
                            as.integer(p), manual_effects[[p]]))
        }
    }
}

# -------------------- Write results --------------------
res <- data.frame(
    pair   = sprintf("pair%04d", pair_idx),
    effect = effects,
    stringsAsFactors = FALSE
)

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)

cat(sprintf("Done: wrote %s\n", out_csv))
cat(sprintf("Summary: %d files processed, %d could not be parsed before manual curation (effect=NA).\n",
            nrow(res), warn_cnt))
