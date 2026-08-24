library(stringr)
# ================================================================
# Function to extract unique sorted integer positions from vir_df columns
# ================================================================
extract_positions <- function(vec) {
  # vec: character vector (possibly length >1), may contain NA or '-'
  vec <- vec[!is.na(vec)]
  if (length(vec) == 0) return(integer(0))

  parts <- unlist(str_split(paste(vec, collapse = ","), ","))
  parts <- str_trim(parts)
  parts <- parts[parts != "" & parts != "-"]

  # convert safely to integer, drop non-numeric entries
  nums <- suppressWarnings(as.integer(parts))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0) return(integer(0))

  sort(unique(nums))
}


# ================================================================
# Cliff's delta (non-parametric effect size)
# - Accepts a grouped tibble/data.frame piped in and a formula: value ~ group
# - Returns a tibble with columns: group1, group2, n1, n2, delta, magnitude
# - Usage: df %>% group_by(Gene) %>% cliff_delta(KL ~ Interface)
# ================================================================
cliff_delta <- function(.data, formula) {
  # lazy load dplyr to avoid hard dependency at package load
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for cliff_delta()")
  }
  # allow passing formula or character
  if (missing(formula)) stop("formula is required, e.g. value ~ group")
  fvars <- all.vars(formula)
  if (length(fvars) < 2) stop("formula must be of the form value ~ group")
  value_name <- fvars[1]
  group_name <- fvars[2]

  # per-group computation using dplyr::group_modify
  dplyr::group_modify(.data, function(df, key) {
    # defensive: ensure columns exist
    if (!(value_name %in% names(df)) || !(group_name %in% names(df))) {
      return(dplyr::tibble(group1 = NA_character_, group2 = NA_character_,
                           n1 = NA_integer_, n2 = NA_integer_,
                           delta = NA_real_, magnitude = NA_character_))
    }

    vals <- df[[value_name]]
    gr <- as.character(df[[group_name]])
    # remove NA values
    keep <- !is.na(vals) & !is.na(gr)
    vals <- vals[keep]
    g <- df[[group_name]][keep]

    levs <- levels(droplevels(g))
    if (length(levs) != 2) {
      # return NA row if not exactly two groups
      return(dplyr::tibble(group1 = NA_character_, group2 = NA_character_,
                           n1 = NA_integer_, n2 = NA_integer_,
                           delta = NA_real_, magnitude = NA_character_))
    }

    g1 <- levs[1]
    g2 <- levs[2]
    x1 <- vals[gr == g1]
    x2 <- vals[gr == g2]
    n1 <- length(x1)
    n2 <- length(x2)
    if (n1 == 0 || n2 == 0) {
      return(dplyr::tibble(group1 = as.character(g1), group2 = as.character(g2),
                           n1 = n1, n2 = n2,
                           delta = NA_real_, magnitude = NA_character_))
    }

    # compute pairwise comparisons: count greater / less
    # use outer with vectorized sign, then sum
    comp <- outer(x1, x2, FUN = function(a, b) sign(a - b))
    n_pos <- sum(comp == 1)
    n_neg <- sum(comp == -1)
    delta <- (n_pos - n_neg) / (n1 * n2)

    absd <- abs(delta)
    magnitude <- if (is.na(absd)) NA_character_
    else if (absd < 0.147) "negligible"
    else if (absd < 0.33) "small"
    else if (absd < 0.474) "medium"
    else "large"

    dplyr::tibble(group1 = as.character(g1), group2 = as.character(g2),
                  n1 = n1, n2 = n2, delta = delta, magnitude = magnitude)
  }) %>% dplyr::ungroup()
}

