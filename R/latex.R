#' @importFrom dplyr case_when mutate mutate_all rename_at select
#' @importFrom purrr map_chr map2_chr possibly
#' @importFrom magrittr %>% extract
#' @importFrom xtable sanitize

#' @export
output_latex <- function(
  reg_table,
  group_names = NULL,
  caption = NULL,
  env = 'longtable',
  repeat_heading = TRUE,
  bottom_rule = TRUE
) {
  # What happens if there's not multiple estimates? We don't have tacked on,
  # .a, .b suffixes, but we do require those to do the rest, so we rename
  # to .a suffix even if it is just one estimate. Note we _need_ to use
  # possibly because rename_at errors out if there is nothing to rename
  # which seems a bit counter-intuitive.
  reg_table <- possibly(rename_at, reg_table)(
    reg_table,
    vars(matches('^estimate$|^std\\.error$|^p\\.value$')),
    ~ paste0(.x, '.a')
  )

  total_rows <- nrow(reg_table)
  non_title_cols <- !grepl('(label|flevels|type)', names(reg_table))
  total_cols <- sum(non_title_cols) + 1

  stat_cols <- names(reg_table) %>% extract(grepl('\\.[a-z]$', .))
  model_groupings <- map_chr(stat_cols, ~ substr(.x, nchar(.x), nchar(.x)))
  header_rle <- rle(model_groupings)
  # We need to fix this to take user input if desired.
  if (is.null(group_names)) {
    group_names <- map_chr(seq_along(header_rle$lengths), ~ sprintf('(%s)', .x))
  } else {
    group_names <- group_names %>% extract(1:length(header_rle$lengths))
  }
  group_names_full <- map2_chr(
    group_names,
    header_rle$lengths,
    ~ sprintf('\\multicolumn{%s}{c}{%s}', .y, .x)
  )
  group_names_full <- paste0(
    ' & ', paste0(group_names_full, collapse = ' & '), '\\\\'
  )

  stat_cols_clean <- stat_cols %>% gsub('\\.[a-z]$', '', .) %>%
    gsub('estimate', 'est', ., fixed = TRUE) %>%
    gsub('std.error', 'se', ., fixed = TRUE) %>%
    gsub('p.value', 'pval', ., fixed = TRUE) %>%
    map_chr(~ sprintf('\\multicolumn{1}{c}{%s}', .x))

  stat_cols_full <- paste0(
    ' & ', paste0(stat_cols_clean, collapse = ' & '), '\\\\'
  )

  table_pream <- c(
    sprintf(
      '\\begin{%s}{l%s}',
      env,
      paste0(rep('D{.}{.}{6}', total_cols - 1), collapse = '')
    ),
    if (!is.null(caption)) sprintf('\\caption{%s}\\\\', caption),
    if (!is.null(caption)) sprintf('\\label{%s}', caption),
    group_names_full,
    if (max(header_rle$lengths) > 1) stat_cols_full,
    '\\midrule',
    if (repeat_heading) {
      c(
        '\\endfirsthead',
        if (!is.null(caption)) {
          sprintf('\\caption*{%s (Continued)}\\\\', caption)
        },
        group_names_full,
        if (max(header_rle$lengths) > 1) stat_cols_full,
        '\\midrule',
        '\\endhead'
      )
    },
    if (bottom_rule) {
      c(
        '\\bottomrule',
        sprintf('\\multicolumn{%s}{r@{}}{continued \\ldots}\\\\', total_cols),
        '\\endfoot',
        '\\endlastfoot'
      )
    }
  )

  ltable <- vector('character', total_rows)
  labels <- names(reg_table)[grepl('label_', names(reg_table))]
  flevels <- names(reg_table)[grepl('flevels_', names(reg_table))]
  for (i in 1:total_rows) {
    label_current <- unlist(reg_table[i, labels], use.names = FALSE)
    label_current <- if (all(is.na(label_current))) {
      NA
    } else {
      discard(label_current, is.na) %>%
        paste(collapse = ' * ')
    }

    label_next <- unlist(reg_table[i + 1, labels], use.names = FALSE)
    label_next <- if (all(is.na(label_next))) {
      NA
    } else {
      discard(label_next, is.na) %>%
        paste(collapse = ' * ')
    }
    type <- reg_table$type[i]
    flevel <- unlist(reg_table[i, flevels], use.names = FALSE)
    flevel <- if (all(is.na(flevel))) {
      NA
    } else {
      discard(flevel, is.na) %>%
        paste(collapse = ' * ')
    }

    if (
      !is.na(label_next) &
      label_current != label_next &
      !(type %in% c('sumstat', 'sumstatN'))
    ) {
      end <- '\\\\[.5em]'
    } else {
      end <- '\\\\'
    }

    if (type == 'omitted') {
      ltable[i] <- sprintf(
        '\\multicolumn{%s}{l}{%s (%s)}\\\\',
        total_cols,
        sanitize(label_current),
        sanitize(flevel)
      )
    } else {
      row_label <- if (is.na(flevel)) label_current else flevel
      sanitized_row <- unlist(
        c(row_label, reg_table[i, non_title_cols]), use.names = FALSE
      ) %>%
        map_chr(~ ifelse(is.na(.x), '', .x))
      sanitized_row[2:length(sanitized_row)] <- gsub(
        '(\\*+)', '^{\\1}', sanitized_row[2:length(sanitized_row)]
      )
      ltable[i] <- paste0(paste0(sanitized_row, collapse = ' & '), end)
    }
  }

  c(
    table_pream,
    ltable,
    '\\bottomrule',
    sprintf('\\end{%s}', env)
  )
}
