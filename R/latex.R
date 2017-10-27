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

  # We need to 'simplify' the labels appropriately.
  reg_table_min <- mutate(
    reg_table,
    label = case_when(
      !is.na(label) & !is.na(flevels) & type %in% c('coef', 'coef/se') ~
        flevels,
      !is.na(label) & !is.na(flevels) & type == 'omitted' ~
        sprintf('%s (%s)', label, flevels),
      !is.na(label) & type == 'se' ~ '',
      TRUE ~ label
    )
  ) %>%
    select(-matches('flevels|level_order|type')) %>%
    mutate_all(~ ifelse(is.na(.x), '', .x))

  total_rows <- nrow(reg_table_min)
  total_cols <- ncol(reg_table_min)

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
  for (i in 1:total_rows) {
    label_current <- reg_table$label[i]
    label_next <- reg_table$label[i + 1]
    type <- reg_table$type[i]

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
        '\\multicolumn{%s}{l}{%s}\\\\', total_cols, label_current
      )
    } else {
      sanitized_row <- sanitize(reg_table_min[i, ])
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
