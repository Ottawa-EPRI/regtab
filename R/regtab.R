retrieve_labels <- function(reg, tidy_reg) {
  reg_factor_levels <- reg$xlevels

  # This has clear problems outside base models, like mlogit that uses ':' for
  # different outcomes. We'll need to think about how to deal with that. Maybe
  # with just a class check, I doubt we want S3 methods all over the place, but
  # maybe.
  terms_no_interactions <- tidy_reg %>%
    filter(!grepl(':', term)) %>%
    select(term)

  if (length(reg_factor_levels) != 0) {
    possible_names <- Map(
      function(x, y) {
        data.frame(
          term = paste0(x, y),
          label = x,
          levels = y
        ) %>%
          mutate(level_order = 1:n()) %>%
          mutate_if(is.factor, as.character)
      },
      names(reg_factor_levels),
      reg_factor_levels
    )
    terms_no_interactions <- bind_rows(possible_names) %>%
      full_join(terms_no_interactions, by = 'term')
  } else {
    terms_no_interactions %<>%
      mutate(
        label = NA_character_, levels = NA_character_, level_order = NA_integer_
      )
  }

  terms_no_interactions %>%
      mutate(label = ifelse(is.na(label), term, label))
}

regtab <- function(
  reg,
  digits = 3,
  pvals = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  sumstat_vars = c('adj.r.squared', 'N')
) {

  sumstats <- t(cbind(broom::glance(reg), N = nobs(reg)))
  tidy_reg <- broom::tidy(reg)

  level_df <- retrieve_labels(reg, tidy_reg)
  if (!is.null(pvals)) {
    sorted_p <- sort(pvals, decreasing = TRUE)
    for (i in seq_along(sorted_p)) {
      tidy_reg %<>% mutate(
        est.sig = if_else(p.value <= sorted_p[i], names(sorted_p)[i], '')
      )
    }
  }

  tidy_reg %<>%
    mutate(
      estimate = formatC(estimate, digits = digits, format = 'f'),
      std.error = formatC(std.error, digits = digits, format = 'f'),
      type = 'coef/se'
    )
  if ('est.sig' %in% names(tidy_reg)) {
    tidy_reg %<>%
      mutate(
        estimate = ifelse(
          is.na(est.sig) | est.sig == '', estimate, paste0(estimate, est.sig)
        )
      )
  }

  sumstats <- rownames_to_column(as.data.frame(sumstats)) %>%
    mutate(
      V1 = paste0(formatC(V1, digits = digits, format = 'f')),
      type = 'sumstat'
    ) %>%
    mutate_all(as.character) %>%
    filter(rowname %in% sumstat_vars) %>%
    rename(term = rowname, estimate = V1)
  table <- bind_rows(tidy_reg, sumstats) %>%
    full_join(level_df, by = 'term')

  # FIXME: What if there are NO LEVELS
  # FIXME X 2
  #if (exists("level_df")) {
  #  table %<>%
  #    full_join(level_df, by = 'term')
  #} else {
  #  table %<>%
  #}

  table %<>%
    mutate(
      label = ifelse(is.na(label), term, label),
      type = ifelse(is.na(type) & level_order == 1, 'omitted', type)
    )

  labels <- unique(table$label)
  table %<>%
    left_join(tibble(label = labels, label_ix = seq_along(labels)),
              by = 'label') %>%
    arrange(label_ix, level_order)

  table %<>%
    select(-matches('term|statistic|p.value|label_ix|est.sig')) %>%
    select(label, levels, level_order, everything())

  table
}

z <- lm(Sepal.Length ~ factor(Sepal.Width), data = iris)
z <- lm(Sepal.Length ~ Sepal.Width, data = iris)
