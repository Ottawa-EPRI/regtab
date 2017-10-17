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
    # FIXME: The issue is HERE. It adds all these extra rows which may not
    # exist because of a full join. We need a join of ONLY what's actually
    # in the regression + any omitteds.
    terms_no_interactions <- bind_rows(possible_names) %>%
      full_join(
        terms_no_interactions %>%
          mutate(present = TRUE),
        by = 'term'
      )
  } else {
    terms_no_interactions %<>%
      mutate(
        label = NA_character_, levels = NA_character_, level_order = NA_integer_
      )
  }

  terms_no_interactions2 <- terms_no_interactions %>%
    mutate(label = ifelse(is.na(label), term, label)) %>%
    group_by(label) %>%
      filter(any(present)) %>%
    ungroup()

  if (!any(grepl(':', tidy_reg$term))) {
    table <- terms_no_interactions
  } else {
    interactors <- tidy_reg$term[grep(':', tidy_reg$term)]
    # FIXME: omitted interactors, level order!
    table_interactors <- Map(
      function(x) {
        int_terms <- Map(
          function(y) {
            left_join(
              data.frame(term = as.character(y)),
              terms_no_interactions,
              by = 'term'
            )
          },
          x
        )
        tibble(
          term = Reduce(function(x, y) paste0(x$term, ':', y$term),
                        int_terms),
          label = Reduce(function(x, y) paste0(x$label, ' * ', y$label),
                         int_terms),
          levels = Reduce(
            function(x, y) {
              if (is.na(x$levels) & is.na(y$levels)) NA_integer_
              else if (!is.na(x$levels) & is.na(y$levels)) paste0(x$levels, ' * ', y$label)
              else if (is.na(x$levels) & !is.na(y$levels)) paste0(x$label, ' * ', y$levels)
              else paste0(x$levels, ' * ', y$levels)
            },
            int_terms
          ),
          omitted = Reduce(
            function(x, y) {
              x_omit <- terms_no_interactions %>% filter(label == x$label, level_order == 1)
              y_omit <- terms_no_interactions %>% filter(label == y$label, level_order == 1)
              paste0(x_omit$levels, ' * ', y_omit$levels)
            },
            int_terms
          )
        )
      },
      strsplit(interactors, ':')
    ) %>%
      bind_rows()

    omitted_table_interactors <- table_interactors %>%
      distinct(omitted, .keep_all = TRUE) %>%
      select(label, levels = omitted) %>%
      mutate(level_order = 1)
    table_interactors %<>%
      select(-omitted) %>%
      full_join(omitted_table_interactors)

    table <- full_join(terms_no_interactions2, table_interactors)
  }

  in_model <- tidy_reg %>%
    select(term) %>%
    left_join(table, by = 'term')

  omitted_levels <- table %>%
    filter(level_order == 1)

  if (nrow(omitted_levels) > 0) {
    in_model %<>% bind_rows(omitted_levels)
    unique_label_order <- in_model %>%
      distinct(label) %>%
      mutate(label_ix = 1:n())
    in_model %<>%
      left_join(unique_label_order, by = 'label') %>%
      arrange(label_ix, level_order) %>%
      select(-label_ix, -present)
  }
  in_model
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
    rename(label = rowname, estimate = V1)
  table <- level_df %>%
    left_join(tidy_reg, by = 'term') %>%
    bind_rows(sumstats) %>%
    select(-matches('term|statistic|p.value|label_ix|est.sig')) %>%
    select(label, levels, level_order, everything())

  table
}

z <- lm(Sepal.Length ~ factor(Sepal.Width), data = iris)
z <- lm(Sepal.Length ~ Sepal.Width, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width) * Species, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width):Species, data = iris)
