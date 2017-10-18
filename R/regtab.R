library(purrr)
library(dplyr)
library(magrittr)
library(tibble)

get_core_levels <- function(xlevels) {
   xlevels <- bind_rows(map2_df(
    names(xlevels), xlevels,
    ~ tibble(term = paste0(.x, .y), label = .x, flevels = .y) %>%
        mutate(level_order = 1:n()) %>%
        mutate_if(is.factor, as.character)
  ))
  if (length(xlevels) == 0) {
    c(chr, int) %<-% list(character(), integer())
    xlevels <- tibble(term = chr, label = chr, flevels = chr, level_order = int)
  }
  xlevels
}

get_interacted_levels <- function(term, xlevels) {
  interactions <- term[grep(':', term)]
  split_interactions <- strsplit(interactions, ':')
  core_levels <- bind_rows(get_core_levels(xlevels))
  interact_tibble <- Map(
    function(x)
    Map(
      function(y) {
        tibble(term = y) %>%
          left_join(core_levels, by = 'term') %>%
          mutate(
            is_factor = !is.na(flevels),
            label = ifelse(is.na(label), term, label),
            flevels = ifelse(is.na(flevels), label, flevels)
          )
      },
      x,
      USE.NAMES = FALSE
    ),
    split_interactions
  )

  interact_reduce <- Map(
    function(x)
      Reduce(
        function(a, b) {
          tibble(
            term = sprintf('%s:%s', a$term, b$term),
            label = sprintf('%s * %s', a$label, b$label),
            flevels = sprintf('%s * %s', a$flevels, b$flevels),
            level_order = NA_integer_,
            is_factor = any(a$is_factor, b$is_factor)
          )
        },
        x
      ),
    interact_tibble
  )
  inter_table <- bind_rows(interact_reduce) %>%
    filter(is_factor != FALSE) %>%
    select(-is_factor)

  if (nrow(inter_table) > 0) {
    inter_table <- inter_table %>%
      group_by(label) %>%
        mutate(level_order = (1:n() + 1)) %>%
      ungroup()
    omitted <- get_interacted_omitted(inter_table, xlevels)
    inter_table <- bind_rows(inter_table, omitted)
  }
  inter_table
}

get_interacted_omitted <- function(inter_table, xlevels) {
  inter_table <- inter_table %>%
    distinct(label)

  table_labels <- strsplit(inter_table$label, ' * ', fixed = TRUE)
  lv <- Map(
    function(x) {
      omitteds <- Map(
        function(y) {
          match_ix <- match(y, names(xlevels))
          if (!is.na(match_ix)) xlevels[[match_ix]][1] else NULL
        },
        x
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
      omitteds <- Filter(function(x) !is.null(x), omitteds)
      paste0(omitteds, collapse = ' * ')
    }, ss
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

  tibble(label = inter_table$label, flevels = unlist(lv), level_order = 1)
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
      select(-label_ix, -present) %>%
      group_by(label) %>%
        mutate(level_order = add_interact_level_order(level_order)) %>%
      ungroup()
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
z <- lm(Sepal.Length ~ Sepal.Width:Petal.Width, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width):Species + Species:Petal.Length + factor(Sepal.Width):Species:Petal.Length, data = iris)
