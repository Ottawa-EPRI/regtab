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
  split_interactions <- term[grep(':', term)] %>%
    strsplit(':', fixed = TRUE)
  core_levels <- get_core_levels(xlevels)
  interact_tibble <- map(
    split_interactions,
    ~ map(
      .x,
      ~ tibble(term = .x) %>%
          left_join(core_levels, by = 'term') %>%
          mutate(
            is_factor = !is.na(flevels),
            label = ifelse(is.na(label), term, label),
            flevels = ifelse(is.na(flevels), label, flevels)
          )
    )
  )

  interact_reduce <- map(
    interact_tibble,
    ~ reduce(
        .x,
        ~ tibble(
            term = sprintf('%s:%s', .x$term, .y$term),
            label = sprintf('%s * %s', .x$label, .y$label),
            flevels = sprintf('%s * %s', .x$flevels, .y$flevels),
            level_order = NA_integer_,
            is_factor = any(.x$is_factor, .y$is_factor)
          )
    )
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

  lv <- map_chr(
    strsplit(inter_table$label, ' * ', fixed = TRUE),
    ~ map(.x, ~ xlevels[[.x]][1]) %>%
        keep(~ !is.null(.x)) %>%
        paste0(collapse = ' * ')
  )

  tibble(label = inter_table$label, flevels = lv, level_order = 1)
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
