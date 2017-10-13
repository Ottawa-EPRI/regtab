retrieve_factor_labels <- function(reg) {
  # FIXME: What happens when there are no levels
  #reg_labels <- attributes(reg$terms)$term.labels
  reg_factor_levels <- reg$xlevels

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
  bind_rows(possible_names)
}

regtab <- function(
  reg,
  digits = 3,
  pvals = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  sumstat_vars = c('adj.r.squared', 'N')
) {
  if (length(reg$xlevels) != 0) {
    levels <- retrieve_factor_labels(reg)
  }

  sumstats <- t(cbind(broom::glance(reg), N = nobs(reg)))
  tidy_reg <- broom::tidy(reg)
  sorted_p <- sort(pvals, decreasing = TRUE)
  for (i in seq_along(sorted_p)) {
    tidy_reg %<>% mutate(
      est.sig = if_else(p.value <= sorted_p[i], names(sorted_p)[i], '')
    )
  }

  tidy_reg %<>%
    mutate(
      estimate = paste0(formatC(estimate, digits = digits, format = 'f'),
                        est.sig),
      std.error = formatC(std.error, digits = digits, format = 'f'),
      type = 'coef/se'
    ) %>%
    select('term', 'estimate', 'std.error', 'type')

  sumstats <- rownames_to_column(as.data.frame(sumstats)) %>%
    mutate(
      V1 = paste0(formatC(V1, digits = digits, format = 'f')),
      type = 'sumstat'
    ) %>%
    mutate_all(as.character) %>%
    filter(rowname %in% sumstat_vars) %>%
    rename(term = rowname, estimate = V1)
  table <- bind_rows(tidy_reg, sumstats)

  if (length(reg$xlevels) != 0) {
    table %<>%
      full_join(levels, by = 'term')
  } else {
    table %<>%
      mutate(
        label = NA_character_,
        levels = NA_character_,
        level_order = NA_integer_
      )
  }

  table %<>%
    mutate(
      label = ifelse(is.na(label), term, label),
      type = ifelse(is.na(type) & level_order == 1, 'omitted', type)
    )

  table %>%
    select(-term)
}


