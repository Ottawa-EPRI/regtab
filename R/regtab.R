library(purrr)
library(dplyr)
library(magrittr)
library(tibble)

get_core_levels <- function(xlevels) {
   c(chr, int) %<-% list(character(), integer())
   bind_rows(
     tibble(term = chr, label = chr, flevels = chr, level_order = int),
     map2_df(
       names(xlevels), xlevels,
       ~ tibble(term = paste0(.x, .y), label = .x, flevels = .y) %>%
           mutate(level_order = 1:n()) %>%
           mutate_if(is.factor, as.character)
     )
   )
}

get_interacted_levels <- function(term, xlevels) {
  split_interactions <- term[grep(':', term)] %>%
    strsplit(':', fixed = TRUE)
  if (length(split_interactions) == 0) {
    return(NULL)
  }

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
  pvals = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  digits = 3,
  exclude_n = TRUE
) {

  # Get the initial model and summary statistics.
  sumstats <- t(cbind(broom::glance(reg), N = nobs(reg)))
  tidy_reg <- broom::tidy(reg)

  # Get all possible levels and interacted levels. Get the list of continous
  # variables as well, and add label to equal to term. This way we have a list
  # of all possible variables. We'll also get a separate df of just the omitted
  # variables.
  core_levels <- get_core_levels(reg$xlevels)
  interacted_levels <- get_interacted_levels(tidy_reg$term, reg$xlevels)
  all_possible_levels <- bind_rows(core_levels, interacted_levels)
  continous_vars <- anti_join(tidy_reg, all_possible_levels, by = 'term') %>%
    select(term) %>%
    mutate(label = term)
  all_possible_vars <- bind_rows(all_possible_levels, continous_vars)
  all_possible_omitted_vars <- all_possible_vars %>%
    filter(level_order == 1)

  # Join the tidy_reg to all possible vars, thereby getting proper labels and
  # level orders. We need to separately look for which omitted terms need to
  # come along, and then bind those together.
  tidy_reg <- tidy_reg %>%
    left_join(all_possible_vars, by = 'term')
  omitted_vars <- semi_join(all_possible_omitted_vars, tidy_reg, by = 'label')
  tidy_table <- bind_rows(tidy_reg, omitted_vars)

  # Since any of the the omitted labels will come at the bottom, we need to
  # have a separate table with the proper label order (i.e., the original
  # tidy_reg label order), and then use the unique order to order the labels,
  # using level order as a second categorization, thus properly ordering the
  # omitted vars (if any). Then we can clean up the column order.
  unique_labels <- unique(tidy_table$label)
  label_order <- tibble(
    label = unique_labels, label_order = seq_along(unique_labels)
  )
  tidy_table <- left_join(tidy_table, label_order, by = 'label') %>%
    arrange(label_order, level_order) %>%
    select(-label_order, -term) %>%
    select(label, flevels, level_order, everything())

  # Add est.sig stars if desired.
  if (!is.null(pvals)) {
    sorted_p <- sort(pvals, decreasing = TRUE)
    for (i in seq_along(sorted_p)) {
      tidy_table <- tidy_table %>% mutate(
        est.sig = if_else(
          p.value <= sorted_p[i], names(sorted_p)[i], NA_character_
        )
      )
    }
  }

  # Properly reform sumstats.
  sumstats <- rownames_to_column(as.data.frame(sumstats)) %>%
    rename(label = rowname, estimate = V1)

  # Return stacked tidy_table and sumstats.
  bind_rows(
    tidy_table %>%
      mutate(
        type = ifelse(
          level_order == 1 & !is.na(level_order), 'omitted', 'coef/se'
        )
      ),
    sumstats %>%
      mutate(type = ifelse(label == 'N', 'sumstatN', 'sumstat'))
  ) %>%
    select(-matches('statistic|is_factor')) %>%
    reg_format(digits = digits, exclude_n = exclude_n)
}

reg_format <- function(
  reg_table,
  digits = 3,
  exclude_n = TRUE
) {
  match_vars <- vars(matches('estimate|std\\.error|p\\.value'))

  na_formatC <- function(x) {
    ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = 'f'))
  }

  if (exclude_n) {
    mutate_at(reg_table, match_vars,
              funs(ifelse(type != 'sumstatN', na_formatC(.), .)))
  } else {
    mutate_at(reg_table, match_vars, na_formatC)
  }
}


z <- lm(Sepal.Length ~ factor(Sepal.Width), data = iris)
z <- lm(Sepal.Length ~ Sepal.Width, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width) * Species, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width):Species, data = iris)
z <- lm(Sepal.Length ~ Sepal.Width:Petal.Width, data = iris)
z <- lm(Sepal.Length ~ factor(Sepal.Width):Species + Species:Petal.Length + factor(Sepal.Width):Species:Petal.Length, data = iris)
