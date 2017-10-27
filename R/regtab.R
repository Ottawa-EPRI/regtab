get_core_levels <- function(xlevels) {
   bind_rows(
     tibble(
       term = character(),
       label = character(),
       flevels = character(),
       level_order = integer()
     ),
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
    select(-label_order, -term)

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
  } else {
    tidy_table %>% mutate(est.sig = NA_character_)
  }

  # Properly reform sumstats.
  sumstats <- rownames_to_column(as.data.frame(sumstats)) %>%
    rename(label = rowname, estimate = V1)

  # Return stacked tidy_table and sumstats. We will format the digits here
  # because this seems the most logical place to do it and only _then_ add the
  # sig stars (if desired), because we can't do it before we format.
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
    reg_format(digits = digits, exclude_n = exclude_n) %>%
    mutate(
      estimate = ifelse(is.na(est.sig), estimate, paste0(estimate, est.sig))
    ) %>%
    select(-matches('statistic|is_factor|est\\.sig')) %>%
    select(label, flevels, level_order, type, everything())
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

reg_combine <- function(reg_list) {
  regs <- map2(
    reg_list, letters[1:length(reg_list)],
    ~ rename_at(
        .x,
        vars(matches('estimate|std\\.error|p\\.value|est\\.sig')),
        function(r) paste0(r, '.', .y)
      )
  )
  if (length(regs) > 1) {
    combined_regs <- reduce(
      regs,
      ~ full_join(.x, .y, by = c('label', 'flevels', 'level_order', 'type'))
    )
  } else {
    combined_regs <- regs
  }
  sumstats <- combined_regs %>%
    filter(type %in% c('sumstat', 'sumstatN'))
  combined_regs %>%
    filter(!type %in% c('sumstat', 'sumstatN')) %>%
    bind_rows(sumstats)
}

reg_bottom_se <- function(reg_table, p.value = FALSE) {
  label_tibble <- tibble(
    label = unique(reg_table$label), order = seq_along(label)
  )

  if (p.value) {
    reg_table <- reg_table %>%
      mutate(
        std.error = ifelse(
          !is.na(std.error) & !is.na(p.value),
          paste0('(', std.error, ', ', p.value, ')'),
          ifelse(!is.na(std.error), paste0('(', std.error, ')'), NA)
        )
      )
  } else {
    reg_table <- reg_table %>%
      mutate(std.error = ifelse(
        !is.na(std.error), paste0('(', std.error, ')'), NA)
      )
  }
  reg_table <- select(reg_table, -matches('p\\.value'))

  gather(
    reg_table, type2, estimate, -label, -flevels, -level_order, -type
  ) %>%
    left_join(label_tibble, by = 'label') %>%
    arrange(order, level_order, type2) %>%
    filter(
      !(type %in% c('omitted', 'sumstat', 'sumstatN') & type2 == 'std.error')
    ) %>%
    mutate(type = case_when(
      type == 'coef/se' & type2 == 'estimate' ~ 'coef',
      type == 'coef/se' & type2 == 'std.error' ~ 'se',
      TRUE ~ type
    )) %>%
    select(-order, -type2)
}

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

  stat_cols <- names(reg_table) %>% `[`(grepl('\\.[a-z]$', .))
  model_groupings <- map_chr(stat_cols, ~ substr(.x, nchar(.x), nchar(.x)))
  header_rle <- rle(model_groupings)
  # We need to fix this to take user input if desired.
  if (is.null(group_names)) {
    group_names <- map_chr(seq_along(header_rle$lengths), ~ sprintf('(%s)', .x))
  } else {
    group_names <- group_names %>% `[`(1:length(header_rle$lengths))
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
