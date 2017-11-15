#' @importFrom tibble tibble rownames_to_column
#' @importFrom dplyr anti_join arrange bind_cols bind_rows case_when distinct
#'                   everything filter full_join funs group_by if_else left_join
#'                   mutate mutate_at matches mutate_if rename rename_at
#'                   rename_all rowwise select semi_join ungroup vars
#' @importFrom tibble rownames_to_column
#' @importFrom purrr keep map map_chr map2 map2_df reduce
#' @importFrom tidyr gather unite
#' @importFrom magrittr %>% extract
#' @importFrom rlang syms

paste_0 <- function(..., collapse = ' * ') {
  dots <- list(...)
  dots <- dots %>% extract(!is.na(.))
  paste0(dots, collapse = collapse)
}

any_NA <- function(...) {
  dots <- list(...)
  if (all(is.na(dots))) return(NA)

  dots <- dots %>% extract(!is.na(.))
  any(as.logical(dots))
}

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
  core_levels <- get_core_levels(xlevels) %>%
    select(-level_order)
  split_interactions <- term[grep(':', term)]
  if (length(split_interactions) == 0) return(NULL)

  inter_table <- tibble(term = split_interactions) %>%
    splitstackshape::cSplit('term', ':', type.convert = FALSE) %>%
    gather(interaction, term) %>%
    left_join(core_levels, by = 'term') %>%
    mutate(
      is_factor = ifelse(!is.na(term), !is.na(flevels), NA),
      label = ifelse(is.na(label), term, label),
      flevels = ifelse(is.na(flevels), label, flevels)
    ) %>%
    split(.$interaction) %>%
    map2(
      1:length(.),
      ~ .x %>%
          select(-interaction) %>%
          rename_all(paste0, '_', .y)
    ) %>%
    bind_cols()

  termnames <- names(inter_table) %>% extract(grepl('term', .))
  labnames <- names(inter_table) %>% extract(grepl('label', .))
  flevelnames <- names(inter_table) %>% extract(grepl('flevels', .))
  is_factornames <- names(inter_table) %>% extract(grepl('is_factor', .))

  inter_table <- inter_table %>%
    rowwise() %>%
      mutate(
        term = paste_0(!!!syms(termnames), collapse = ':'),
        label = paste_0(!!!syms(labnames)),
        flevels = paste_0(!!!syms(flevelnames)),
        is_factor = any_NA(!!!syms(is_factornames))
      ) %>%
    ungroup() %>%
    unite(flevel_factors, !!!syms(is_factornames)) %>%
    select(-matches('_[0-9]')) %>%
    filter(is_factor) %>%
    select(-is_factor)

  if (nrow(inter_table) > 0) {
    inter_table <- inter_table %>%
      group_by(label) %>%
        mutate(level_order = (1L:n() + 1L)) %>%
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

#' @export
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

#' @export
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

#' @export
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
    reg_table, type2, estimate, -label, -flevels, -level_order, -type,
    -matches('flevel_factors')
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

#' @export
reg_remove_base <- function(reg_table, when = 'binary') {
  if (when == 'always') {
    filter(reg_table, type != 'omitted')
  } else if (when == 'binary') {
    reg_table %>%
      group_by(label) %>%
        filter(
          !(all(!is.na(flevels)) & max(level_order) == 2 & type == 'omitted')
        ) %>%
      ungroup()
  }
}

#' @export
reg_add_labels <- function(reg_table, label_list) {
  labels <- map(strsplit(reg_table$label, '*', fixed = TRUE), trimws)
  label_matches <- map(labels, match, names(label_list))
  label_list_flat <- set_names(unlist(label_list), NULL)
  labels_changed <- map(label_matches, ~ label_list_flat[.x])
  labels_new <- map2(labels_changed, labels, ~ ifelse(is.na(.x), .y, .x))
  labels_rebuild <- map_chr(labels_new, paste0, collapse = ' * ')
  reg_table$label <- labels_rebuild

  if (is.null(reg_table$flevel_factors)) return(reg_table)

  levels <- map(strsplit(reg_table$flevels, '*', fixed = TRUE), trimws)
  level_factors <-  map(
    strsplit(reg_table$flevel_factors, "_", fixed = TRUE), as.logical
  )
  levels_new <- pmap(
    list(level_factors, levels, labels_new), ~ ifelse(..1, ..2, ..3)
  )

  levels_rebuild <- map_chr(
    levels_new,
    ~ {
      modify_vector <- keep(.x, negate(is.na))
      ifelse(
        length(modify_vector) != 0, paste0(modify_vector, collapse = ' * '), NA
      )
    }
  )
  mutate(reg_table,
         flevels = ifelse(type != 'omitted', levels_rebuild, flevels))
}
