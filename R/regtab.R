#' @importFrom tibble tibble rownames_to_column
#' @importFrom dplyr anti_join arrange bind_cols bind_rows case_when distinct
#'                   everything filter full_join funs group_by if_else left_join
#'                   mutate mutate_at matches mutate_if rename rename_at
#'                   rename_all rowwise select semi_join starts_with ungroup
#'                   vars group_indices slice right_join
#' @importFrom tibble rownames_to_column
#' @importFrom purrr keep map map_chr map2 map2_df reduce discard
#' @importFrom tidyr gather
#' @importFrom magrittr %>% extract
#' @importFrom rlang sym syms
#' @importFrom splitstackshape cSplit expandRows

get_core_levels <- function(xlevels) {
  if (length(xlevels) == 0) return(NULL)

  map2_df(
    names(xlevels), xlevels,
    ~ tibble(
        term = paste0(.x, .y),
        label_1 = .x,
        flevels_1 = .y
      ) %>%
        mutate_if(is.factor, as.character)
  ) %>%
    group_by(label_1) %>%
      mutate(type = ifelse(1:n() == 1, 'omitted', 'coef/se')) %>%
    ungroup()
}

get_interacted_levels <- function(term, xlevels) {
  split_interactions <- term[grep(':', term)]
  if (length(split_interactions) == 0) return(NULL)

  core_levels <- get_core_levels(xlevels) %>%
    select(-type)

  inter_table <- tibble(term = split_interactions) %>%
    cSplit('term', ':', type.convert = FALSE) %>%
    gather(interaction, term) %>%
    left_join(
      rename(core_levels, label = label_1, flevels = flevels_1),
      by = 'term'
    ) %>%
    mutate(label = ifelse(is.na(label), term, label)) %>%
    split(.$interaction) %>%
    map2(
      1:length(.),
      ~ .x %>%
          select(-interaction) %>%
          rename_all(paste0, '_', .y)
    ) %>%
    bind_cols() %>%
    select(starts_with('label'), starts_with('flevels')) %>%
    mutate(term = split_interactions) %>%
    select(term, everything())

  # Omitted interactions
  # Note: The split screws up the order, but I do not think it actually matters.
  omitted_core <- core_levels %>%
    group_by(label_1) %>%
    distinct(label_1, .keep_all = TRUE)
  labels <- names(inter_table)[grep('label_', names(inter_table))]
  seek_omit <- split(inter_table, group_indices(inter_table, !!!syms(labels)))

  inter_table <- map(
    seek_omit,
    ~ {
      o <- slice(.x, 1)
      ix <- 1
      for (l in labels) {
        match_omit <- match(o[, l], omitted_core$label_1)
        if (!is.na(match_omit)) {
          o[[paste0('flevels_', ix)]] <- omitted_core$flevels_1[[match_omit]]
        } else {
          o[[paste0('flevels_', ix)]] <- NA
        }
        ix <- ix + 1
      }
      term <- map_chr(
        1:length(labels),
        ~ {
          label <- paste0('label_', .x)
          level <- paste0('flevels_', .x)
          if (is.na(o[1, level])) {
            o[[label]]
          } else {
            paste0(o[[label]], o[[level]])
          }
        }
      ) %>%
        discard(~ is.na(.x)) %>%
        paste0(collapse = ':')
      o$term <- term
      o$type <- 'omitted'

      if (identical(o$term[1], .x$term[1])) .x else bind_rows(o, .x)
    }
  ) %>%
    bind_rows() %>%
    mutate(type = ifelse(is.na(type), 'coef/se', type))

  inter_table
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

  all_possible_vars <- if (nrow(all_possible_levels) == 0) {
    tidy_reg %>%
      mutate(
        label_1 = term,
        flevels_1 = NA_character_,
        type = 'coef/se'
      ) %>%
      select(term, label_1, flevels_1, type)
  } else {
    bind_rows(
      all_possible_levels,
      anti_join(tidy_reg, all_possible_levels, by = 'term') %>%
        select(term) %>%
        mutate(label_1 = term, type = 'coef/se')
    )
  }

  tidy_reg <- left_join(
    all_possible_vars,
    tidy_reg %>%
      mutate(order = 1:n()),
    by = 'term'
  ) %>%
    arrange(order)

  labels <- names(tidy_reg) %>% extract(grepl('label_', .))
  tidy_table <- filter(tidy_reg, !is.na(order)) %>%
    distinct(!!!syms(labels)) %>%
    mutate(unique_order = 1:n()) %>%
    right_join(tidy_reg, by = labels) %>%
    filter(!is.na(unique_order)) %>%
    mutate(order = ifelse(is.na(order), 0, order)) %>%
    arrange(unique_order, order) %>%
    select(-order, -unique_order)

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
    tidy_table <- tidy_table %>% mutate(est.sig = NA_character_)
  }

  # Properly reform sumstats.
  sumstats <- rownames_to_column(as.data.frame(sumstats)) %>%
    rename(label_1 = rowname, estimate = V1)

  # Return stacked tidy_table and sumstats. We will format the digits here
  # because this seems the most logical place to do it and only _then_ add the
  # sig stars (if desired), because we can't do it before we format.
  bind_rows(
    tidy_table,
    sumstats %>%
      mutate(type = ifelse(label_1 == 'N', 'sumstatN', 'sumstat'))
  ) %>%
    reg_format(digits = digits, exclude_n = exclude_n) %>%
    mutate(
      estimate = ifelse(is.na(est.sig), estimate, paste0(estimate, est.sig))
    ) %>%
    select(-matches('statistic|is_factor|est\\.sig')) %>%
    select(
      starts_with('label_'),
      starts_with('flevels_'),
      type,
      everything()
    )
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
      ~ {
        label_x <- names(.x)[grepl('label_', names(.x))]
        label_y <- names(.y)[grepl('label_', names(.y))]

        flevels_x <- names(.x)[grepl('flevels_', names(.x))]
        flevels_y <- names(.y)[grepl('flevels_', names(.y))]

        common_labels <- intersect(label_x, label_y)
        common_flevels <- intersect(flevels_x, flevels_y)

        full_join(.x, .y, by = c(common_labels, common_flevels, 'type'))
      }
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
  reg_table <- select(reg_table, -matches('p\\.value')) %>%
    mutate(count = ifelse(type == 'coef/se', 2, 1)) %>%
    expandRows('count')

  local({
    t <- vector(mode = 'character', nrow(reg_table))
    i <- 1
    while (TRUE) {
      if (reg_table$type[i] == 'coef/se') {
        t[i] <- 'coef'
        t[i + 1] <- 'se'
        i <- i + 2
      } else {
        t[i] <- reg_table$type[i]
        i <- i + 1
      }
      if (i > nrow(reg_table)) break
    }

    reg_table$type <<- t
  })

  reg_table %>%
    mutate(
      estimate = case_when(
        type == 'coef' ~ estimate, type == 'se' ~ std.error, TRUE ~ estimate
      )
    ) %>%
    select(-std.error)
}

#' @export
reg_remove_base <- function(reg_table, when = 'binary') {
  if (when == 'always') {
    filter(reg_table, type != 'omitted')
  } else if (when == 'binary') {
    labels <- names(reg_table)[grepl('label_', names(reg_table), fixed = TRUE)]
    reg_table %>%
      group_by(!!!syms(labels)) %>%
        mutate(
          bottom = any(type == 'se'),
          count = (n() - 1) / ifelse(bottom, 2, 1)
        ) %>%
        filter(!(type == 'omitted' & count == 1)) %>%
      ungroup() %>%
      select(-bottom, -count)
  }
}

