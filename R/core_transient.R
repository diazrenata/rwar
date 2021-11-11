#' Remove transient and/or non-core species
#'
#' @param dat matss-style
#' @param core_only if T, returns only core (present in > 2/3 of samples). if F, removes transients (present in < 1/3 of samples).
#'
#' @return matss-style but filtered
#' @export
#'
#' @importFrom dplyr mutate filter
core_transient <- function(dat, core_only = F){


  sp_occurrences <- dat$abundance > 0

  sp_presences <- data.frame(
    sp = colnames(sp_occurrences),
    noccurrences = colSums(sp_occurrences),
    nyears = nrow(sp_occurrences)
  ) %>%
    dplyr::mutate(prop_occurrences = noccurrences / nyears) %>%
    dplyr::mutate(is_transient = prop_occurrences < 1/3,
                  is_core = prop_occurrences > 2/3) # following Coyle et al.

  not_transient <- sp_presences %>%
    dplyr::filter(!is_transient)

  core <- sp_presences %>%
    dplyr::filter(is_core)


  if(core_only) {

    sp_to_keep <- core$sp

  } else {

    sp_to_keep <- not_transient$sp
  }

  dat_out <- dat

  dat_out$abundance <- dat_out$abundance[ , sp_to_keep]
  dat_out$metadata$species_table <- dat_out$metadata$species_table %>%
    dplyr::filter(id %in% sp_to_keep)

  dat_out
}

#' Remove transient and/or non-core species
#'
#' @param dat matss-style
#' @param core_only if T, returns only core (present in > 2/3 of samples). if F, removes transients (present in < 1/3 of samples).
#' @param null_seed 1989
#'
#' @return matss-style but filtered
#' @export
#'
#' @importFrom dplyr mutate filter
core_transient_null <- function(dat, core_only = F, null_seed = 1989){


  sp_occurrences <- dat$abundance > 0

  sp_presences <- data.frame(
    sp = colnames(sp_occurrences),
    noccurrences = colSums(sp_occurrences),
    nyears = nrow(sp_occurrences)
  ) %>%
    dplyr::mutate(prop_occurrences = noccurrences / nyears) %>%
    dplyr::mutate(is_transient = prop_occurrences < 1/3,
                  is_core = prop_occurrences > 2/3) # following Coyle et al.

  # now randomly reassign species labels

  set.seed(null_seed)
  sp_presences$sp <- sample(sp_presences$sp, size = length(sp_presences$sp), replace = F)
  set.seed(NULL)

  not_transient <- sp_presences %>%
    dplyr::filter(!is_transient)

  core <- sp_presences %>%
    dplyr::filter(is_core)


  if(core_only) {

    sp_to_keep <- core$sp

  } else {

    sp_to_keep <- not_transient$sp
  }

  dat_out <- dat

  dat_out$abundance <- dat_out$abundance[ , sp_to_keep]
  dat_out$metadata$species_table <- dat_out$metadata$species_table %>%
    dplyr::filter(id %in% sp_to_keep)

  dat_out
}


#' Get only transient or non-core speices
#'
#' @param dat matss-style
#' @param transient_only if T, returns only transient (present in <1/3 of samples). if F, removes core (present in > 2/3 of samples).
#'
#' @return matss-style but filtered
#' @export
#'
#' @importFrom dplyr mutate filter
just_transient <- function(dat, transient_only = T){


  sp_occurrences <- dat$abundance > 0

  sp_presences <- data.frame(
    sp = colnames(sp_occurrences),
    noccurrences = colSums(sp_occurrences),
    nyears = nrow(sp_occurrences)
  ) %>%
    dplyr::mutate(prop_occurrences = noccurrences / nyears) %>%
    dplyr::mutate(is_transient = prop_occurrences < 1/3,
                  is_core = prop_occurrences > 2/3) # following Coyle et al.

  not_core <- sp_presences %>%
    dplyr::filter(!is_core)

  transient <- sp_presences %>%
    dplyr::filter(is_transient)


  if(transient_only) {

    sp_to_keep <- transient$sp

  } else {

    sp_to_keep <- not_core$sp
  }

  dat_out <- dat

  dat_out$abundance <- dat_out$abundance[ , sp_to_keep]
  dat_out$metadata$species_table <- dat_out$metadata$species_table %>%
    dplyr::filter(id %in% sp_to_keep)

  dat_out
}

