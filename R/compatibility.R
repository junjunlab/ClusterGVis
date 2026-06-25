pop_legacy_arg <- function(value, dots, legacy, current_missing) {
  if (legacy %in% names(dots)) {
    if (current_missing) {
      value <- dots[[legacy]]
    }

    dots[[legacy]] <- NULL
  }

  list(value = value, dots = dots)
}
