#' @keywords internal
meta_names <- function(component) {
  names_list <- list(
    data = c("yi", "yif", "vi", "affirm", "cluster"),
    values = c("selection_ratio", "selection_tails", "model_type",
               "favor_positive", "alpha_select", "ci_level", "small", "k",
               "k_affirmative", "k_nonaffirmative"),
    stats = c("estimate", "se", "ci_lower", "ci_upper", "p_value"))
  names_list[[component]]
}

#' @keywords internal
meta_names_str <- function(component) {
  cnames <- meta_names(component)
  paste(paste0("`", cnames, "`"), collapse = ", ")
}

#' @keywords internal
svalue_names <- function(component) {
  names_list <- list(
    data = c("yi", "vi", "affirm", "cluster"),
    values = c("q", "model_type", "favor_positive", "alpha_select",
               "ci_level", "small", "selection_ratio_max", "k",
               "k_affirmative", "k_nonaffirmative"),
    stats = c("sval_est", "sval_ci"))
  names_list[[component]]
}

#' @keywords internal
svalue_names_str <- function(component) {
  cnames <- svalue_names(component)
  paste(paste0("`", cnames, "`"), collapse = ", ")
}
