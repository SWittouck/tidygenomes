filter_orthogroups <- function(tog, ...) {
  
  tog %>%
    modify_at("orthogroups", filter, ...)
  
}