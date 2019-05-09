process_new_orthogroup_names <- function(tog) {

  orthogroups_names <- tog$orthogroups %>%
    select(orthogroup, orthogroup_new)

  tog$genes <- tog$genes %>%
    left_join(orthogroups_names) %>%
    mutate(orthogroup = orthogroup_new) %>%
    select(- orthogroup_new)

  tog$orthogroups <- tog$orthogroups %>%
    mutate(orthogroup = orthogroup_new) %>%
    select(- orthogroup_new)

  tog

}