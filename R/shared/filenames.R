# Returns the filename to use for a figure and given project and plot name.
figure_filename <- function(project, topic, name) {
  return(paste0('output/figures/', project, '/', topic, '/', name, '.png'))
}

# Returns the filename to use for a table (CSV) and given project, topic and name.
table_filename <- function(project, topic, name) {
  return(paste0('output/tables/', project, '/', topic, '/', name, '.csv'))
}
