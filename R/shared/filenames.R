# Returns the filename to use for a figure and given project and plot name.
figure_filename <- function(project, topic, name) {
  filename = paste0('output/figures/', project, '/', topic, '/', name, '.png')
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  return(filename)
}

# Returns the filename to use for a table (CSV) and given project, topic and name.
table_filename <- function(project, topic, name) {
  filename = paste0('output/tables/', project, '/', topic, '/', name, '.csv')
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  return(filename)
}
