# Returns the filename to use for a figure and given project and plot name.
figure_filename <- function(project, name) {
  return(paste0('output/figures/', project, '/', name, '.png'))
}
