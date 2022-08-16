# R packages required for  the plot of number of hyperedges against clust. coef in project 1
packages = c('RcppArmadillo','Rcpp', 'geometry', 'mnormt', 'igraph', 'CVXR','cluster', 'ggplot2','ggpubr','factorextra','viridis','gridExtra','bezier')

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}



