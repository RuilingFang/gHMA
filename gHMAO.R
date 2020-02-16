#############################################################################################################################
################     An omnibus test for gene-based high-dimensional mediation analysis (gHMA-O)
#############################################################################################################################
#' An omnibus test for gene-based high-dimensional mediation analysis (gHMA-O)
#'
#' @param p1 an output object of the gHMAL function.
#' @param p2 an output object of the gHMANL function.
#'
#' @return A combined p-value. 
#' 

gHMAO <- function(p1, p2)
{
  pval <- t(c(p1$pcom, p2$pcom))
  
  pcom <- 1-pcauchy(sum(tan(pi*(0.5-pval))/ncol(pval)))
  
  return(pcom)
}