\name{Continuous}
\alias{Continuous}
\title{Run BayesTraits Continuous}
\usage{
Continuous(tree, data, mode = "ML", regression = FALSE, directional = FALSE, lambda = 1, kappa = 1, delta = 1, 
ou = 0, force.correlation.to.0 = FALSE, mlt = 10, it = 100000, bi = 5000, sa = 100, silent=TRUE, rm=TRUE)
}
\description{
	Runs BayesTraits Continuous models.
}
\arguments{
  \item{tree}{an object of class phylo or multiPhylo}
  \item{data}{an object of class data.frame. The first column contains species names and subsequent columns contain continuous trait values. Species names must match those in the phylogeny.}
  \item{mode}{takes one of two values, "ML" to run models in maximum likelihood mode or "Bayesian" to run models in Bayesian MCMC mode. Defaults to "ML".}
  \item{regression}{logical indicating whether to fit a regression model. Only valid with multiple traits. Assumes first trait in the dataframe is the response variable. If this is TRUE, the directional parameter is overridden. Defaults to FALSE.}
  \item{directional}{logical indicating whether to fit directional evolution parameter. Defaults to FALSE.}
  \item{lambda}{value from 0 to 1 for the lambda parameter, or "ML" to estimate the maximum likelihood value lambda. Defaults to 1.}
  \item{kappa}{value from 0 to 3 for the kappa parameter, or "ML" to estimate the maximum likelihood value kappa. Defaults to 1.}
  \item{delta}{value from 0 to 3 for the delta parameter, or "ML" to estimate the maximum likelihood value delta. Defaults to 1.}
  \item{ou}{value of the Ornstein-Uhlenbeck "rubber band" parameter, reflecting the strength of stabilizing selection, or "ML" to estimate the maximum likelihood value of the parameter. Defaults to 0.}
  \item{force.correlation.to.0}{logical indicated whether to force the correlation between traits to be 0 for the purpose of model testing. Ignored if the number of traits is less than or greater than 2. Defaults to FALSE.}
  \item{mlt}{integer specifying the number of runs of the optimization algorithm in ML mode. Ignored if Bayesian mode is selected. Defaults to 10.}
  \item{it}{integer specifying the number of MCMC iterations in Bayesian mode. Default is 100000.}
  \item{bi}{integer specifying the number of iterations to discard as burn-in in Bayesian mode. Default is 5000.}
  \item{sa}{integer specifying the number of iterations to skip between samples in Bayesian mode. Default is 100.}
  \item{silent}{logical indicating whether to print BayesTraits output. Defaults to FALSE} 
    \item{rm}{logical indicating whether to clean up and remove all temporary files. Defaults to TRUE} 
}
\details{
	Runs BayesTraits continous models for one or more continuous traits. See BayesTraits documentation for more information: http://www.evolution.rdg.ac.uk/BayesTraits.html
	
	This function works by using `System` to write input files, run BayesTraits on your system, import the output into R, and finally delete both the input and output files from your system. When everything runs smoothly, there should be no extraneous files left on your system. However, if an error occurs, it is possible that some files will not be deleted. They will be overwritten the next time you run this function. 
	
	WARNING: The default settings for the number of iterations and burn-in period for Bayesian analysis is wayyyyy too low for a real analysis. 
}
\value{
	Results from BayesTraits output in a dataframe. 
}
\author{
	Randi H. Griffin 
}
\references{
	Pagel, M. (1997). Inferring evolutionary processes from phylogenies." Zoologica Scripta 26(4): 331-348. 

	Pagel, M. (1999). Inferring the historical patterns of biological evolution." Nature 401(6756): 877-884. 
}


