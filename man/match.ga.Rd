\name{match.ga}
\alias{match.ga}
\title{
	Non-bipartite matching using the Genetic Algorithm (GA)
}
\description{
	An implementation of the Genetic Algorithm for solving non-bipartite matching tasks with customizable evolutionary events and parameters
}
\usage{
match.ga(d, g, 
	pops, 
	generations = 100, 
	popsize = 100, 
	nmutate = 100, 
	ndeath = 30, 
	type = "min", 
	mutate = hamlet:::.ga.mutate, 
	breed = hamlet:::.ga.breed, 
	weight = hamlet:::.ga.weight, 
	fitness = hamlet:::.ga.fitness, 
	step = hamlet:::.ga.step, 
	initialize = hamlet:::.ga.init, 
	progplot = T, 
	plot = T, 
	verb = 0, 
	progress = 500, 
	...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
A distance/dissimilarity matrix 'd'
}
  \item{g}{
The size in submatches, as in how many observations are always connected
}
  \item{pops}{
If user wants to specify starting populations, they can be provided here as a matrix. Each row correspondings to the observations, while columns are the different solutions (population in the GA). For example, a 10 row 100 column pops-matrix would be 100 different matching solutions of 10 observations. Each number in the matrix indicates a different submatch.
}
  \item{generations}{
Number of simulations to run in the GA. In each step, mutations, breeding and breeding occur according to user's specified settings, and a new generation is created out of this.
}
  \item{popsize}{
Number of solutions (='individuals') to have in each step of the algorithm.
}
  \item{nmutate}{
Number of mutations to occur in each step. Individuals are sampled with replacement, and then given the corresponding number of mutations.
}
  \item{ndeath}{
Number of deaths to occur in each step. Each dead solution (='individual') is then replaced by breeding suitable parents (probability of being a parent weighted by fitness).
}
  \item{type}{
Type of optimization, can be 'min' or 'max'.
}
  \item{mutate}{
Mutation function; by default the hamlet internal function '.ga.mutate' is used. This function takes in solution vector 'x'. Two random positions are then swapped, which could be seen as a form of a point mutation.
}
  \item{breed}{
Breeding function; by default the hamlet internal function '.ga.breed' is used. This function takes in solution vectors 'x' and 'y' ,which will be the parents, and the distance matrix 'd'. The products x*d and y*d are computed, and row-wise differences are computed between the two matrices. The row with the highest difference indicates where one of the parents can be most improved, and this trait is inherited from the other parent.
}
  \item{weight}{
Weighting function; by default the hamlet internal function '.ga.weight' is used. This weight should be correspond to probabilities that the corresponding individuals will undergo some sort of event (i.e. mutation, death) or participate in producing offspring (i.e. breed). This probability weight is computed according to ranks of fitnesses computed in the
}
  \item{fitness}{

Fitness function; by default the hamlet internal function '.ga.fitness' is used. This should yield the numeric fitness for a solution, indicating how viable the solution is in relation to the others. In a minimization task the lower fitness indicates better viability.
}
  \item{step}{
A step function; by default the hamlet internal function '.ga.step' is used. The step function which combines all operations in the GA, in order to produce the next generation of solutions given the previous one.
}
  \item{initialize}{
Initialization function; by default the hamlet internal function '.ga.initialize' is used. This function should format a set of valid solutions to produce the first generation in the beginning of the GA.
}
  \item{progplot}{
Should progress be plotted. If true, in every generation index dividable by the parameter 'progress', a function of fitnesses over the generations is plotted. The plot shows development of solution cost quantiles over time.
}
  \item{plot}{
Should the function plot the final quantiles over all the generations.
}
  \item{verb}{
Level of verbosity; -1 indicates omitting of verbal output, 0 indicates normal level, and +1 indicates debugging/additional information.
}
  \item{progress}{
How often should the function plot and print intermediate information on the progress.
}
  \item{\dots}{
Additional parameters for the internal GA functions.
}
}
\details{
The Genetic Algorithm (GA) is a form of an evolutionary optimization algorithm, where a population (a group of solutions to an optimization tasks) reproduce among themselves, die, mutate, and live on in a simulated environment. As the GA is a generic framework of solution approaches, it has many adjustable parameters and user may wish to explore many different options for the populations (for example in population size, mutation frequencies, fitness functions, drift etc) and also the evolutionary mechanics (such as breeding technique, types of mutations, and suitability for reproducing). Here, general default options and mechanics are provided, but it is advisable to explore different parameters for the particular optimization task in hand to find optimal solutions. If the user wishes to explore the implementation of the default mechanics, the function implementations are internally available in the hamlet-package. For example, the mutation function is accessible with the command: ' hamlet::.ga.mutate '.
}
\value{
The returned list compromises of:
	\itemize{
		\item A list of solutions; a matrix 'pops' which contains the population of solutions in the final generation of the algorithm, a vector 'fitnesses' which portrait the corresponding fitnesses to the columns of 'pops', and 'weights' which were the corresponding probabilities to events in the GA.
		\item A vector 'bestsol', for which the fitness function obtained minimum (or maximum) value during the algorithm.
		\item A value 'best', which is the optimum solution cost value observed during the algorithm.
	}
}
\note{
	Notice that end quality of the matching based allocation is heavily dependent on providing a feasible matrix D. One should carefully consider choice and tuning of the similarity metric. For example, Euclidean distance without standardization is often not a good choice as it does not normalize the variance of each variable and emphasis is on baseline variables that have a large relative variance.
	
	Note that the R-package 'GA' offers a wide range of generalized GA-related tools.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{match.bb}}
}
\examples{
# Set up a distance matrix and add dummies, then run GA
data(vcapwide)

# Construct an Euclidean distance example distance matrix using 15 observations from the VCaP study
d <- as.matrix(dist(vcapwide[1:15,c("PSAWeek10", "BWWeek10")]))
# Or rather, z-score transform all input variables first
d2 <- as.matrix(dist(scale(vcapwide[1:15,c("PSAWeek10", "BWWeek10")])))

# Notice that random simulations take place, so we will fix the RNG seed for reproducibility
set.seed(1)
# Resulting genetic algorithm progression is plotted by default
ga <- match.ga(d2, g=3, generations=60)
str(ga)
# Submatches, i.e. similar individuals that ought to be allocated to separate groups
ga[[2]]

}


\keyword{ ga }
