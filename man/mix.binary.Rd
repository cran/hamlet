\name{mix.binary}
\alias{mix.binary}
\title{
Binary coding of categorical variables
}
\description{
This function encodes categorical variables (e.g. columns of type 'factor' or 'character'). U new columns are created per each such column, where U is the number of unique instances of that column. The new columns are named OriginalColumnName_U1, OriginalColumnName_U2, etc.
}
\usage{
mix.binary(x)
}
\arguments{
  \item{x}{
A data.frame or a matrix where categorical columns are to be binary coded. Categorical columns are assumed to be all non-numeric fields.
}
}
\details{
A function that codes categorical variables in a dataset into binary variables. This is done in the following manner:
e.g. x = {red, green, blue, green} --> x_new = {{1,0,0}, {0,1,0}, {0,0,1}, {0,1,0}} where the dimensions in x_new are is_red, is_green and is_blue
}
\value{
The function returns a data.frame, where categorical variables have been replaced with 0/1-binary fields, and numeric fields have been left untouched. Notice that the order of the columns may not be the original.
}
\author{
Teemu Daniel Laajala <teelaa@utu.fi>
}
\examples{
data(vcapwide)

ex <- mix.binary(vcapwide[,c("Group", "CastrationDate")])
apply(ex, MARGIN=1, FUN=sum) 
# Notice that each row sums to 2, as two categorical variables were binary coded 
# and no missing values were present

mix.binary(vcapwide[,c("PSAWeek4", "Group", "CastrationDate")]) 
# Binary coding is only applied to non-numeric fields
}
\keyword{ manip }
