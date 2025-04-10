\name{riboflavin}
\alias{riboflavin}
\docType{data}
\title{Riboflavin Production Data}
\description{
This dataset contains measurements of riboflavin (vitamin B2) production by Bacillus subtilis, a Gram-positive bacterium commonly used in industrial fermentation processes. The dataset includes
\eqn{n = 71} observations with \eqn{p = 4088} predictors, representing the logarithm of the expression levels of 4088 genes. The response variable is the log-transformed riboflavin production rate.
}
\usage{
data(riboflavin)
}
\format{
\describe{
\item{y}{Log-transformed riboflavin production rate (original name: \code{q_RIBFLV}). This is a continuous variable indicating the efficiency of riboflavin production by the bacterial strain.
}
\item{x}{A matrix of dimension \eqn{71 \times 4088} containing the logarithm of the expression levels of 4088 genes. Each column corresponds to a gene, and each row corresponds to an observation (experimental condition or time point).
}
}
}
\details{
The riboflavin dataset is a high-dimensional dataset commonly used in statistical research, particularly in the fields of bioinformatics and systems biology. It was originally collected to study the genetic regulation of riboflavin biosynthesis in Bacillus subtilis. The data were generated using DNA microarray technology to measure gene expression levels under various experimental conditions.
}
\note{
The dataset is provided by DSM Nutritional Products Ltd., a leading company in the field of nutritional ingredients. The data have been preprocessed and normalized to account for technical variations in the microarray measurements.
}
\source{
DSM Nutritional Products Ltd., Basel, Switzerland.
}
\references{
\itemize{
\item Bühlmann, P., Kalisch, M., & Meier, L. (2014). 'High-dimensional statistics with a view towards applications in biology.' \emph{Annual Review of Statistics and its Applications}, \bold{1}, 255--278.
\item DSM Nutritional Products Ltd. (2005). 'Genome-scale analysis of Bacillus subtilis riboflavin production.' \emph{Internal Report}.
}
}
\examples{
# Load the riboflavin dataset
data(riboflavin)

# Display the dimensions of the dataset
print(dim(riboflavin$x))
print(length(riboflavin$y))

# Summary statistics for the response variable
summary(riboflavin$y)
}
\keyword{datasets}
