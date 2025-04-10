\name{riboflavinv100}
\alias{riboflavinv100}
\docType{data}
\title{Riboflavin Production Data (Top 100 Genes)}
\description{
This dataset is a subset of the riboflavin production data by Bacillus subtilis, containing \eqn{n = 71} observations. It includes the response variable (log-transformed riboflavin production rate) and the 100 genes with the largest empirical variances from the original dataset.
}
\usage{
data(riboflavinv100)
}
\format{
\describe{
\item{y}{Log-transformed riboflavin production rate (original name: \code{q_RIBFLV}). This is a continuous variable indicating the efficiency of riboflavin production by the bacterial strain.
}
\item{x}{A matrix of dimension \eqn{71 \times 100} containing the logarithm of the expression levels of the 100 genes with the largest empirical variances.
}
}
}
\details{
This dataset is derived from the original riboflavin dataset, which contains 4088 gene expressions. The riboflavinV100 dataset is created for ease of reproduction in examples and contains only the 100 genes with the largest empirical variances. It is commonly used in statistical research for high-dimensional data analysis.
}
\note{
The dataset is provided by DSM Nutritional Products Ltd., a leading company in the field of nutritional ingredients. The data have been preprocessed and normalized.
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
# Load the riboflavinv100 dataset
data(riboflavinv100)

# Display the dimensions of the dataset
print(dim(riboflavinv100$x))
print(length(riboflavinv100$y))

# Summary statistics for the response variable
summary(riboflavinv100$y)
}
\keyword{datasets}
