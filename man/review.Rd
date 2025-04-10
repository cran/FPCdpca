\name{review}
\alias{review}
\docType{data}
\title{Review}
\description{
This dataset contains travel reviews from TripAdvisor.com, covering destinations in 11 categories across East Asia. Each traveler's rating is mapped to a scale from Terrible (0) to Excellent (4), and the average rating for each category per user is provided.
}
\usage{
review
}
\format{
A data frame with multiple rows and 12 columns.
\itemize{
\item \code{User_ID}: Unique identifier for each user (Categorical).
\item \code{Art_Galleries}: Average user feedback on art galleries.
\item \code{Dance_Clubs}: Average user feedback on dance clubs.
\item \code{Juice_Bars}: Average user feedback on juice bars.
\item \code{Restaurants}: Average user feedback on restaurants.
\item \code{Museums}: Average user feedback on museums.
\item \code{Resorts}: Average user feedback on resorts.
\item \code{Parks_Picnic_Spots}: Average user feedback on parks and picnic spots.
\item \code{Beaches}: Average user feedback on beaches.
\item \code{Theaters}: Average user feedback on theaters.
\item \code{Religious_Institutions}: Average user feedback on religious institutions.
}
}
\details{
The dataset is populated by crawling TripAdvisor.com and includes reviews on destinations in 11 categories across East Asia. Each traveler's rating is mapped as follows: Excellent (4), Very Good (3), Average (2), Poor (1), and Terrible (0). The average rating for each category per user is used.
}
\note{
This dataset is licensed under a Creative Commons Attribution 4.0 International (CC BY 4.0) license, which allows for sharing and adaptation of the data for any purpose, provided that appropriate credit is given.
}
\source{
UCI Machine Learning Repositor}

\examples{
# Load the dataset
data(review)

# Print the first few rows of the dataset
print(head(review))

# Summary statistics for the ratings
summary(review[, 2:11])
}
\keyword{datasets}
