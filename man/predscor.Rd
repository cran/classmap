\name{predscor}
\alias{predscor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draws a predictions correlation plot, which visualizes the correlations between the prediction terms in a regression fit.
}
\description{
Computes the correlations between the prediction terms in a regression fit, and displays them graphically in a way that takes the standard deviations of the prediction terms into account. The input variables of the regression can be numerical, categorical, logical and character, and the regression model can be linear or generalized linear. The regression formula in \code{\link[stats:lm]{lm}} or \code{\link[stats:glm]{glm}} may contain transformations and interactions.
}
\usage{
predscor(fit, maxnpreds = 8, sort.by.stdev = TRUE, adj.order = FALSE,
cell.length = "stdev", plot.abs.cor = FALSE, palette = NULL,
diagonalcolor = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{fit}{an output object of \code{\link[stats:lm]{lm}} or \code{\link[stats:glm]{glm}}.}
 \item{maxnpreds}{the maximal number of prediction terms
 to plot. When there are more prediction
 terms than this, those with smallest
 standard deviations are combined.}
 \item{sort.by.stdev}{if \code{TRUE}, sorts the prediction terms by
 decreasing standard deviation.}
 \item{adj.order}{if \code{TRUE}, this modifies the order
 of the prediction terms in an attempt
 to bring highly correlated
 prediction terms close to each other
 in the \code{\link{predscor}} display. If
 \code{sort.by.stdev} is \code{TRUE}, this happens
 after that sorting.}
 \item{cell.length}{if "\code{stdev}", the sides of the square
 cells on the diagonal of the correlation
 matrix are proportional to the standard
 deviation of their prediction term.
 If "\code{sqrt}" they are proportional to the
 square root of the standard deviation.
 If "\code{equal}" all sides are the same.}
 \item{plot.abs.cor}{if \code{FALSE}, the default, positive and
 negative correlations are shown in
 different colors, typically red and
 blue. If \code{TRUE} the absolute values
 of the correlations are shown.}
 \item{palette}{a vector with colors to display
 correlations ranging from -1 to 1.
 If \code{NULL}, the default palette shows
 positive correlations in red,
 negative correlations in blue, and
 uses white for correlation zero.}
 \item{diagonalcolor}{color of the cells on the diagonal
 of the correlation matrix. The
 default is "black".}

}

\value{
A list containing
  \item{cormat}{the correlation matrix of the prediction terms.}
  \item{predterms}{matrix of cases by prediction terms.}
  \item{predsummary}{data frame with the standard deviation of each prediction term and the total linear prediction.}
}
\references{
Rousseeuw, P.J. (2025). Explainable Linear and Generalized Linear Models by the Predictions Plot
\url{https://arxiv.org/abs/2412.16980v2} (open access).
}
\author{
 Rousseeuw, P.J.
}



\seealso{
\code{\link{predsplot}}
}
\examples{
data(data_titanic)
attach(data_titanic)
Pclass = factor(Pclass, unique(Pclass))
Sex = factor(Sex, labels = c("F","M"))
fit <- glm(y ~ Sex + Age + SibSp + Parch + Pclass, family=binomial)
predscor(fit)

# For more examples, we refer to the vignette:
\dontrun{
vignette("predsplot_examples")
}
}
