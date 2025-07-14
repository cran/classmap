\name{predsplot}
\alias{predsplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a predictions plot
}
\description{
Plots the prediction terms of a regression, together with the total prediction.  The input variables of the regression can be numerical, categorical, logical and character, and are visualized by histograms, densities, or bar plots. The regression model can be linear or generalized linear. The regression formula in lm() or glm() may contain transformations and interactions.  The plot shows all cases in the data, but can also be made for a single case, either in-sample or out-of-sample, to explain its prediction.
}
\usage{
predsplot(fit, maxnpreds = 8, sort.by.stdev = TRUE, displaytype = "histogram",
totalpred.type = "response", trunc.totalpred = TRUE, casetoshow = NULL,
staircase = FALSE, verbose = TRUE, maxchar.level = 5, nfact = 8, main = NULL,
cex.main = 1, xlab = NULL, ylab = NULL, vlinewidth = 0.25, hlinewidth = 0.25,
drawborder = TRUE, borderwidth = 1.0, densitycolors = NULL, predupcolor = "red",
preddowncolor = "blue", predpointsize = 3, draw.segments = TRUE,
predsegmentwidth = 0.8, profile = FALSE, bw = "nrd0", adjust = 1)
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
 \item{displaytype}{when \code{"histogram"}, the distributions
of the numerical prediction terms and
the total prediction are displayed as
histograms. When \code{"density"}, the
default density estimate of R is
plotted instead.}
 \item{totalpred.type}{when fit is a \code{glm} object, option
\code{"response"} plots labels of the total
prediction that are obtained by the
inverse link function. Option
\code{"linear"} plots labels of the total
linear prediction. When fit is an
\code{lm} object, argument \code{totalpred.type}
is ignored.}
 \item{trunc.totalpred}{if \code{TRUE}, the default, the range of
the total prediction is truncated
so ymin and ymax are determined by
the individual prediction terms. \code{FALSE}
may make the prediction terms look small
in the plot.}
 \item{casetoshow}{if not \code{NULL}, the particular case
to be displayed. This can be a case
number or row name of a case in the
dataset, or a list or vector with
input values of a new case.}
 \item{staircase}{if \code{TRUE} and casetoshow is not
\code{NULL}, the prediction for the case
is shown in staircase style.}
 \item{verbose}{if \code{TRUE}, some intermediate results are
shown on the console.}
 \item{maxchar.level}{only the first maxchar.level
characters of the categorical levels
are displayed.}
 \item{nfact}{if a numeric input variable has at most
\code{nfact} unique values, it will be
displayed as if it were a factor.
This may be useful since R considers
a binary variable to be numeric.}
 \item{main}{main title of the plot.}
 \item{cex.main}{its size.}
 \item{xlab}{horizontal legend. Its size can be changed by setting
the height and width of the plot.}
 \item{ylab}{vertical legend.
Its size can be changed by setting
the height and width of the plot.}
 \item{vlinewidth}{width of the vertical lines.}
 \item{hlinewidth}{width of the horizontal line.}
 \item{drawborder}{if \code{TRUE}, draws a box outside the entire
plot.}
\item{borderwidth}{width of the border. Defaults to 1.}
 \item{densitycolors}{a vector with 4 colors. The first is
for numeric input variables pointing up,
the second for numeric input variables
pointing down, the third for
prediction terms without orientation
and the total prediction, and the
fourth for factors. If \code{NULL}, the
default colors are used.}
 \item{predupcolor}{the color of a positive prediction
for casetoshow.}
 \item{preddowncolor}{the color of a negative prediction
for casetoshow.}
 \item{predpointsize}{the size of the points displaying the
predictions for casetoshow.}
 \item{draw.segments}{if \code{TRUE}, also plots a line segment
from the center of each prediction term
to the prediction term for casetoshow,
in the same color as the prediction.}
 \item{predsegmentwidth}{the width of that line segment.}
 \item{profile}{when casetoshow is not \code{NULL} and staircase
is \code{FALSE}, this plots the profile of the
case with a feint grey line.}
 \item{bw}{the bandwidth of the density estimation,
only used when \code{displaytype = "density"}.
This is the argument \code{'bw'} of the function
\code{\link[stats:density]{density}}.}
 \item{adjust}{multiplier of the bandwidth of the density
estimation, only used when \code{displaytype =
"density"}. This is the argument \code{'adjust'}
of the function \code{\link[stats:density]{density}}.}
}

\value{
A list with items
 \item{p}{the predictions plot, which is a ggplot2
 object.}
 \item{totpred}{vector with the total linear
 prediction of all cases for which
 it can be computed.}
 \item{centercept }{the centercept, which is the total
 linear prediction when all prediction
 terms have their average value.}
 \item{predterms}{matrix of cases by prediction terms.}
 \item{predsummary}{data frame with the standard deviation of
 each prediction term and the total linear
 prediction.}
 \item{casetotpred}{a number, the total linear prediction
 for casetoshow if one is given.}
 \item{casepredterms}{a vector with the values of the
 prediction terms for casetoshow.}
 \item{casesummary}{data frame which shows all prediction
 terms for casetoshow together with
 the centercept, total linear
 prediction, and for a glm fit also
 the total prediction in response units.}
}
\references{
Rousseeuw, P.J. (2025). Explainable Linear and Generalized Linear Models by the Predictions Plot
\url{https://arxiv.org/abs/2412.16980v2} (open access).
}
\author{
 Rousseeuw, P.J.
}

\seealso{
\code{\link{predscor}}
}
\examples{
data(data_titanic)
attach(data_titanic)
survival = y
Pclass = factor(Pclass, unique(Pclass))
Sex = factor(Sex, labels = c("F","M"))
fit <- glm(survival ~ Sex + Age + SibSp + Parch + Pclass, family=binomial)
predsplot(fit, main = "Titanic data", displaytype = "density")

# For more examples, we refer to the vignette:
\dontrun{
vignette("predsplot_examples")
}
}
