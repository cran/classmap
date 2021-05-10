\name{confmat.vcr}
\alias{confmat.vcr}
\title{
Build a confusion matrix from the output of a function \code{vcr.*.*}.
}

\description{
Build a confusion matrix from the output of a function \code{vcr.*.*}.
Optionally, a separate column for outliers can be added to the confusion matrix.
}

\usage{
confmat.vcr(vcrout, cutoff = 0.99, showClassNumbers = FALSE,
            showOutliers = TRUE, silent = FALSE)
}

\arguments{
 \item{vcrout}{output of \code{vcr.*.train} or \code{vcr.*.newdata}.}
 \item{cutoff}{cases with overall farness \code{vcrout$ofarness} > \code{cutoff} are flagged as outliers.}
 \item{showClassNumbers}{if \code{TRUE}, the row and column names are the number of each level instead of the level itself. Useful for long level names.}
 \item{showOutliers}{if \code{TRUE} and some points were flagged as outliers, it adds an extra column on the right of the confusion matrix for these outliers, with label "outl".}
 \item{silent}{if \code{FALSE}, the confusion matrix and accuracy are shown on the screen.}
}

\value{
A confusion matrix
}

\references{
Raymaekers J., Rousseeuw P.J., Hubert M. (2021). Class maps for visualizing
classification results. \emph{Technometrics}, forthcoming. \href{https://arxiv.org/abs/2007.14495}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J.
}

\seealso{
\code{\link{vcr.da.train}}, \code{\link{vcr.da.newdata}}, \code{\link{vcr.knn.train}}, \code{\link{vcr.knn.newdata}}, \code{\link{vcr.svm.train}}, \cr
\code{\link{vcr.svm.newdata}}
}

\examples{

vcrout = vcr.knn.train(scale(iris[,1:4]), iris[,5], k=5)
# The usual confusion matrix:
confmat.vcr(vcrout, showOutliers = FALSE)

# Cases with ofarness > cutoff are flagged as outliers:
confmat.vcr(vcrout, cutoff = 0.98)

# With the default cutoff = 0.99 only one case is flagged here:
confmat.vcr(vcrout)
# Note that the accuracy is computed before any cases
# are flagged, so it does not depend on the cutoff.

confmat.vcr(vcrout,showClassNumbers=TRUE)
# Shows class numbers instead of labels. This option can
# be useful for long level names.

# For more examples, we refer to the vignettes:
vignette("Discriminant_analysis_examples")
vignette("K_nearest_neighbors_examples")
vignette("Support_vector_machine_examples")
}
