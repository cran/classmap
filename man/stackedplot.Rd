\name{stackedplot}
\alias{stackedplot}
\title{
Make a vertically stacked mosaic plot of class predictions.
}
\description{
Make a vertically stacked mosaic plot of class predictions from the output of
\code{vcr.*.train} or \code{vcr.*.newdata}. Optionally,
the outliers for each class can be shown as a gray rectangle at the top.
}
\usage{
stackedplot(vcrout, cutoff = 0.99, classCols = NULL, classLabels = NULL,
            separSize=1, minSize=1.5, showOutliers = TRUE, showLegend = FALSE,
            htitle = NULL, vtitle = NULL)
}

\arguments{
 \item{vcrout}{output of \code{vcr.*.train} or \code{vcr.*.newdata}.}
 \item{cutoff}{cases with overall farness \code{vcrout$ofarness} > \code{cutoff} are flagged as outliers.}
\item{classCols}{user-specified colors for the classes. If \code{NULL} a default palette is used.}
 \item{classLabels}{names of given labels. If \code{NULL} they are taken from \code{vcrout}.}
 \item{separSize}{how much white between rectangles.}
 \item{minSize}{rectangles describing less than \code{minSize} percent of the data, are shown as \code{minSize} percent.}
 \item{showOutliers}{if \code{TRUE}, shows a separate class in gray with the outliers, always at the top.}
 \item{showLegend}{if \code{TRUE}, legend is shown to the right of the plot. Default \code{FALSE}, since the legend is not necessary as the colors are already visible in the bottom part of each stack.}
 \item{htitle}{title for horizontal axis (given labels). If  \code{NULL}, a default title is shown.}
 \item{vtitle}{title for vertical axis (predicted labels). If  \code{NULL}, a default title is shown.}
}

\value{
A ggplot object.
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

data("data_floralbuds")
X = data_floralbuds[,1:6]; y = data_floralbuds[,7]
vcrout = vcr.da.train(X,y)
cols = c("saddlebrown","orange","olivedrab4","royalblue3")
stackedplot(vcrout, classCols = cols, showLegend=TRUE)

# The legend is not really needed, since we can read the
# color of a class from the bottom of its vertical bar:
stackedplot(vcrout, classCols = cols)

# If we do not wish to show outliers:
stackedplot(vcrout, classCols = cols, showOutliers=FALSE)

# For more examples, we refer to the vignettes:
vignette("Discriminant_analysis_examples")
vignette("K_nearest_neighbors_examples")
vignette("Support_vector_machine_examples")
}

