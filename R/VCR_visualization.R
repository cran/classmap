

confmat.vcr = function(vcrout, cutoff = 0.99,
                       showClassNumbers = FALSE,
                       showOutliers = TRUE, silent = FALSE) {
  # Builds a confusion matrix from a vcr object.
  #
  # Arguments:
  # vcrout           : output of vcr.*.train() or vcr.*.newdata().
  # cutoff           : the objects with overall farness
  #                    ofarness > cutoff are flagged as outliers.
  # showClassNumbers : if TRUE, the row and column names are the
  #                    number of each level instead of the level
  #                    itself. Useful for long level names.
  # showOutliers     : if TRUE and some points were flagged as
  #                    outliers, it adds an extra column on the
  #                    right of the confusion matrix for these
  #                    outliers, with label "outl".
  # silent           : if FALSE, the confusion matrix and
  #                    accuracy are shown on the screen.
  #
  # Returns:
  # confMat          : the confusion matrix.
  #
  levels = vcrout$levels
  nlab   = length(levels)
  if(is.null(vcrout$yintnew)){
    if(is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew") }
    given = vcrout$yint
    training = TRUE
  } else {
    given = vcrout$yintnew
    training = FALSE
  }
  whichyint = if(training) {"yint"} else {"yintnew"}
  indsv = which(!is.na(given))
  if(length(indsv) == 0){
    stop(paste0(whichyint," has no available values,",
                " so no confusion matrix can be made.")) }
  given = given[indsv]
  if(sum(!(given %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ",whichyint,"[indsv] lie in 1:",nlab)) }
  avobs = sort(unique(given)) # available given values
  predicted = vcrout$predint[indsv] # corresponding predictions
  if(sum(is.na(predicted)) > 0){
    stop("\n predint[indsv] has missing values.")}
  if(sum(!(predicted %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all predint[indsv] lie in 1:",nlab)) }
  # compute the accuracy before taking outliers into account:
  accuracy = sum(given == predicted)/length(indsv)
  myrownames = levels[avobs]
  ofarness = vcrout$ofarness[indsv]
  outliers = (ofarness > cutoff)
  if(showOutliers & sum(outliers) > 0) {
    predicted[which(outliers)] = nlab+1
    confMat = table(given, predicted)
    mycolints = sort(unique(predicted))
    mycolnames = c(levels,"outl")[mycolints]
    mycolints[which(mycolints == (nlab+1))] = "outl"
  } else {
    confMat = table(given, predicted)
    mycolints = sort(unique(predicted))
    mycolnames = levels[mycolints]
  }
  if(showClassNumbers){
    rownames(confMat) = avobs
    colnames(confMat) = mycolints
  } else {
    rownames(confMat) = myrownames
    colnames(confMat) = mycolnames
  }
  if(silent != TRUE){
    wnq(paste0("\nConfusion matrix:"))
    pnq(confMat)
    wnq(paste0("\nThe accuracy is ",round(100*accuracy,2),"%."))
  }
  invisible(confMat)
}



stackedplot = function(vcrout, cutoff = 0.99, classCols = NULL,
                       classLabels = NULL,
                       separSize=1, minSize=1.5,
                       showOutliers = TRUE, showLegend = FALSE,
                       htitle = NULL, vtitle = NULL) {
  # Vertically stacked mosaic plot of class predictions.
  #
  # Arguments:
  #   vcrout       : output of DA2vcr(), knn2vcr(), svm2vcr(),
  #                  neuralvcr(), treevcr(), or forestvcr().
  #   cutoff       : the cases with ofarness > cutoff will be
  #                  considered outliers.
  #   classCols    : user-specified colors for the classes.
  #                  If NULL a default palette is used.
  #   classLabels  : names of given labels. If NULL they are taken
  #                  from vcrout.
  #   separSize    : how much white between rectangles.
  #   minSize      : rectangles < minSize % are shown as minSize %
  #   showoUtliers : if TRUE, shows a separate class in gray with
  #                  the outliers, always at the top.
  #   htitle       : title for horizontal axis (given labels)
  #   vtitle       : title for vertical axis (predicted labels)
  #   showLegend   : if T, legend is shown to the right of the plot.
  #                  Default F, since the legend is not necessary
  #                  as the colors are already visible in the
  #                  bottom part of each stack.
  #
  # Returns:
  #   gg           : the plot object.
  #
  if(is.null(vcrout$yintnew)){
    if(is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew") }
    whichdata = "[training]"
  } else { whichdata = "[newdata]" }
  if(is.null(htitle)) htitle = paste0(whichdata," given class")
  if(is.null(vtitle)) vtitle = paste0(whichdata," predicted class")
  #
  nlab = length(vcrout$levels)
  if(is.null(classLabels)) { # no classLabels are given
    lvls = vcrout$levels
    if("outl" %in% lvls) stop(
      "Please rename the level currently called 'outl'.")
  } else { # classLabels are given
    if(!is.vector(classLabels)){
      stop("\n classLabels should be a vector") }
    lvls = classLabels
    if("outl" %in% lvls) stop(
      "Do not use 'outl' as a class label.")
    if(length(lvls) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ",nlab,".")) }
  }
  if(is.null(classCols)) {
    mosaicColors = rainbow(nlab)
  } else {
    if(!is.vector(classCols)){
      stop("\n classCols should be a vector") }
    if(length(classCols) < nlab){
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ",nlab,".")) }
    mosaicColors = classCols[seq_len(nlab)]
  }
  confusi  = confmat.vcr(vcrout, cutoff = cutoff,
                         showOutliers = showOutliers,
                         showClassNumbers = TRUE, silent = TRUE)
  pconfusi = confusi/rowSums(confusi) * 100
  # = confusion matrix in percentages
  ayint = as.numeric(rownames(confusi)) # available yint
  nlaba = length(ayint) # nlab of available yint
  if(nlaba < nlab){
    wnq(paste0("\nNot all classes occur in these data.",
               " The classes to plot are:"))
    pnq(ayint)
  }
  colints = colnames(confusi) # integers and possibly "outl"
  colints[which(colints == "outl")] = nlab+1 # possibly none
  for(j in seq_len(nlaba)) { # avoid plotting tiny frequencies
    confvec  = pconfusi[j, ] # confusion vector
    tooSmall = which((confvec > 0) & (confvec < minSize))
    if (length(tooSmall) > 0) {
      confvec[tooSmall] = minSize
      confvec = confvec/sum(confvec)  * 100
      pconfusi[j, ] = confvec
    }
  }
  hasOutliers = 0 # by default the class of outliers pi_0 is empty
  if(showOutliers){
    if("outl" %in% colnames(confusi)) {
      lvls = c(lvls, "outliers")
      hasOutliers = 1
      mosaicColors = c(mosaicColors[seq_len(nlab)], "gray20")
    }
  }
  # dataframe for the first available class:
  confRow = rep(0,length(lvls))
  colints = as.numeric(colints)
  confRow[colints] = pconfusi[1, ]
  df = data.frame(
    "g1_fill"  = factor(lvls[seq_len(nlab + hasOutliers)], lvls),
    "g1_count" = confRow # counts of the predictions
  )
  # dataframe for other classes:
  if(nlaba > 1){
    for (j in 2:nlaba) { # j=2
      # j is the row number of confusi
      g                = ayint[j] # actual class number
      confRow          = rep(0,length(lvls))
      colints          = as.numeric(colints)
      confRow[colints] = pconfusi[j, ]
      fillname         = paste0("g", j, "_fill")
      perm             = c(g, (seq_len(nlab + hasOutliers))[-g])
      # permutes the colors, putting the class itself first
      df[[fillname]]   = factor(lvls[perm], lvls[perm])
      # has added a new column to df
      countname        = paste0("g", j, "_count")
      df[[countname]]  = confRow[perm]
    }
  }
  widths = table(vcrout$yint)/length(vcrout$yint)
  # widths of the given rectangles
  pos = 0.5*(cumsum(widths)+cumsum(c(0,widths[-length(widths)])))
  # positions of the centers of the given rectangles
  # Build ggplot:
  gg = ggplot(df)
  for (i in seq_len(nlaba)) { # fills column by column
    fillname  = paste0("g", i, "_fill")
    countname = paste0("g", i, "_count")
    gg = gg + geom_bar(mapping = aes_string(x = pos[i],
                                            y = countname,
                                            fill = fillname),
                       width = widths[i], stat = "identity",
                       size = separSize, colour = "white",
                       show.legend = showLegend)
  }
  gg = gg + scale_x_continuous(htitle, breaks = pos,
                               labels = lvls[ayint],
                               expand = c(0.02, 0))
  # adds htitle (horizontal title) and its available levels
  gg = gg + scale_y_continuous(vtitle, trans = "reverse",
                               expand = c(0.02, 0))
  gg = gg + theme(axis.text.y = element_blank(), # was axis.text.x
                  axis.ticks.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank())

  # adds vtitle (vertical title) and reverses vertical direction
  if(showLegend) { gg = gg + labs(fill = "Classes") }
  gg = gg + theme(axis.ticks.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank())
  # removes axis tickmarks and numbers
  gg = gg + scale_fill_manual(values = mosaicColors)
  # replaces default colors by selected colors. This even
  # works when not all available mosaicColors occur here.
  return(gg)
}




classmap = function(vcrout, whichclass, classLabels = NULL,
                    classCols = NULL, main = NULL, cutoff = 0.99,
                    plotcutoff = TRUE, identify = FALSE, cex = 1,
                    cex.main=1.2, cex.lab=NULL, cex.axis=NULL,
                    opacity = 1, squareplot = TRUE,
                    maxprob = NULL, maxfactor=NULL) {
  #
  # Draws class maps to visualize classification results.
  #
  # Arguments:
  # vcrout       : the result of a function vcr.*.*()
  # whichclass   : the class to be displayed. Required. The colors
  #                of the points are those of their predicted
  #                labels. Required.
  # classLabels  : the labels (levels) of the classes. If NULL, they
  #                are taken from vcrout.
  # classCols    : a list of colors for the class labels 1,... There
  #                should be at least as many as there are levels.
  #                If NULL the classCols are taken as 2,...
  # main         : title for the plot.
  # cutoff       : objects with overall farness > cutoff will be
  #                flagged as outliers and plotted with a black
  #                border around them.
  # plotcutoff   : if TRUE, a dashed vertical line is drawn at
  #                the cutoff.
  # identify     : if TRUE, left-click on a point to get its number,
  #                then ESC to exit.
  # cex          : passed on to plot().
  # cex.main     : same, for title.
  # cex.lab      : same, for labels on horizontal and vertical axes.
  # cex.axis     : same, for axes.
  # opacity      : determines opacity of plotted dots.
  #                Value between 0 and 1, where 0 is transparent
  #                and 1 is opaque.
  # squareplot   : if TRUE, the horizontal and vertical axis will
  #                get the same length.
  # maxprob      : draws the farness axis at least upto probability
  #                maxprob.
  #                If NULL, the limits are obtained automatically.
  # maxfactor    : if not NULL, a number slightly higher than 1
  #                to increase the space at the right hand side of
  #                the plot, to make room for marking points.
  #
  # Returns:
  # coordinates  : a matrix with 2 columns containing the
  #                coordinates of the plotted points. The
  #                first coordinate is the quantile of the
  #                farness probability.

  # Auxiliary function:
  qfunc = function(probs){
    # quantile function of N(0,1) restricted to [0,4].
    if(min(probs) < 0) stop("There are negative probs")
    if(max(probs) > 1) stop("there are probs > 1")
    qnorm(probs*(pnorm(4)-0.5) + 0.5)
  }

  # The main function starts here:
  nlab = length(vcrout$levels)
  if(is.null(classLabels)) { # no classLabels are given
    levels = vcrout$levels
  } else { # classLabels are given
    if(!is.vector(classLabels)){
      stop("\n classlabels should be a vector") }
    levels = classLabels
    if(length(levels) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ",nlab,".")) }
  }
  if(is.null(classCols)) {
    classCols = seq_len(nlab) + 1 # red, green, blue,...
  } else {
    if(!is.vector(classCols)){
      stop("\n classCols should be a vector") }
    if(length(classCols) < nlab){
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ",nlab,".")) }
    classCols = classCols[seq_len(nlab)]
  }
  if(is.null(vcrout$yintnew)){
    if(is.null(vcrout$yint)) stop("vcrout$yint is missing")
    yint = vcrout$yint
    training = TRUE
  } else {
    yint = vcrout$yintnew
    training = FALSE
  }
  whichdata = if(training) {"[training]"} else {"[newdata]"}
  whichyint = if(training) {"yint"} else {"yintnew"}
  indsv = which(!is.na(yint))
  if(length(indsv) == 0){
    stop(paste0(whichyint," has no available values,",
                " so the plot cannot be made.")) }
  yint     = vcrout$yint[indsv]
  predint  = vcrout$predint[indsv]
  PAC      = vcrout$PAC[indsv]
  ofarness = vcrout$ofarness[indsv]
  outliers = (ofarness > cutoff)
  # In the classmap we will plot farness at quantile positions:
  qfarness = qfunc(vcrout$farness[indsv])
  coordinates = cbind(qfarness,PAC)
  rownames(coordinates) = indsv # case numbers in entire set.
  #
  fardashed = qfunc(cutoff)
  if(!is.null(maxprob)) maxfar = min(maxprob,1)
  if(is.null(maxfactor)) maxfactor = 1
  maxfar = max(max(qfarness),fardashed)
  if(sum(!(yint %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ",whichyint," lie in 1:",nlab)) }
  ayint = sort(unique(yint)) # available yint values
  #
  if(sum(!(predint %in% seq_len(nlab))) > 0) {
    stop(paste0(whichdata," Not all predint lie in 1:",nlab)) }
  #
  if(is.numeric(whichclass)){
    if(!(whichclass %in% seq_len(nlab))) stop(paste0(
      "whichclass must be in 1:",nlab,", the number of classes."))
  } else { #
    classnumber = match(whichclass, vcrout$levels)
    if(is.na(classnumber)) stop(paste0(
      'level "',whichclass,'" does not exist.'))
    whichclass = classnumber
  }
  if(!(whichclass %in% ayint)) stop(paste0(
    "Class number ",whichclass," with label ",
    levels[whichclass]," has no objects to visualize."))
  if(is.null(main)) {
    main = paste0(whichdata," predictions of class ",
                  levels[whichclass]) }
  clinds   = which(yint == whichclass) # class indices
  predint  = predint[clinds]
  casecols = classCols[predint]
  coordinates = coordinates[clinds,]
  outliers = which(outliers[clinds])
  #
  # Prepare for the plot:
  #
  PAClim = PACsolid = c(0,1)
  mypch = rep(19,length(predint)) # disk (16 is a bit smaller)
  mypch[outliers] = 21 # filled circle
  mycex = rep(cex,length(predint))
  mycex[outliers] = 1.1*cex
  mycol = mybg = adjustcolor(casecols, alpha.f = opacity)
  mycol[outliers] = "black"
  mycex.main = 1
  if(!is.null(cex.main)) { mycex.main = cex.main }
  mycex.lab = 1
  if(!is.null(cex.lab)) { mycex.lab = cex.lab }
  mycex.axis = 1
  if(!is.null(cex.axis)) { mycex.axis = cex.axis }
  #
  # Make actual plot
  #
  if(squareplot == TRUE) par(pty="s")
  # makes the axes of the plot equally long
  farlim = c(0,maxfar*maxfactor)
  plot(0, 0, type= "n", ann = FALSE, axes = FALSE,
       xlim = farlim, ylim = PAClim)
  # makes an empty plot
  title(main = main, line = 1, cex.main = mycex.main)
  title(ylab =  "P[alternative class]", line = 2.3,
        cex.lab = mycex.lab)
  title(xlab = "farness from given class",
        line = 2.3, cex.lab = mycex.lab)
  rect(farlim[1]-farlim[2],-2,farlim[2]*2,0.5,col="lightgray",
       border=FALSE) # adds rectangle
  # col can be modified to other gray intensities
  par(new=TRUE) # overlay next plot on this one
  plot(coordinates, xlab = "", ylab = "", xlim = farlim,
       ylim = PAClim, col = mycol, bg = mybg, pch = mypch,
       cex = mycex, cex.axis = mycex.axis, xaxt = "n")
  probs = c(0, 0.5, 0.75, 0.9, 0.99, 0.999, 1)
  axis(1, at = qfunc(probs), labels = probs,
       cex.axis = mycex.axis)
  if(plotcutoff) abline(v = fardashed, lty=4)
  abline(h = PACsolid)
  attr(coordinates, "ids") = integer(0)
  if(identify) {
    message("Press the escape key to stop identifying.")
    iout = identify(coordinates, order = TRUE)
    ids = iout$ind[order(iout$order)]
    if(length(ids) > 0) {
      wnq("Identified point(s): ")
      pnq(ids)
      attr(coordinates, "ids") = ids
    } else { wnq("No points were identified.") }
  }
  invisible(coordinates)
}



silplot = function(vcrout, classLabels = NULL, classCols = NULL,
                   showLegend = TRUE, showClassNumbers = FALSE,
                   showCases = FALSE, drawLineAtAverage = FALSE,
                   topdown = TRUE, main = NULL, summary = TRUE){
  #
  # Draws probability-based silhouette plot of a classification.
  #
  # Arguments:
  # vcrout            : output of vcr.*.*
  # classCols         : user-specified colors for the classes.
  #                     If NULL a default palette is used.
  # classLabels       : names of given labels. If NULL they are the
  #                     levels taken from vcrout.
  #                     This option can be useful if we want to show
  #                     shorter labels in the plot.
  # showLegend        : if T, legend is shown to the right of the plot.
  # showClassNumbers  : if T, the legend will show the class numbers
  #                     instead of the class labels.
  # showcases         : if T, the plot shows the numbers of the cases.
  #                     They are only readable when the number of cases
  #                     is relatively small.
  # topdown           : if TRUE (the default), the silhouettes are
  #                     plotted from top to bottom. Otherwise they
  #                     are plotted from left to right.
  # drawLineAtAverage : if TRUE, drwas a line at the average value
  #                     of the s(i).
  # main              : title for the plot
  # summary           : if TRUE, puts a ummary table on the screen
  #                     with for each class its number, label, size,
  #                     and the average of its s(i).
  #
  # Returns:
  # gg                : the plot object.
  #
  nlab = length(vcrout$levels)
  if(is.null(classLabels)) { # no classLabels are given
    lvls = vcrout$levels
  } else { # classLabels are given
    if(!is.vector(classLabels)){
      stop("\n classLabels should be a vector") }
    lvls = classLabels
    if(length(lvls) != nlab) {
      stop(paste0("\n The number of classLabels should equal",
                  " length(vcrout$levels) = ",nlab,".")) }
  }
  if(is.null(classCols)) {
    classCols = rainbow(nlab)
  } else {
    if(!is.vector(classCols)){
      stop("\n classCols should be a vector") }
    if(length(classCols) < nlab){
      stop(paste0("\n The number of classCols should be at",
                  " least length(vcrout$levels) = ",nlab,".")) }
    classCols = classCols[seq_len(nlab)]
  }
  if(is.null(vcrout$yintnew)){
    if(is.null(vcrout$yint)) {
      stop("there is no vcrout$yint or vcrout$yintnew") }
    yintv = vcrout$yint
    training = TRUE
  } else {
    yintv = vcrout$yintnew
    training = FALSE
  }
  whichdata = if(training) {"[training]"} else {"[newdata]"}
  whichyint = if(training) {"yint"} else {"yintnew"}
  indsv = which(!is.na(yintv))
  if(length(indsv) == 0){
    stop(paste0(whichyint," has no available values,",
                " so no silhouette plot can be made.")) }
  if(length(indsv) < 2){
    stop(paste0("At least 2 cases with non-missing ",whichyint,
                " are required.")) }
  yintv = yintv[indsv]
  ayint = sort(unique(yintv))
  if(sum(!(yintv %in% seq_len(nlab))) > 0) {
    stop(paste0("Not all ",whichyint,"[indsv] lie in 1:",nlab)) }
  PAC = vcrout$PAC[indsv]
  if(sum(is.na(PAC)) > 0) stop("PAC[indsv] has missing values.")
  si = 1 - 2*PAC # = s(i) based on PAC
  if(is.null(main)) {
    main = paste0(whichdata, " Silhouette plot of classification") }
  df = as.data.frame(cbind(class=yintv, si, name=indsv))
  # head(df)
  if(topdown) { df = df[order(-df$class, df$si), ]
  } else { df = df[order(df$class, -df$si), ] }
  # sorts first by class number, then by si
  avswidth = mean(df$si) # overall average silhouette width
  df$name  = factor(df$name, levels = df$name)
  # df$class = as.factor(lvls[df$class])
  df$class = as.factor(df$class)
  # makes class a factor
  if(showClassNumbers){
    mapping = aes_string(x = "name", y = "si",
                         color = "class", fill = "class")
  } else {
    df$label = factor(lvls[df$class], levels=lvls[seq_len(nlab)])
    mapping  = aes_string(x = "name", y = "si",
                          color = "label", fill = "label")
  }
  gg = ggplot(df, mapping) +
    geom_bar(stat = "identity", show.legend = showLegend,
             size = 0.05, width = 0.75) +
    labs(y = paste0("Silhouette width s(i)"),
         caption = paste0("\nOverall average silhouette width: ",
                          round(avswidth,2)), x = "",
         title = paste0(main))
  gg = gg + scale_fill_manual(values = classCols[ayint],
                              aesthetics = c("colour", "fill"))
  if(topdown) gg = gg + coord_flip()
  # switches the horizontal and vertical axes
  if(drawLineAtAverage){
    gg = gg + geom_hline(yintercept = avswidth,
                         linetype = "dashed", color = "red")
  }
  if(!showCases){
    if(topdown){
      gg = gg + theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())
    } else {
      gg = gg + theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank()) }
  } else {
    if(topdown){
      gg = gg + theme(axis.text.y = element_text(angle = 0))
    } else {
      gg = gg + theme(axis.text.x = element_text(angle = 90)) }
  }
  gg = gg + theme(plot.title = element_text(hjust = 0.5)) # center title
  gg = gg + theme(axis.line.x = element_line(size = 0.25))
  gg = gg + scale_x_discrete(expand=c(0,0.013*length(indsv)))
  gg = gg + scale_y_continuous(expand=c(0,0), limits = c(-1, 1))
  gg = gg +  theme(panel.background = element_blank())
  silwidths_perclass = round(sapply(unique(df$class), function(y) mean(df$si[which(df$class == y)])), 2)

  if (topdown) {
    xlocs <- cumsum(c(0, rev(table(df$class)[-1]))) + rev(table(df$class))*5/6
  } else {
    xlocs <- rev(cumsum(c(0, rev(table(df$class)[-1]))) + rev(table(df$class))*1/2) # 2/6)
  }
  for (i in 1:length(silwidths_perclass)) {
    templabel = paste('paste(bar(s), " = ",', silwidths_perclass[i], ")")
    gg = gg + annotate("text", y = -0.5, x = xlocs[i],
                       label = templabel, parse = TRUE,
                       color = rev(classCols[ayint])[i])
  }
  gg = gg + theme(plot.title.position = "plot") # center title over whole figure
  gg = gg + theme(axis.title=element_text(size=10))
  gg = gg + theme(plot.caption = element_text(size=10,
                                              hjust = 0,
                                              margin=ggplot2::margin(0,0,0,0)))
  if(summary){
    ave = tapply(df$si, df$class, mean) # average s(i) by class
    n = tapply(df$class, df$class, length) # cardinality of each class
    sil.sum = data.frame(classNumber = names(ave),
                         classLabel = lvls[as.numeric(names(ave))],
                         classSize = n,
                         classAveSi = round(ave, 2),
                         stringsAsFactors = TRUE)
    print(sil.sum, row.names=FALSE)
  }
  return(gg)
}

