

predscomp <- function(fit, maxnpreds = 8,
                      sort.by.stdev = TRUE,
                      adj.order = FALSE,
                      casetoshow = NULL,
                      verbose = TRUE){
  # This is an auxiliary function which is called by
  # both predsplot() and predscor(). It does not
  # need to be visible to the casual user of the
  # package.
  #
  # Arguments:
  #
  # fit           : an output object of lm() or glm().
  # maxnpreds     : the maximal number of prediction terms to
  #                 plot. When there are more prediction terms
  #                 than this, those with smallest
  #                 standard deviations are combined.
  # sort.by.stdev : if TRUE, sorts the prediction terms by
  #                 decreasing standard deviation.
  # adj.order     : if TRUE, this modifies the order
  #                 of the prediction terms in an attempt
  #                 to bring highly correlated
  #                 prediction terms close to each other
  #                 in the predscor() display. If
  #                 sort.by.stdev is TRUE, this happens
  #                 after that sorting.
  # casetoshow    : if not NULL, the particular case
  #                 to be displayed by predsplot().
  #                 This can be a case number or
  #                 row name of a case in the fitted
  #                 dataset, or a list or vector with
  #                 input values of a new case.
  # verbose       : if TRUE, some computational results are
  #                 shown on the console.

  # Check whether fit is an lm or glm object:
  myclass = class(fit)
  mymodeltype = "lm"
  # check whether the model only contains "lm" and "glm"
  if(sum(!(myclass %in% c("glm","lm"))) > 0) stop(paste0(
    "Your fit object is of class ",myclass," but this",
    " tool requires the fit to be from lm() or glm(). \n"))
  #
  # Compute total prediction:
  #
  if("glm" %in% myclass) { # then is a glm
    totpred = predict(fit, type="link", na.action = na.omit)
    # this is the linear total prediction
    mymodeltype = "glm"
  } else {
    totpred = predict(fit, type="response", na.action = na.omit)
  }
  if(is.vector(totpred) == FALSE){
    if(is.matrix(totpred) == TRUE){
      stop(paste0("This model predicts a multivariate response.\n",
                  "Please rerun the fit for one response at a time."))
    } else {
      stop(paste0("Running predict(fit, type=\"response\"",
                  " did not return a vector."))
    }
  }
  centercept = mean(totpred, na.rm = TRUE)
  totpred = totpred - centercept
  totpredscale = sd(totpred, na.rm = TRUE)

  # Extract response name:
  attformul = attributes(terms.formula(fit$call$formula))
  yposit    = attformul$response
  yname     = rownames(attformul$factors)[yposit]

  # Compute prediction terms:
  predterms = predict(fit, type="terms", na.action = na.omit)
  npreds = ncol(predterms)
  colm = apply(predterms, 2, mean, na.rm=T)
  # We submean the prediction terms, which is
  # needed for regression without intercept:
  predterms = scale(predterms, center = colm, scale = FALSE)
  predscales = apply(predterms, 2, sd, na.rm=T)
  regressortypes = attributes(fit$terms)$dataClasses[-1]
  mycoeffs = fit$coefficients
  if(sum(is.na(mycoeffs)) > 0) stop(paste0(
    "This fit contains coefficients that are NA.\n",
    "  This implies that one or more coefficients are aliased",
    " in your model."))
  intercept = NULL
  if("(Intercept)" %in% names(mycoeffs)){
    jj = which(names(mycoeffs) == "(Intercept)")
    intercept = mycoeffs[jj]
  }

  combinedpreds = FALSE
  predtermnames = colnames(predterms)
  if(!(is.numeric(maxnpreds))) stop("maxnpreds must be numeric")
  if(length(maxnpreds)>1) stop("maxnpreds must have length 1")
  if(maxnpreds < 2) stop("maxnpreds must be at least 2")
  if(maxnpreds < npreds){ # for npreds <= 2 we indeed do nothing
    rankpreds = rank(predscales, ties.method = "last")
    tocombine = which(rankpreds < (npreds - maxnpreds + 2))
    remterm   = rowSums(predterms[,tocombine])
    tokeep    = which(rankpreds > (npreds - maxnpreds + 1))
    predterms = cbind(predterms[,tokeep],remterm)
    regressortypes = c(regressortypes[tokeep],"combination")
    colnames(predterms)[maxnpreds] = "remainder"
    predtermnames = colnames(predterms) # reordered
    predscales = apply(predterms, 2, sd, na.rm=T)
    npreds     = ncol(predterms)
    combinedpreds = TRUE
  }

  # determine the directions:
  directions = myslopes = rep(0, npreds)
  predtypes  = rep("none", npreds)
  names(predtypes) = predtermnames
  for (jj in seq_len(npreds)){ # jj=4
    # Note that interactions aren't in the regressortypes
    # list, only in the list of prediction terms. So for
    # each prediction term, we first check whether it is
    # in the regressortypes list.
    varname = predtermnames[jj]
    if(varname %in% names(regressortypes)){
      ll = which(names(regressortypes) == varname)
      if(regressortypes[ll] == "factor"){
        predtypes[jj] = "factor"
      }
      if(regressortypes[ll] == "numeric"){
        if(varname %in% names(mycoeffs)){
          predtypes[jj] = "numeric" # it was already
          kk = which(names(mycoeffs) == varname)
          myslopes[jj] = mycoeffs[kk]
          if(mycoeffs[kk] > 0) directions[jj] = 1
          if(mycoeffs[kk] < 0) directions[jj] = -1
          # any zero coeffs keep directions = 0
        }
      }
    }
  }
  names(predscales) = predtermnames
  names(directions) = predtermnames
  names(myslopes)   = predtermnames
  updown = rep("",npreds)
  for(jj in seq_len(npreds)){
    if(directions[jj] ==  1) updown[jj] = "up"
    if(directions[jj] == -1) updown[jj] = "down"
  }
  if(mymodeltype == "lm"){
    tpredname = paste0("Total prediction of ", yname)
  } else {
    tpredname = paste0("Total linear prediction of ", yname)
  }
  termnames = c(predtermnames, tpredname)
  sdevs = c(predscales, totpredscale)
  sdevs = signif(sdevs,4)
  updown = c(updown, "")
  predsummary = data.frame(termnames,
                           sdevs,
                           updown,
                           row.names = NULL)
  colnames(predsummary) = c("prediction term",
                            "stdev", "up/down")
  if(verbose == TRUE & is.null(casetoshow)){
    cat("\n")
    cat("Prediction terms in the model:\n")
    print(predsummary, row.names = FALSE)
    cat("\n")
  }

  ## Prediction for one case:
  ##
  showcasename = ""
  casesummary = NULL
  casetotpred = 0.0; casepredterms = predscales*0.0
  casenames = rownames(predterms)
  if(!is.null(casetoshow)){
    showcasename = " for the new case"
    if(!is.data.frame(casetoshow) & !is.list(casetoshow)){
      # then it is not a set of input variables
      if(length(casetoshow) != 1) stop(paste0(
        "casetoshow is not a list or dataframe of input variables,\n",
        "  a single case number, or a single case name."))
      if(casetoshow %in% casenames){
        # (This can be a case number or a name.)
        ii = which(casenames == casetoshow)
        showcasename = paste0(" for case ",casetoshow)
        casetotpred = totpred[ii]
        casepredterms = predterms[ii,,drop=F]
        # These prediction terms were already submeaned before.
        if(mymodeltype == "glm") {
          casetpred2 = predict(fit, type="response",
                               na.action = na.omit)[ii]
        }
      } else { stop(paste0(
        "casetoshow = ", casetoshow,
        " is not a valid case in the dataset."))
      }
    } else { # here it should be a set of input variables,
      # either in the form of a dataframe or a list.
      if(verbose == TRUE){
        cat("\n The case to show is: \n")
        print(casetoshow)
      }
      if(mymodeltype == "glm") {
        casetotpred = predict(fit, casetoshow, type="link",
                              na.action = na.omit)
        # = linear totpred
        casetpred2 = predict(fit, casetoshow, type="response",
                             na.action = na.omit)
        # = totpred in response units
      } else {
        casetotpred = predict(fit, casetoshow, type="response",
                              na.action = na.omit)
      }
      # This will correctly stop when some input variables
      # are missing that are needed in the model.
      casetotpred = casetotpred - centercept
      casepredterms = predict(fit, casetoshow, type="terms",
                              na.action = na.pass)
      # Here we do need to submean:
      casepredterms = casepredterms - colm
    }
    colnames(casepredterms) = colnames(predterms)
    listtermnames = c(predtermnames, "SUM", "centercept",
                      tpredname)
    listpreds1 = c(casepredterms, casetotpred)
    listpreds2 = c(centercept,casetotpred+centercept)
    if(mymodeltype == "glm") {
      listpreds2 = c(listpreds2, casetpred2)
      totprdname = paste0("Total prediction of ", yname,
                          " in response units")
      listtermnames = c(listtermnames, totprdname)
    }
    listpredterms = c(formatC(listpreds1, format="f", digits=5, flag="+"), formatC(listpreds2, format="f", digits=5))
    casesummary = data.frame(listtermnames,
                             listpredterms,
                             row.names = NULL)
    colnames(casesummary) = c("prediction term","value")
    # rownames(casesummary) = NULL
    if(verbose == TRUE){
      cat("\n ")
      cat(paste0("Prediction terms",showcasename,":\n"))
      print(casesummary, row.names = FALSE)
      cat("\n")
    }
  }

  ## Reordering the prediction terms:
  ##
  orderpreds = seq_len(npreds)
  if(sort.by.stdev == TRUE){
    if(combinedpreds == FALSE){
      orderpreds = order(predscales, decreasing = TRUE)
    } else { # we need to keep the remainder at the end:
      orderpreds = order(predscales[-npreds], decreasing = TRUE)
      orderpreds = c(orderpreds, npreds)
    }
  }
  regressortypes = regressortypes[orderpreds]
  predtypes      = predtypes[orderpreds]
  predscales     = predscales[orderpreds]
  predterms      = predterms[, orderpreds, drop=F]
  casepredterms  = casepredterms[orderpreds]
  directions     = directions[orderpreds]
  myslopes       = myslopes[orderpreds]

  dendro = NULL
  orderpreds = seq_len(npreds)
  if(adj.order == TRUE){
    if(combinedpreds == FALSE){
      if(npreds > 2){ # otherwise nothing to do
        cormat = cor(predterms)
        dendro <- as.dendrogram(hclust(as.dist(1 - abs(cormat)),
                                       method = "average"))
        # hclust does method = "complete" by default, but we
        # can also take "single", "average" (= UPGMA) etc.
        dendro <- reorder(dendro, seq_len(npreds))
        # Due to seq_len(npreds) this reorders the dendrogram
        # in a way that is closest to the previous order of
        # the variables.
        orderpreds <- order.dendrogram(dendro)
        if(length(orderpreds) != npreds)
          stop("dendrogram ordering gave wrong order length")
      }
    } else { # we need to keep the remainder last
      if(npreds > 3){ # otherwise nothing to do
        cormat = cor(predterms[-npreds])
        dendro <- as.dendrogram(hclust(as.dist(1 - abs(cormat)),
                                       method = "average"))
        dendro <- reorder(dendro, seq_len(npreds))
        orderpreds <- c(order.dendrogram(dendro), npreds)
        if(length(orderpreds) != npreds)
          stop("dendrogram ordering gave wrong order length")
      }
    }
  }
  regressortypes = regressortypes[orderpreds]
  predtypes      = predtypes[orderpreds]
  predscales     = predscales[orderpreds]
  predterms      = predterms[, orderpreds, drop=F]
  casepredterms  = casepredterms[orderpreds]
  directions     = directions[orderpreds]
  myslopes       = myslopes[orderpreds]
  attr(predterms, "names") = NULL

  outpt = list(mymodeltype = mymodeltype,
               yname = yname,
               totpred = totpred,
               centercept = centercept,
               intercept = intercept,
               totpredscale = totpredscale,
               predterms = predterms,
               myslopes = myslopes,
               predtypes = predtypes,
               regressortypes = regressortypes,
               predscales = predscales,
               directions = directions,
               predsummary = predsummary,
               showcasename = showcasename,
               casetotpred = casetotpred,
               casepredterms = casepredterms,
               casesummary = casesummary)
  return(outpt)
}


predsplot = function(fit,
                     maxnpreds = 8,
                     sort.by.stdev = TRUE,
                     displaytype = "histogram",
                     totalpred.type = "response",
                     trunc.totalpred = TRUE,
                     casetoshow = NULL,
                     staircase = FALSE,
                     verbose = TRUE,
                     maxchar.level = 5,
                     nfact = 8,
                     main = NULL,
                     cex.main = 1,
                     xlab = NULL,
                     ylab = NULL,
                     vlinewidth = 0.25,
                     hlinewidth = 0.25,
                     drawborder = TRUE,
                     borderwidth = 1.0,
                     densitycolors = NULL,
                     predupcolor = "red",
                     preddowncolor = "blue",
                     predpointsize = 3,
                     draw.segments = TRUE,
                     predsegmentwidth = 0.8,
                     profile = FALSE,
                     bw = "nrd0",
                     adjust = 1){

  # Arguments:
  #
  # fit              : an output object of lm() or glm().
  # maxnpreds        : the maximal number of prediction terms
  #                    to plot. When there are more prediction
  #                    terms than this, those with smallest
  #                    standard deviations are combined.
  # sort.by.stdev    : if TRUE, sorts the prediction terms by
  #                    decreasing standard deviation.
  # displaytype      : when "histogram", the distributions
  #                    of the numerical prediction terms and
  #                    the total prediction are displayed as
  #                    histograms. When "density", the
  #                    default density estimate of R is
  #                    plotted instead.
  # totalpred.type   : when fit is a glm object, option
  #                    "response" plots labels of the total
  #                    prediction that are obtained by the
  #                    inverse link function. Option
  #                    "linear" plots labels of the total
  #                    linear prediction. When fit is an
  #                    lm object, argument totalpred.type
  #                    is ignored.
  # trunc.totalpred  : if TRUE, the default, the range of
  #                    the total prediction is truncated
  #                    so ymin and ymax are determined by
  #                    the individual prediction terms. FALSE
  #                    may make the prediction terms look small
  #                    in the plot.
  # casetoshow       : if not NULL, the particular case
  #                    to be displayed. This can be a case
  #                    number or row name of a case in the
  #                    dataset, or a list or vector with
  #                    input values of a new case.
  # staircase        : if TRUE and casetoshow is not
  #                    NULL, the prediction for the case
  #                    is shown in staircase style.
  # verbose          : if TRUE, some intermediate results are
  #                    shown on the console.
  # maxchar.level    : only the first maxchar.level
  #                    characters of the categorical levels
  #                    are displayed.
  # nfact            : if a numeric input variable has at most
  #                    nfact unique values, it will be
  #                    displayed as if it were a factor.
  #                    This may be useful since R considers
  #                    a binary variable to be numeric.
  # main             : main title of the grill plot.
  # cex.main         : its size.
  # xlab, ylab       : horizontal and vertical legend.
  #                    Their size can be changed by setting
  #                    the height and width of the plot.
  # vlinewidth       : width of the vertical lines.
  # hlinewidth       : width of the horizontal line.
  # drawborder       : if TRUE, draws a box outside the entire
  #                    plot.
  # densitycolors    : a vector with 4 colors. The first is
  #                    for numeric input variables pointing up,
  #                    the second for numeric input variables
  #                    pointing down, the third for
  #                    prediction terms without orientation
  #                    and the total prediction, and the
  #                    fourth for factors. If NULL, the
  #                    default colors are used.
  # predupcolor      : the color of a positive prediction
  #                    for casetoshow.
  # preddowncolor    : the color of a negative prediction
  #                    for casetoshow.
  # predpointsize    : the size of the points displaying the
  #                    predictions for casetoshow.
  # draw.segments    : if TRUE, also plots a line segment
  #                    from the center of each prediction term
  #                    to the prediction term for casetoshow,
  #                    in the same color as the prediction.
  # predsegmentwidth : the width of that line segment.
  # profile          : when casetoshow is not NULL and staircase
  #                    is FALSE, this plots the profile of the
  #                    case with a feint grey line.
  # bw               : the bandwidth of the density estimation,
  #                    only used when displaytype = "density".
  #                    This is the argument 'bw' of the function
  #                    density().
  # adjust           : multiplier of the bandwidth of the density
  #                    estimation, only used when displaytype =
  #                    "density". This is the argument 'adjust'
  #                    of the function density().
  #
  # Values:
  #
  # p               : the predictions plot, which is a ggplot2
  #                   object.
  # totpred         : vector with the total linear
  #                   prediction of all cases for which
  #                   it can be computed.
  # centercept      : the centercept, which is the total
  #                   linear prediction when all prediction
  #                   terms have their average value.
  # predterms       : matrix of cases by prediction terms.
  # predsummary     : data frame with the standard deviation of
  #                   each prediction term and the total linear
  #                   prediction.
  # casetotpred     : a number, the total linear prediction
  #                   for casetoshow if one is given.
  # casepredterms   : a vector with the values of the
  #                   prediction terms for casetoshow.
  # casesummary     : data frame which shows all prediction
  #                   terms for casetoshow together with
  #                   the centercept, total linear
  #                   prediction, and for a glm fit also
  #                   the total prediction in response units.

  draw_densities <- function(plotdata, predterms, xcoord,
                             displaytype, height = 0.7,
                             bw = "nrd0", adjust = 1) {
    # Auxiliary function
    varnames <- plotdata$varnames
    nvar <- length(varnames)
    centers <- plotdata$centers
    scales <- plotdata$scales
    directions <- plotdata$directions
    distcolor <- plotdata$distcolor
    out <- as.list(rep(NA, nvar))
    for (jj in seq_len(nvar)) {
      miny = plotdata$yminvals[jj]
      maxy = plotdata$ymaxvals[jj]
      if(displaytype == "histogram" ||
         plotdata$showtypes[jj] %in% c("factor","num2fac")){
        histdata <- predterms[,jj] + centers[jj]
        #
        # If categorical variable and few levels,
        # draw a kind of barplot:
        upreds = unique(histdata)
        if(length(upreds) < 8){
          upreds = sort(upreds)
          delta = min(abs(diff(upreds)))/3
          delta = max(delta, (maxy-miny)/500)
          delta = min(delta, (maxy-miny)/200)
          breaks = sort(c(upreds - delta, upreds + delta))
          myhist = hist(histdata, plot = FALSE,
                        breaks = breaks)
          eps = 0
        } else {
          myhist = hist(histdata, plot = FALSE, breaks=50)
          eps    = (maxy-miny)/1000
        }
        breaks = myhist$breaks
        nbr    = length(breaks)
        allbreaks = c(breaks, breaks[2:(nbr-1)] - eps)
        allbreaks = sort(allbreaks)
        allbreaks = c(min(allbreaks)-eps, allbreaks,
                      max(allbreaks)+eps)
        mygrid = c(miny, allbreaks, maxy)
        densit = myhist$density
        maxd = max(densit)
        for(ii in seq_along(densit)){
          # Hack to make outliers visible:
          if(densit[ii] > 0) densit[ii] = max(densit[ii],maxd/6)
        }
        alldensity = rep(NA, 2*(nbr-1))
        for(ii in seq_along(densit)){
          alldensity[2*ii - 1] = densit[ii]
          alldensity[2*ii] = densit[ii]
        }
        curve = c(0, 0, alldensity, 0, 0)
        length(mygrid); length(curve)
      }
      if(displaytype == "density" &
         showtypes[jj] == "numeric"){
        predvals <- predterms[,jj] + centers[jj]
        dens = stats::density(predvals, from = miny, to = maxy,
                              bw = bw, adjust = adjust)
        mygrid = c(miny, dens$x, maxy)
        curve = c(0, dens$y, 0)
      }
      ngrid = length(mygrid)
      predterm = rep(varnames[jj], ngrid)
      denscolor = rep(distcolor[jj], ngrid)
      out[[jj]] <- data.frame(predterm = predterm,
                              denscolor = denscolor,
                              gridvals = mygrid,
                              curve = curve)
      cmax <- max(out[[jj]]$curve)
      out[[jj]]$curve <- out[[jj]]$curve*height/cmax + xcoord[jj]
    }
    out[[nvar]]$curve <- 2*xcoord[nvar] - out[[nvar]]$curve
    out <- do.call("rbind", out)
    return(out)
  }

  ## Here the main function starts
  #
  # Check whether fit$model is available:
  exist = !(is.null(fit$model))
  if(exist == FALSE) stop(paste0(
    "There is no $model in your fit object. Please rerun",
    "\n  the fit with the default option model = TRUE."))
  #
  if (!is.numeric(vlinewidth)) stop("vlinewidth must be a number.")
  if (!is.numeric(hlinewidth)) stop("hlinewidth must be a number.")
  if(!(displaytype) %in% c("histogram","density")) stop(
    "displaytype must be either \"histogram\" ot \"density\"")

  out = predscomp(fit=fit, maxnpreds = maxnpreds,
                  sort.by.stdev = sort.by.stdev,
                  casetoshow = casetoshow,
                  verbose = verbose)
  mymodeltype    = out$mymodeltype
  yname          = out$yname
  regressortypes = out$regressortypes
  predtypes      = out$predtypes
  totpred        = out$totpred
  centercept     = out$centercept
  predterms      = out$predterms
  scales         = out$predscales
  directions     = out$directions
  myslopes       = out$myslopes
  nvar           = length(scales)
  varnames       = names(scales)
  xcoord         = seq_len(nvar)
  xcasepos       = xcoord
  predsummary    = out$predsummary
  showcasename   = out$showcasename
  casepredterms  = out$casepredterms
  casesummary    = out$casesummary

  showregmat <- predterms
  showtypes = predtypes
  modmat <- model.matrix(fit) # these are all numbers
  # Obtain the original values of the input variables
  # (possibly after log or X^2 or x2*x3) before
  # they were multiplied by the coefficients:
  datset = fit$model # may contain character strings
  classes = sapply(datset, class)
  datacolnames = colnames(datset)
  valid.rows = rownames(predterms)
  datset = datset[valid.rows,,drop=F]
  truelevels = as.list(rep(NA, nvar))
  #
  for (jj in seq_len(nvar)){ # jj=4
    varname = varnames[jj]
    # Note that interactions aren't in the regressortypes
    # list, only in the list of prediction terms. So for
    # each prediction term, we first check whether it is
    # in the regressortypes list.
    if(varname %in% names(regressortypes)){
      ll = which(names(regressortypes) == varname)
      if(predtypes[ll] == "numeric"){
        # This excludes "none","remainder" as it should.
        if(varname %in% colnames(modmat)){
          kk = which(colnames(modmat) == varname)
          showregmat[,jj] = as.vector(modmat[,kk])
          truelevels[[jj]] = "none"
        }
      }
      if(predtypes[ll] != "numeric"){
        if(varname %in% datacolnames){
          kk = which(datacolnames == varname)
          if("factor" %in% classes[[kk]]) {
            # this includes ordinal factors!
            showtypes[jj] = "factor"
            showregmat[,jj]  = datset[[kk]]
            truelevels[[jj]] = levels(datset[[kk]])
          }
          if(classes[[kk]] == "character") {
            levs = unique(datset[[kk]])
            datset[[kk]] = factor(datset[[kk]], levels = levs)
            showtypes[jj] = "factor"
            showregmat[,jj]  = datset[[kk]]
            truelevels[[jj]] = levs
          }
          if(classes[[kk]] == "logical") {
            showtypes[jj] = "factor"
            showregmat[,jj][datset[[kk]] == FALSE] = 1
            showregmat[,jj][datset[[kk]] == TRUE] = 2
            truelevels[[jj]] = c("F","T")
          }
        }
      }
    } else {
      showtypes[jj] = "none"
      truelevels[[jj]] = "none"
    }
  }

  for (jj in seq_len(nvar)){ # jj=4
    # If a numeric prediction term takes fewer than nfact
    # values, we will plot it like a factor.
    if(showtypes[jj] %in% c("numeric","none")){
      # so not "factor"
      uniq = unique(as.vector(showregmat[,jj]))
      if(length(uniq) <= nfact){
        showtypes[jj] = "num2fac"
        truelevels[[jj]] = sort(uniq)
      }
    }
  }

  if(is.null(casetoshow) == TRUE) staircase = FALSE

  if(is.null(xlab)){
    xlab = "prediction terms and total prediction"
    xlab = paste0(xlab,showcasename)
  }

  if(is.null(ylab)){
    if(mymodeltype == "lm"){
      ylab = "contributions to the prediction"
      if(staircase == TRUE){
        ylab = "cumulative prediction" }
    } else {
      ylab = "contributions to the linear prediction"
      if(staircase == TRUE){
        ylab = "cumulative linear prediction" }
    }
  }
  ylab = paste0(ylab, " of ", yname)

  totalpredlabel = "linear"
  ymaxpreds = apply(predterms, 2, max, na.rm=T)
  yminpreds = apply(predterms, 2, min, na.rm=T)
  # Append total prediction below and subtract the
  # centercept (we give the right labels later).
  # Instead of its ymaxvals take the max of the
  # ymaxvals of the prediction terms (so when a case
  # is shown, we can get wider xlim and ylim.)
  nvar                 = nvar + 1
  varnames             = c(varnames, "totalpred")
  scales               = c(scales, out$totpredscale)
  predterms            = cbind(predterms, totpred)
  showtypes            = c(showtypes, "numeric")
  casepredterms        = c(casepredterms, out$casetotpred)
  directions           = c(directions, 0)
  names(scales)        = varnames
  names(predterms)     = varnames
  names(casepredterms) = varnames
  names(directions)    = varnames
  xcoord         = c(xcoord, nvar + 0.6)
  xcasepos       = xcoord
  xcasepos[nvar] = xcoord[nvar] - 0.008*nvar
  if(mymodeltype == "glm"){
    if(!(totalpred.type %in% c("linear","response"))){
      stop("totalpred.type must be \"linear\" or \"response\".")}
    totalpredlabel = totalpred.type
  }

  centers  = rep(0, nvar)
  ymaxvals = apply(predterms, 2, max, na.rm=T)
  yminvals = apply(predterms, 2, min, na.rm=T)

  if(is.null(densitycolors)){
    densitycolors = c("olivedrab3", "chocolate",
                      "grey60", "black")
  }
  distcolor = rep(densitycolors[3], nvar)
  for(jj in seq_len(nvar-1)){
    if(directions[jj] == 1) distcolor[jj] = densitycolors[1]
    if(directions[jj] == -1) distcolor[jj] = densitycolors[2]
    if(directions[jj] == 0) distcolor[jj] = densitycolors[3]
    if(showtypes[jj] == "factor") distcolor[jj] = densitycolors[4]
  }
  distcolor[nvar] = densitycolors[3]
  if(!(is.null(casetoshow) == TRUE)) distcolor = rep("grey60", nvar)

  casepos = as.vector(casepredterms) # all 0 if no casetoshow
  predpointcolor = rep(predupcolor,nvar)
  if(!is.null(preddowncolor)){
    for(j in seq_len(nvar)){
      if(casepos[j] < 0) predpointcolor[j] = preddowncolor
    }
  }

  if(staircase == TRUE){
    partsum = cumsum(c(0.0,casepos))
    centers = centers + partsum[-nvar]
    centers[nvar] = 0.0
    casepos   = casepos + centers
    yminvals  = yminvals + centers
    ymaxvals  = ymaxvals + centers
    yminpreds = yminpreds + centers[seq_along(yminpreds)]
    ymaxpreds = ymaxpreds + centers[seq_along(yminpreds)]
  }

  names(centers) = varnames
  names(casepos) = varnames
  if(trunc.totalpred == TRUE){
    globalmax = max(c(ymaxpreds,casepos))
    globalmin = min(c(yminpreds,casepos))
  } else {
    globalmax = max(c(ymaxvals,casepos))
    globalmin = min(c(yminvals,casepos))
  }
  globalran = globalmax - globalmin # range
  ymaxvals  = ymaxvals + 0.06*globalran
  yminvals  = yminvals - 0.06*globalran
  globalmax = globalmax + 0.06*globalran
  globalmin = globalmin - 0.25*globalran
  allymax = rep(globalmax, nvar)
  allymin = rep(globalmin, nvar)
  vlinewidth = rep(vlinewidth, nvar)
  vlinewidth[nvar] = 0.0
  eps = globalran/50
  arrowstart = arrowend = rep(0, nvar)
  for(jj in seq_len(nvar)){
    if(directions[jj] == 1){
      arrowstart[jj] = ymaxvals[jj] - eps
      arrowend[jj]   = ymaxvals[jj]
    }
    if(directions[jj] == -1){
      arrowstart[jj] = yminvals[jj] + eps
      arrowend[jj]   = yminvals[jj]
    }
  }

  plotdata = data.frame(varnames, centers, scales, showtypes,
                        yminvals, ymaxvals, directions,
                        casepos, xcoord, allymin, allymax,
                        vlinewidth, arrowstart, arrowend,
                        distcolor, predpointcolor,
                        row.names = NULL, check.rows = TRUE)
  height = 0.6
  densities <- draw_densities(plotdata = plotdata,
                              predterms = predterms,
                              xcoord = xcoord,
                              displaytype = displaytype,
                              height = height, bw = bw,
                              adjust = adjust)
  height = rep(height, nvar)
  height[nvar] = -height[nvar]
  plotdata$height = height

  p <- ggplot(data = plotdata,
              aes(x = xcoord,
                  y = centers,
                  ymin = yminvals,
                  ymax = ymaxvals)) + theme_classic() +
    theme(axis.line = element_line(color='white'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

  p <- p + geom_polygon(data = densities, aes(y = .data$gridvals, x = .data$curve, group = .data$predterm, fill = .data$denscolor), alpha = 0.7, show.legend = FALSE, inherit.aes = FALSE) + scale_fill_identity()

  ## Add vertical line segments:
  p = p + ggplot2::geom_segment(aes(x = xcoord, xend = xcoord, y = yminvals, yend = ymaxvals), colour = "black", linewidth = vlinewidth, show.legend = FALSE)

  ## Add arrowheads in the right directions.
  # (If directions==0 no arrowhead is drawn.)
  p = p + ggplot2::geom_segment(aes(x = xcoord, xend = xcoord, y = arrowstart, yend = arrowend), colour = "black", linewidth = vlinewidth, arrow = arrow(angle = 20, length = unit(abs(directions)*0.2, "cm"), ends = "last", type="open"), show.legend = FALSE)

  # Draw the horizontal line through zero:
  if(!(staircase == TRUE) | nvar == 2){
    p <- p + geom_hline(yintercept = 0, linetype = 1, linewidth = hlinewidth, color = "black") }

  # Draw the profile if requested:
  if(!(staircase == TRUE) & profile == TRUE & nvar > 2){
    p <- p + geom_segment(aes(x = c(xcoord[1],xcoord[1:nvar-1]), y = c(casepos[1],casepos[1:nvar-1]), xend = c(xcoord[2],xcoord[2:(nvar-1)],xcoord[nvar-1]), yend = c(casepos[2],casepos[2:(nvar-1)],casepos[nvar-1])), linewidth = 2*hlinewidth, colour = "grey65", linetype=1) }

  ## Draw horizontal line segments for prediction staircase:
  if(staircase == TRUE){
    # Plot a point at the center:
    p <- p + geom_point(aes(x = xcasepos[1:nvar], y = centers[1:nvar]), size = 1.5, colour = "black")
    # Plot the steps of the staircase:
    p <- p + geom_segment(aes(x = c(xcoord[1],xcoord[1:nvar-1]), y = c(casepos[1],casepos[1:nvar-1]), xend = c(xcoord[1],xcoord[2:nvar]), yend = c(casepos[1],casepos[1:nvar-1])), linewidth = hlinewidth, colour = "black", linetype=1) }

  ## Add colored points at predictions:
  if(!is.null(casetoshow)){
    bla = c("red","blue")
    if(draw.segments){
      p <- p + ggplot2::geom_pointrange(aes(x = xcasepos, y = casepos, ymin = pmin(casepos,centers), ymax = pmax(casepos,centers)), position = ggplot2::position_dodge(width = 0.5), fill = "white", fatten = predpointsize, colour = predpointcolor, linewidth = predsegmentwidth, show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_pointrange(aes(x = xcoord, y = casepos, ymin = casepos, ymax = casepos), position = ggplot2::position_dodge(width = 0.5), fill = "white", fatten = predpointsize, colour = predpointcolor, linewidth = 0, show.legend = FALSE)
    }
  }

  ## add vertical label
  p <- p + ylab(ylab)

  ## Add variable names and overall axis label
  p <- p + scale_x_continuous(breaks = plotdata$xcoord, labels = NULL, name = xlab)
  p <- p + xlab(xlab)

  ## Set the xlim and ylim:
  xleft = min(xcoord)-0.7; xright = max(xcoord)
  ybottom = 1.01*globalmin; ytop = 1.01*globalmax
  p <- p + coord_cartesian(xlim = c(xleft, xright),
                           ylim = c(ybottom,ytop),
                           expand = FALSE)

  ## Remove x-ticks
  p <- p + theme(axis.ticks.x = element_blank())

  ## draw the vertical axis and set size of break labels.
  myxrange = max(xcoord) - min(xcoord)
  mysize  = 2.4*min(c(21/(myxrange), 4.5))
  p <- p + theme(axis.line.y = element_line(linewidth = vlinewidth, colour = "black"), axis.text.y = element_text(size = mysize, color = "black") )

  ## Add second axis:
  if(totalpredlabel == "linear"){
    if(mymodeltype == "lm"){
      tpredname = "total prediction"
    } else {
      tpredname = "total linear prediction"
    }
    p <- p + scale_y_continuous(sec.axis = ggplot2::sec_axis(~. + centercept, name = paste0(tpredname," (including centercept = ",signif(centercept,3),")"), breaks = scales::breaks_extended(n = 9), guide = guide_axis(theme=theme(axis.line = element_line(linewidth = 0, colour = "white")))))
  }
  #
  if(totalpredlabel == "response"){
    # Now we have to transform the linear totpred
    # by the inverse of the link function.
    myident   = function(y) { y + centercept }
    myexpit   = function(y) { 1/(1 + exp(-(y + centercept))) }
    myinverse = function(y) { 1/(y + centercept) }
    myexp     = function(y) { exp(y + centercept) }
    my1osqr   = function(y) { 1/sqrt(y + centercept) }
    #
    mylink = fit$family$link
    if(mylink == "identity") mytransf = myident
    if(mylink == "logit")    mytransf = myexpit
    if(mylink == "inverse")  mytransf = myinverse
    if(mylink == "log")      mytransf = myexp
    if(mylink == "1/mu^2")   mytransf = my1osqr
    if(!(mylink %in% c("identity","logit","inverse",
                       "log","1/mu^2"))) stop(paste0(
                         "The link ",mylink," has not been implemented,",
                         "\n  only identity, logit, inverse, log, 1/mu^2."))
    cat(paste0("The glm fit used the ",mylink," link function. \n"))
    p <- p + scale_y_continuous(sec.axis = ggplot2::sec_axis(transform = mytransf, name = "total prediction in response units", breaks = scales::breaks_extended(n = 9), guide = guide_axis(theme=theme(axis.line = element_line(linewidth = 0, color ="white")),check.overlap = TRUE)))
  }
  # flip the right y-axis title to match the left y-axis:
  p <- p + theme(axis.title.y.right = element_text(angle=90))

  ## Annotations on the vertical lines:
  lennames = max(nchar(varnames))+2
  varsize = min(c(30/lennames, 3))
  for (jj in seq_len(ncol(showregmat))){
    varname = varnames[jj]
    p <- p + annotate("text", x = xcoord[jj], angle=90,
                      y = yminvals[jj],
                      label = paste0(varname,"  "),
                      hjust = 1, vjust = 0.3,
                      size = varsize)
    if(showtypes[jj] == "numeric"){
      showvals = as.numeric(showregmat[,jj])
      mybreaks = pretty(showvals, n=5,
                        min.n = 4, bounds = FALSE)
      rawmin = min(showvals); rawmax = max(showvals)
      tokeep = which(mybreaks >= rawmin & mybreaks <= rawmax)
      mybreaks = mybreaks[tokeep]
      mylabels = mybreaks
      for(pp in seq_along(mylabels)){
        mylabels[pp] = paste0(mylabels[pp],"-") }
      if(length(mybreaks) < 2) stop("under 2 breaks left")
      mymean  = mean(showvals, na.rm = TRUE)
      myslope = myslopes[jj]
      lchar   = max(nchar(mybreaks))
      mysize  = 0.9*min(c(21/(myxrange), 14/lchar, 4.5))
      ypos    = myslope*(mybreaks - mymean) + centers[jj]
      p <- p + annotate("text", x = xcoord[jj],
                        y = ypos, label = mylabels,
                        hjust = 1, size = mysize)
    }
    if(showtypes[jj] == "factor"){
      levls = truelevels[[jj]]
      unt = sort(unique(predterms[,jj]))
      mypreds = mylabs = rep(NA,length(levls))
      for(kk in seq_along(unt)){
        mm = which(predterms[,jj] == unt[kk])[1]
        numlevl = showregmat[mm, jj]
        mylabs[kk] = levls[numlevl]
        mylabs[kk] = strtrim(mylabs[kk], maxchar.level)
        mylabs[kk] = paste0(mylabs[kk],"-")
      }
      lchar = max(nchar(mylabs))
      mysize = 0.9*min(c(21/(myxrange), 14/lchar, 4.5))
      ypos = unt + centers[jj]
      p <- p + annotate("text", x = xcoord[jj],
                        y = ypos, label = mylabs,
                        hjust = 1, size = mysize)
    }
    if(showtypes[jj] == "num2fac"){
      levls = truelevels[[jj]]
      unt = rep(NA,length(levls))
      mypreds = mylabs = rep(NA,length(levls))
      for(kk in seq_along(levls)){
        mm = which(showregmat[,jj] == levls[kk])[1]
        unt[kk] = predterms[mm,jj]
        mylabs[kk] = levls[kk]
        mylabs[kk] = strtrim(mylabs[kk], maxchar.level)
        mylabs[kk] = paste0(mylabs[kk],"-")
      }
      lchar = max(nchar(mylabs))
      mysize = 0.9*min(c(21/(myxrange), 14/lchar, 4.5))
      ypos = unt + centers[jj]
      p <- p + annotate("text", x = xcoord[jj],
                        y = ypos, label = mylabs,
                        hjust = 1, size = mysize)
    }
  }

  if(!is.null(main)){
    p <- p + theme(plot.title = element_text(hjust = 0.5, size = rel(cex.main)))
    p <- p + ggtitle(main)
  }

  if(drawborder == TRUE){
    p <- p + theme(plot.background = element_rect(color = "black", fill = NA, linewidth = borderwidth))
  }
  plot(p)
  attr(predterms, "names") = NULL
  outp = list(p               = p,
              totpred         = totpred,
              centercept      = centercept,
              predterms       = predterms,
              predsummary     = out$predsummary,
              casetotpred     = out$casetotpred,
              casepredterms   = casepredterms,
              casesummary     = out$casesummary)
  invisible(outp)
}



predscor <- function(fit,
                     maxnpreds = 8,
                     sort.by.stdev = TRUE,
                     adj.order = FALSE,
                     cell.length = "stdev",
                     plot.abs.cor = FALSE,
                     palette = NULL,
                     diagonalcolor = "black"){
  #
  # The function predscor() computes the correlations
  # between the prediction terms in a regression fit, and
  # displays them graphically in a way that takes the
  # standard deviations of the prediction terms into account.
  #
  # Arguments:
  #
  # fit           : an output object of lm() or glm().
  # maxnpreds     : the maximal number of prediction terms
  #                 to plot. When there are more prediction
  #                 terms than this, those with smallest
  #                 standard deviations are combined.
  # sort.by.stdev : if TRUE, sorts the prediction terms by
  #                 decreasing standard deviation.
  # adj.order     : if TRUE, this modifies the order
  #                 of the prediction terms in an attempt
  #                 to bring highly correlated
  #                 prediction terms close to each other
  #                 in the predscor() display. If
  #                 sort.by.stdev is TRUE, this happens
  #                 after that sorting.
  # cell.length   : if "stdev", the sides of the square
  #                 cells on the diagonal of the correlation
  #                 matrix are proportional to the standard
  #                 deviation of their prediction term.
  #                 If "sqrt" they are proportional to the
  #                 square root of the standard deviation.
  #                 If "equal" all sides are the same.
  # plot.abs.cor  : if FALSE, the default, positive and
  #                 negative correlations are shown in
  #                 different colors, typically red and
  #                 blue. If TRUE the absolute values
  #                 of the correlations are shown.
  # palette       : a vector with colors to display
  #                 correlations ranging from -1 to 1.
  #                 If NULL, the default palette shows
  #                 positive correlations in red,
  #                 negative correlations in blue, and
  #                 uses white for correlation zero.
  # diagonalcolor : color of the cells on the diagonal
  #                 of the correlation matrix. The
  #                 default is "black".
  #
  # Values:
  #
  # cormat      : the correlation matrix of the
  #               prediction terms.
  # predterms   : matrix of cases by prediction terms.
  # predsummary : data frame with the standard deviation of each
  #               prediction term and the total linear
  #               prediction.

  if(!(cell.length) %in% c("stdev","sqrt","equal")) stop(
    "cell.length must be \"stdev\", \"sqrt\" or \"equal\"")
  out = predscomp(fit = fit,
                  maxnpreds = maxnpreds,
                  sort.by.stdev = sort.by.stdev,
                  adj.order = adj.order,
                  casetoshow = NULL,
                  verbose = FALSE)
  predterms = out$predterms
  scales    = out$predscales
  nvar      = length(scales)
  varnames  = names(scales)
  cormat    = cor(predterms)
  cat("Correlation matrix of the prediction terms:\n")
  print(round(cormat,2))
  if(nvar < 2){ cat(paste0(
    "\nNo correlation plot is drawn because there is only",
    " one prediction term."))
  } else {
    if(is.null(palette)) palette = colorRampPalette(c("blue", "white", "red"))(500)
    palette[length(palette)] = diagonalcolor
    # The palette is long on purpose, so that we only
    # get black if corr > 0.996.
    lmat = lwid = lhei = 1
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    par(mar = c(5,10,10,5))
    z <- cormat[, nvar:1] # had to flip for image().
    # transformation to emphasize big abs(cor):
    z = z*(abs(z)<=0.5) + (1.5*z^2+0.125)*sign(z)*(abs(z)>0.5)
    if(plot.abs.cor == TRUE) z <- abs(z)
    mywidths = c(0,scales)
    if(cell.length == "sqrt") mywidths = c(0,sqrt(scales))
    if(cell.length == "equal") mywidths = c(0,rep(1,nvar))
    mywidths = pmax(mywidths, 0.02*max(mywidths))
    # print(mywidths)
    hbreaks = cumsum(mywidths)
    hbreaks = 15*hbreaks/max(hbreaks)
    vbreaks = max(hbreaks) - rev(hbreaks)
    #
    # Compute required width of margins:
    char_width_inch <- par("cin")[1]
    # = width of a single character in inches
    label_width_inch = max(strwidth(varnames,
                                    units = "inches"))
    margin_lines = label_width_inch/char_width_inch
    par(mar = c(0.3, margin_lines, margin_lines, 0.3))
    #        bottom       left         top      right
    #
    image(x=hbreaks, y=vbreaks, z=z,
          xlim = range(hbreaks), ylim = range(vbreaks),
          zlim = c(-1.625,1.625), axes = FALSE,
          xlab = "", ylab = "", col = palette)
    hmidpoints = (hbreaks[-(nvar+1)] + hbreaks[-1])/2
    vmidpoints = (vbreaks[-(nvar+1)] + vbreaks[-1])/2
    axis(side = 3, at = hmidpoints, labels = varnames,
         tick = FALSE, line = -0.7, cex.axis = 1,
         gap.axis = -1, las = 2)
    axis(side = 2, at = vmidpoints, labels = rev(varnames),
         tick = FALSE, line = -0.7, cex.axis = 1,
         gap.axis = -1, las = 2)
    abline(v = hbreaks)
    abline(h = vbreaks)
  }
  invis = list(cormat = cormat,
               predterms = predterms,
               predsummary = out$predsummary)
  invisible(invis)
}
