
#######################
# Auxiliary functions #
#######################


# Hidden from the user:
wnq = function(string,qwrite=T){ # auxiliary function
  # writes a line without quotes
  if(qwrite) write(noquote(string),file="",ncolumns=100)
}


# Hidden from the user:
pnq = function(string,qwrite=T){ # auxiliary function
  # prints without quotes
  if(qwrite) print(noquote(string))
}


# Hidden from the user:
checkLabels = function(y, n, training=T, levels=NULL){
  if(training == T){ # for training data
    if(!is.null(levels)) {
      stop("For training data, the levels argument is not used.") }
    whichy = "y"
  } else { # for new data
    if(is.null(levels)) {
      stop("For new data, the argument levels is required") }
    whichy = "ynew"
  }
  if(!(is.factor(y))) {
    stop(paste0("\n The response ",whichy," is not a factor.")) }
  if(is.null(attr(y,"levels"))) {
    stop(paste0("\n The response factor ",whichy," has no levels.",
                "\n Please use as.factor() or factor() first.")) }
  levelsattr = attr(y,"levels")
  if(sum(is.na(levelsattr)) > 0){
    stop(paste0("\n The levels attribute of ",whichy,
                " has at least one NA.")) }
  dupls = duplicated(levelsattr)
  if(sum(dupls) > 0) stop(paste0(
    "\n The label ",levelsattr[dupls == T]," occurs more than ",
    " once in the levels attribute of ",whichy))
  if(length(y) != n) {
    stop(paste0("\n The response y should have length ",n,
                ", the number of cases.")) }
  yv = as.character(y) # this yv is not a factor any more.
  # When the factor y is malformed, the function stops here.
  # This happens when not all entries of y belong to levelsattr.
  indsv = which(!is.na(yv)) # INDiceS of y that can be Visualized
  if(training == T) { # for training data:
    if(length(indsv)==0) {
      stop(paste0("The response factor y only contains NA's,"
                  ,"\n so training is not possible.")) }
    yv = yv[indsv] # thus has at least 1 entry
    uniqy = sort(unique(yv)); uniqy
    for(g in seq_len(length(uniqy))){ # g=1
      if(!(uniqy[g] %in% levelsattr)) {
        stop(paste0("\n The response label ",uniqy[g]," does not",
                    " appear in the levels attribute of y.")) }
    }
    # From here on we know that all the yv are in levelsattr.
    if(length(uniqy) < 2){
      stop(paste0("\n Only a single level occurs ( ", uniqy[1],
                  " ) so training is not possible.")) }
    for(g in seq_len(length(levelsattr))){ # g=1
      if(!(levelsattr[g] %in% uniqy)){
        wnq(paste0("\nThe level ",levelsattr[g]," does not occur",
                   " among the training cases. A model trained",
                   "\non these data will not be able to",
                   " assign any objects to this class.\n")) }
    }
    # Create array "levels" with the labels that actually occur:
    levels = levelsattr[which(levelsattr %in% uniqy)]; levels
    #
  } else { # for new data:
    if(length(indsv) > 0){
      yv = yv[indsv] # thus has at least 1 entry
      uniqy = sort(unique(yv)); uniqy
      for(g in seq_len(length(uniqy))){ # g=1
        if(!(uniqy[g] %in% levelsattr)) {
          stop(paste0("\n The response label ",uniqy[g]," does not",
                      " appear in the levels attribute of ynew.")) }
      }
      # From here on we know that all the yv are in levelsattr.
      badlevels = uniqy[which(!(uniqy %in% levels))]; badlevels
      if(length(badlevels) > 0) { wnq(paste0(
        "\n The level ", badlevels," occurs in ynew but",
        " was not in the training data.",
        "\n Such cases will be treated as if their ynew is NA",
        " so their response","\n can be predicted, but",
        " these cases will not be in the class map.\n"))
        indsv = indsv[which(!(yv %in% badlevels))]
        # take those out too
        # the resulting indsv may again be empty
      }
    }
  }
  xlevelz = levels # for use in the following function:
  lab2int = function(labs){
    # labs is a vector of labels, that may contain NA's
    ints = rep(NA,length(labs))
    indv = which(!is.na(labs))
    labv = labs[indv] # uses the final indsv above
    for(g in seq_len(length(xlevelz))){ # g=1
      clinds = indv[which(labv == xlevelz[g])] # in all of labs
      ints[clinds] = g
    }
    ints
  }
  return(list(lab2int = lab2int, levels=levels, indsv=indsv))
}


# Hidden from the user:
meancov = function(X) {
  # Wrapper around colMeans and classical covariance matrix
  return(list(m = colMeans(X), S = cov(X)))
}


# Hidden from the user:
DetMCD = function(X) {
  # Wrapper around robustbase::covMCD
  detmcd.out = robustbase::covMcd(X,nsamp="deterministic",alpha=0.5)
  return(list(m = detmcd.out$center, S = detmcd.out$cov))
}


# Hidden from the user:
numericalCond = function(S) {
  # Computes the inverse condition number of a matrix, which
  # is the ratio: smallest eigenvalue/largest eigenvalue .
  # Unlike the CN, this avoids division by (near) zero.
  # I call this the "numerical condition" (NC).
  #
  S = as.matrix(S)
  if(length(which(is.na(S))) > 0) stop(" S contains NA's.")
  d = nrow(S)
  if(ncol(S) != d) stop(" S is not square.")
  maxS = max(abs(as.vector(S)))
  if(!isSymmetric(S, tol=1e-10*maxS)) stop(" S is not symmetric.")
  eigvals = eigen(S, only.values=T)$values
  eigvals[d]/eigvals[1]
}


# Visible to the user (is called in example script):
makeKernel = function(X1, X2 = NULL, svfit){
  # Computes kernel value or kernel matrix, where the
  # kernel type is extracted from svfit.
  #
  # Arguments:
  # X1    : first data (coordinate) matrix or vector.
  # X2    : if not NULL, second data matrix or vector.
  # svfit : an object returned by e1071::svm(), or a list
  #         with the members $scaled, $kernel, $degree,
  #         $gamma, and $coef0 .
  #         Here $kernel may be in c(0,1,2,3) or in
  #         c("linear","polynomial","radial","sigmoid").
  #
  # A heuristic choice for gamma in radial and sigmoid is
  #    0.5/median(as.vector(as.matrix(dist(Xmat)))^2)
  #
  # Value:
  # kmat  : the kernel matrix or kernel value.
  #
  if(is.vector(X1)) { X1 = t(as.matrix(X1)) } else {
    X1 = as.matrix(X1) }
  d = ncol(X1)
  if(!is.null(X2)){
    if(is.vector(X2)) { X2 = t(as.matrix(X2)) } else {
      X2 = as.matrix(X2) }
    if(ncol(X2) != d) stop(" X1 and X2 have different dimensions")
  }
  if(is.null(svfit$kernel)) {
    stop(" svfit$kernel is missing") } else {
      k = svfit$kernel }
  if(!(k==0 | k=="linear" | k==1 | k=="polynomial" |
       k==2 | k=="radial" | k==3 | k=="sigmoid")) {
    stop(paste(" Unknown kernel type: ",k,sep="")) }
  if(is.null(svfit$scaled)) {
    stop(" svfit$scaled is missing") } else {
      scaled = svfit$scaled }
  if(is.null(svfit$degree)) {
    stop(" svfit$degree is missing") } else {
      degree = svfit$degree }
  if(is.null(svfit$gamma)) {
    stop(" svfit$gamma is missing") } else {
      gamma = svfit$gamma }
  if(is.null(svfit$coef0)) {
    stop(" svfit$coef0 is missing") } else {
      coef0 = svfit$coef0 }
  if(is.null(svfit$cost)) {
    stop(" svfit$cost is missing") } else {
      cost = svfit$cost }
  #
  # scaling if this was done in call to svm()
  #
  if(length(scaled) == 1) { scaling = rep(scaled,d) }
  if(length(scaled) != d) {
    stop(paste(" svfit$scaled should have length ",d,sep="")) }
  # Now we scale the corresponding columns:
  sind = which(scaled == T)
  if(length(sind) > 0){
    if(is.null(svfit$x.scale)) {
      stop(" svfit$x.scale is missing")
    } else {
      ctrs = svfit$x.scale[[1]]
      scls = svfit$x.scale[[2]]
      X1[,sind] = scale(X1[,sind],center=ctrs,scale=scls)
      if(!is.null(X2)){
        X2[,sind] = scale(X2[,sind],center=ctrs,scale=scls) }
    }
  }
  #
  # kernel definitions
  #
  if(k==0 | k=="linear"){ kern = function(xi,yj) {
    sum(xi*yj) } }
  if(k==1 | k=="polynomial"){ kern = function(xi,yj) {
    (gamma*sum(xi*yj) + coef0)^degree } }
  if(k==2 | k=="radial"){ kern = function(xi,yj) {
    exp(-gamma*sum((xi-yj)*(xi-yj))) } }
  if(k==3 | k=="sigmoid"){ kern = function(xi,yj) {
    tanh(gamma*sum(xi*yj) + coef0) } }
  n1 = nrow(X1)
  if(is.null(X2)) { # compute n1 by n1 matrix kmat
    kmat = matrix(rep(0,n1*n1),nrow=n1)
    dim(kmat)
    # Only fill triangular matrix and then put in mirror
    # image to obtain whole symmetric matrix:
    for(i in seq_len(n1)){
      for(j in seq_len(i)) {
        kmat[i,j] = kern(X1[i,],X1[j,]) }
    }
    kmat = kmat + t(kmat)
    diag(kmat) = 0.5*diag(kmat)
  }
  if(!is.null(X2)) { # compute n1 by n2 matrix kmat
    n2 = nrow(X2)
    kmat = matrix(rep(0,n1*n2),nrow=n1)
    for(i in seq_len(n1)){
      for(j in seq_len(n2)) { kmat[i,j] = kern(X1[i,],X2[j,]) }
    }
  }
  if(nrow(kmat)==1 & ncol(kmat)==1) kmat = kmat[1,1]
  return(kmat)
}


# Visible to the user:
makeFV = function(kmat,transfmat=NULL,precS=1e-12){
  #
  # Constructs feature vectors from a kernel matrix.
  #
  # Arguments:
  # kmat      : a kernel matrix. If transfmat is NULL, we are dealing
  #             with training data and then kmat must be a square
  #             kernel matrix (of size n by n when there are n cases).
  #             Such a PSD matrix kmat can e.g. be produced by
  #             makeKernel() or by kernlab::kernelMatrix().
  #             If on the other hand transfmat is not NULL, we are
  #             dealing with a test set. Denote the number of cases
  #             in the test set by m >= 1.
  #             Each row of kmat of the test set then must contain
  #             the kernel values of a new object with all objects
  #             in the training set. Therefore the kernel matrix
  #             kmat must have dimensions m by n. The matrix kmat
  #             can e.g. be produced by makeKernel(X,Z). It can
  #             also be obtained by running kernlab::kernelMatrix()
  #             on the union of the training set and the test set,
  #             yielding an (n+m) by (n+m) matrix from which one
  #             then takes the [(n+1):m , 1:n] submatrix.
  # transfmat : transformation matrix. If not NULL, it is the value
  #             $transfmat of makeFV on training data. It has to
  #             be a square matrix, with as many rows as there
  #             were training data.
  # precS     : if not NULL, eigenvalues (possibly negative) below
  #             precS will be set equal to precS.
  #
  # Value:
  # Xf        : When makeKV() is applied to the training set, Xf has
  #             coordinates of n points (vectors), the plain inner
  #             products of which equal the kernel matrix of the
  #             training set. That is, kmat = Xf%*%t(Xf).
  #             When makeKV() is applied to a test set, Xf_test are
  #             coordinates of the feature vectors of the test set
  #             in the same space as those of the training set, and
  #             then kmat ~= (Xf_test) %*% t(Xf_training).
  # transfmat : square matrix for transforming kmat to Xf.
  #
  kmat = as.matrix(kmat)
  if(length(which(is.na(kmat))) > 0) {
    stop(" The kernel matrix kmat contains NA's.") }
  n = nrow(kmat)
  maxK = max(abs(as.vector(kmat)))
  if(is.null(transfmat)){
    if(ncol(kmat) != n) {
      stop(" The training kernel matrix kmat is not square.") }
    if(maxK < 1e-3) {
      stop(" The kernel matrix is too close to zero.") }
    if(max(abs(as.vector(kmat - t(kmat)))) > 1e-10*maxK) {
      stop(" The training kernel matrix kmat is not symmetric.") }
    defaultPrecS = 1e-16 # was 0
    if(is.null(precS)){ precS = defaultPrecS }
    if(precS < defaultPrecS){ precS = defaultPrecS }
    # There is no previous transfmat, so we construct one.
    eig = eigen(kmat)
    eigvals = eig$values
    eig = eig$vectors # overwrite to save space
    V = eig; rm(eig) # R's rename
    keep = length(which(eigvals > precS))
    if(keep < ncol(kmat)) {
      wnq(paste0(" The kernel matrix has ",ncol(kmat)-keep,
                 " eigenvalues below precS = ",precS,"."))
      wnq(" Will take this into account.") }
    eigvals = pmax(eigvals,precS)
    scals = sqrt(eigvals) # all strictly positive
    V = (V%*%diag(1/scals)%*%t(V))
    transfmat = V; rm(V) # R's rename
    rm(eigvals) # we'll output scals
    Xf = kmat%*%transfmat
    # This Xf is already a valid representation of the feature
    # space of the training set. This result is only determined
    # up to an orthogonal transformation.
    # As a particular orthogonal transformation we express the
    # feature vectors in the coordinate system of their principal
    # components, to increase the chances that plots of the first
    # few coordinates show interesting patterns.
    #
    eig = eigen(t(Xf)%*%Xf)
    eig = as.matrix(eig$vectors) # overwrite to save space
    loadings = eig; rm(eig) # R's rename
    colnames(loadings) = NULL # to avoid long colnames later
    flipcolumn = function(x){
      if(x[which.max(abs(x))] < 0) {-x} else {x} }
    loadings = apply(loadings,2L,FUN=flipcolumn) # is unique
    loadings = transfmat%*%loadings
    transfmat = loadings; rm(loadings) # to save space
    attr(transfmat,"keep") = keep
    Xf = (kmat%*%transfmat)[,seq_len(keep)]
  } else {
    # There is a previous transfmat available that we can use.
    dim = nrow(transfmat) # dimension of training kernel matrix
    if(ncol(kmat)!=dim) {
      stop(paste(" The matrix kmat has ",ncol(kmat)," columns",
                 " but the argument transfmat was made for ",dim,
                 " columns.",sep="")) }
    if(is.null(attr(transfmat,"keep"))) { keep = dim } else {
      keep = attr(transfmat,"keep") }
    Xf = (kmat%*%transfmat)[,seq_len(keep)]
    scals = NULL
  }
  attr(Xf,"origin") = "makeFV"
  return(list(Xf = Xf, transfmat = transfmat, scals = scals))
}


# Hidden from the user:
compFarness = function(type="affine", testdata=F, yint, nlab,
                       X=NULL, fig=NULL, d=NULL, PCAfits=NULL,
                       keepPCA=F, figparams = NULL){
  #
  # Computes farness measures of different types.
  #
  # Arguments:
  #   type  : when "affine" this will standardize the Mahalanobis
  #           distances relative to each class.
  #           When "pca" it will compute the score distances and
  #           orthogonal distances to the PCA of each class.
  #   yint  : the class of each object, coded as 1, 2, ...
  #   X     : if not NULL, a data matrix with the x_i as rows.
  #   fig   : if not NULL, initial farness of each object i
  #           to each class g.
  #   d     : if not NULL, the dimension of the x_i.

  # Auxiliary functions:

  far2probability = function(farness, trainInds = NULL) {
    #
    # converts farness values into probabilities
    # args:
    #   farness: a vector of farness values
    #   trainInds: if not NULL, these will be used in the estimation
    # returns:
    #   probs: the farness values converted to probabilities
    #   tfunc: a function which transforms the farness according to
    #   the estimated transformation (including all standardizations)
    #

    YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
      indlow <- which(y < 0)
      indhigh <- which(y >= 0)
      if (lambda != 0) {
        y[indhigh] = ((1 + y[indhigh])^(lambda) - 1)/lambda
      }
      else {
        y[indhigh] = log(1 + y[indhigh])
      }
      if (lambda != 2) {
        y[indlow] = -((1 - y[indlow])^(2 - lambda) - 1)/(2 -
                                                           lambda)
      }
      else {
        y[indlow] = -log(1 - y[indlow])
      }
      if (stdToo) {
        if (length(y) > 1) {
          locScale <- cellWise::estLocScale(matrix(y, ncol = 1),
                                            type = "hubhub")
          zt <- (y - locScale$loc)/locScale$scale
        }
        else {
          zt <- y
        }
      }
      else {
        zt <- NULL
      }
      return(list(yt = y, zt = zt))
    }

    origfarness = farness
    if (!is.null(trainInds)) {
      farness = farness[trainInds]
    }
    indnz  = which(farness > 1e-10)
    indnz
    farnz  = farness[indnz]
    farnz
    farloc = median(farnz, na.rm = T)
    farloc
    farsca = mad(farnz, na.rm = T)
    farsca
    if (farsca < 1e-10) { farsca = sd(farnz, na.rm = T) }
    sfar   = scale(farnz, center = farloc, scale = farsca)
    sfar
    # Fit a Yeo-Johnson transformation:
    YJ.out  = cellWise::transfo(X = sfar, robust = T,
                                prestandardize = F,
                                checkPars = list(silent = T))
    xt      = YJ.out$Xt
    tfarloc = median(xt, na.rm = T)
    tfarsca = mad(xt, na.rm = T)
    lambda  = YJ.out$lambdahats
    #
    # Now transform origfarness (to be used in fig):
    origIndnz = which(origfarness > 1e-10)
    origIndnz
    origFarnz = origfarness[origIndnz]
    origFarnz
    zt        = scale(YJ(scale(origFarnz, farloc, farsca),
                                    lambda = lambda, stdToo = FALSE)$yt,
                      tfarloc, tfarsca)
    probs     = rep(0,length(origfarness))
    probs[origIndnz] = pnorm(zt)
    probs
    tfunc = function(qs) {
      #
      # In Titanic data, 52% of casualties had duplicated features.
      # Therefore, we take out zero farness values in YJ computation:
      #
      qsIndnz = which(qs > 1e-10)
      qsnz = qs[qsIndnz]
      zn = scale(YJ(scale(qsnz, farloc, farsca),
                               lambda = lambda, stdToo = FALSE)$yt,
                 tfarloc, tfarsca)
      pn = rep(0,length(qs))
      pn[qsIndnz] = pnorm(zn)
      return(pn)
    }
    return(list(probs = probs, far2prob = tfunc))
  }


  transformYJ <- function(X, lambda) {
    # Yeo-Johnson transformation with given parameter lambda.
    # The argument must be a vector or a number, without NA's.
    #
    if(!is.vector(X)){
      stop("transformYJ: the data should be a vector or a number") }
    if(sum(is.na(X)) > 0) stop("transformYJ: the data contains NA's")
    Xt = rep(NA, length(X))
    indlow  <- which(X < 0)
    indhigh <- which(X >= 0)
    if (lambda != 0) {
      Xt[indhigh] = ((1 + X[indhigh])^(lambda) - 1) / lambda
    } else {
      Xt[indhigh] = log(1 + X[indhigh])
    }
    if (lambda != 2) {
      Xt[indlow] = -((1 - X[indlow])^(2 - lambda) - 1) / (2 - lambda)
    } else {
      Xt[indlow] = -log(1 - X[indlow])
    }
    return(Xt)
  }


  invYJ = function(y, lambda) {
    # Inverse of the Yeo-Johnson transformation.
    if(y >= 0) {
      if (lambda != 0) {
        inv = (1 + lambda*y)^(1/lambda) - 1
      } else { inv = exp(y) - 1 }
    }
    if(y < 0){
      if (lambda != 2) {
        inv = -((1+(lambda-2)*y)^(1/(2-lambda))) + 1
      } else { inv = 1 - exp(-y) }
    }
    return(inv)
  }


  mednz = function(x, prec=1e-12){
    # Median of values > 0. Can be useful in denominators.
    #
    if(!(prec > 0)) stop("mednz: prec = ",prec," should be > 0")
    xx = x[which(!is.na(x))]
    xx = xx[which(xx > prec)]
    if(length(xx) == 0) { mednzx = prec } else { mednzx = median(xx) }
    # By construction, always mednz >= prec .
    mednzx
  }


  # Here the main function compFarness starts:
  #
  if(!(type %in% c("affine","knn","pca"))){
    stop(' type must be either "affine" or "knn" or "pca".') }
  n          = length(yint)
  indsv      = which(!is.na(yint))
  yintv      = yint[indsv] # yint of cases that can be visualized
  ayint      = sort(unique(yint[indsv])) # the values it takes
  farness    = rep(NA,n) # will remain NA for cases outside indsv
  classSizes = rep(0,nlab)
  for(g in seq_len(nlab)){ classSizes[g] = sum(yintv == g) }
  medOD = medSD = SDs = ODs = NULL
  #
  if(type == "affine"){
    if(is.null(fig)){ stop("argument fig is missing") }
    #
    if(testdata == F){ # this is for training data
      for(g in ayint){
        clinds = which(yint == g) # indices in all of 1,...,n
        # Indices for which yint is NA do not get in clinds.
        farness[clinds] = fig[clinds,g]
      }
      omedian = mednz(farness, prec=1e-08) # overall median
      #
      # Renormalize the classes:
      farscales = rep(NA, nlab)
      for(g in ayint){
        clinds          = which(yint == g) # in 1,...,n
        gfarness        = farness[clinds] # all non-NA
        farscales[g]    = mednz(gfarness, prec=1e-08)/omedian
        farness[clinds] = gfarness/farscales[g]
      }
      #
      # Fit a distribution to the homogenized farness measures:
      mfar    = median(farness, na.rm=T)
      sfar    = max(mad(farness, na.rm=T),1e-08) # no division by 0
      farness = as.vector(scale(farness, center=mfar, scale=sfar))
      tout    = transfo(farness, type="YJ", robust=T,
                        prestandardize=F, checkPars=list(silent=T))
      lambda    = tout$lambdahat[1]
      muhat     = tout$muhat
      sigmahat  = tout$sigmahat
      figparams = list(farscales = farscales,
                       mfar = mfar,
                       sfar = sfar,
                       lambda = lambda,
                       muhat = muhat,
                       sigmahat = sigmahat)
    } else { # This is for new data
      # for new data, figparams should be given:
      if(is.null(figparams)) stop("Argument figparams is missing")
      farscales = figparams$farscales
      if(is.null(farscales)) stop("figparams$farscales is missing")
      mfar     = figparams$mfar
      sfar     = figparams$sfar
      lambda   = figparams$lambda
      muhat    = figparams$muhat
      sigmahat = figparams$sigmahat
    }
    fig = scale(fig,center=F,scale=farscales) # for training and newdata
    #
    # Now convert these distances to probabilities:
    for(g in seq_len(nlab)){ # cycle over columns of fig[,]
      fig[,g] = as.vector(scale(fig[,g],center=mfar,scale=sfar))
      # transform by YJ using lambda:
      fig[,g] = transformYJ(fig[,g],lambda)
      fig[,g] = as.vector(scale(fig[,g],center=muhat,scale=sigmahat))
      fig[,g] = pnorm(fig[,g])
    }
    for(g in ayint){ # for training and newdata
      clinds = which(yint == g)
      farness[clinds] = fig[clinds,g]
      # so it remains NA outside of indsv
    }
  } # end of type == "affine"
  #
  if(type == "knn"){
    if(is.null(fig)){ stop("argument fig is missing") }
    #
    if(testdata == F){ # this is for training data
      # compute initial farness and overall median
      # farness = rep(NA,n) # already defined above
      for (g in ayint) {
        clinds = which(yint == g) # indices in 1,...,n
        farness[clinds] = fig[clinds,g]
      }
      omedian = mednz(farness, prec=1e-08) # overall median
      #
      # Renormalize the classes:
      farscales = rep(NA,nlab)
      for(g in ayint){
        clinds       = which(yint == g) # in all of 1,...,n
        gfarness     = farness[clinds] # all non-NA
        farscales[g] = mednz(gfarness, prec=1e-08)/omedian
        # farness[clinds] = gfarness/farscales[g]
      }
      # wnq("farscales:")
      # pnq(farscales)
      fig = scale(fig, center = F, scale = farscales)
      farness = rep(NA,n)
      far2probFuncs = list()
      # list of estimated functions to transform raw farness values
      probs = matrix(NA, nrow = nrow(fig), ncol = ncol(fig))
      for (g in seq_len(nlab)) {
        clinds             = which(yint == g) # class indices
        ftemp              = fig[, g]
        far2prob.out       = far2probability(ftemp, trainInds = clinds)
        fig[, g]           = far2prob.out$probs
        farness[clinds]    = fig[clinds, g]
        far2probFuncs[[g]] = far2prob.out$far2prob
      }
      figparams = list(farscales = farscales,
                       far2probFuncs = far2probFuncs)
    } else { # this is for new data
      # for new data, figparams should be given:
      if(is.null(figparams)) stop("Argument figparams is missing")
      farscales = figparams$farscales
      if(is.null(farscales)) stop("figparams$farscales is missing")
      fig = scale(fig, center=F, scale=farscales)
      #
      # transform the fig values using the far2prob functions:
      far2probFuncs = figparams$far2probFuncs
      for(g in seq_len(nlab)){
        far2probTemp = far2probFuncs[[g]]
        fig[,g]      = far2probTemp(fig[, g])
      }
      # now compute farness for all new data in indsv:
      for(g in ayint){ # g=3 # loop over classes in new data
        clinds = which(yint == g)
        farness[clinds] = fig[clinds,g]
        # so it remains NA outside of indsv
      }
    } # end of new data
  } # end of type == "knn"
  #
  if(type == "pca"){
    if(is.null(X)){ stop("argument X is missing") }
    if(is.vector(X)) { X = matrix(X, ncol = 1) }
    X = as.matrix(X)
    d = ncol(X)
    if(nrow(X) != n){ stop(" X and yint have incompatible sizes") }
    fig = SDs = ODs = matrix(rep(NA,n*nlab),ncol=nlab)
    #
    if(testdata == F){ # this is for training data
      if(n < 2) stop(" X should have more than one row")
      PCAfits = NULL
      if(keepPCA) PCAfits = list()
      medSD = medOD = rep(NA,nlab)
      medscores = madscores = list()
      for(j in seq_len(nlab)){ # Loop over the classes # j=1
        clinds = which(yint == j) # class indices
        # Cases with missing yint do not enter any clinds.
        Xj = X[clinds,,drop=F] # class j might contain only 1 case
        nXj = length(clinds) # there can be small classes
        nXju = sum(duplicated(Xj)==F)
        cntr = colMeans(Xj) # full dimension
        outpca = prcomp(sweep(Xj,2,cntr),scale=F)
        if(nXju == 1){ outpca$rotation = outpca$rotation*0 }
        sdev = outpca$sdev # square roots of eigenvalues
        tokeep = length(which(sdev > 1e-8))
        ncomps = max(tokeep,1)
        qmax = NULL
        if(!is.null(qmax)) ncomps = min(ncomps, qmax)
        # wnq(paste(" Keeping ",ncomps,
        #           " principal components in class ",j,sep=""))
        sdev = sdev[seq_len(ncomps)]
        loadings = outpca$rotation[,seq_len(ncomps),drop=F]
        # columns are eigenvectors
        # store loadings, cntr, and sdev per class:
        if(keepPCA) PCAfits[[j]] = list(lj=loadings,cj=cntr,sj=sdev)
        #
        dim(Xj); length(cntr); dim(loadings)
        pX = sweep(X,2,cntr)%*%loadings # projects all cases of X
        pX = matrix(pX,nrow = n) # so it stays a matrix when ncomps=1
        dim(pX)
        #
        ## For orthogonal distances:
        if(tokeep == d){ # then no orthogonal distances exist
          ODs[,j] = rep(0,n) # for cases in all classes
        } else {
          Xdiffs = sweep(X,2,cntr) - pX%*%t(loadings)
          # differences in full-dimensional space, for all cases
          ODs[,j] = sqrt(rowSums(Xdiffs^2)) # Euclidean
          ODs[clinds,j] = rep(0,length(clinds))
          # cases have no orthogonal distance to their own class
        }
        medOD[j] = mednz(ODs[setdiff(indsv,clinds),j],prec=1e-8)
        ODs[,j] = ODs[,j]/medOD[j]
        # for all cases, but =zero in class j
        #
        ## For score distances:
        scores = pX[clinds,,drop=F]
        # projected, ONLY for cases in class j
        dim(scores)
        # distances from robustly scaled scores:
        medscores[[j]] = apply(scores,2,median)
        madscores[[j]] = apply(scores,2,mad)
        madscores[[j]] = pmax(madscores[[j]],1e-8)
        # clips at same precision as above, avoids dividing by 0
        scscores = scale(pX,center=medscores[[j]],scale=madscores[[j]])
        # we will use this for all cases, not only in class j
        SDs[,j]  = sqrt(rowSums(scscores^2))
        medSD[j] = mednz(SDs[clinds,j], prec=1e-8)
        # only over cases IN class j
        SDs[,j] = SDs[,j]/medSD[j] # again for all cases
        # Special situation where sdev[1] ~ 0:
        if(tokeep == 0) SDs[,j] = rep(1,length(SDs[,j]))
        ## Now combine SD with OD:
        fig[,j] = sqrt(SDs[,j]^2 + ODs[,j]^2)
      } # ends loop over classes j
      #
      for(g in seq_len(nlab)){ # g=1
        clinds = which(yint == g) # indices in 1,...,n
        farness[clinds] = fig[clinds,g] # only for available yint
      }
      omedian = mednz(farness[indsv], prec=1e-08) # overall median
      # omedian # 1 # This is because the farness is SDs^2 at this
      # stage, as ODs=0 in each class, and the SDs have already
      # been standardized to median 1 in each class. If no
      # farness is below 1e-08, mednz is the median, and if all
      # classes have odd cardinality we get omedian=1.
      #
      # Renormalize the classes to bring farness on the same scale:
      farscales = rep(NA,nlab)
      for(g in ayint){
        clinds       = which(yint == g) # in 1,...,n
        gfarness     = farness[clinds] # all non-NA
        farscales[g] = mednz(gfarness, prec=1e-08)/omedian
        # These will also be 1 when no gfarness is below 1e-08,
        # otherwise not.
        farness[clinds] = gfarness/farscales[g]
      }
      # Fit a distribution to the homogenized farness measures:
      mfar      = median(farness[indsv], na.rm=T)
      sfar      = mad(farness[indsv], na.rm=T)
      sfar      = max(sfar,1e-08) # avoids division by 0
      farness   = as.vector(scale(farness, center=mfar, scale=sfar))
      tout      = transfo(farness[indsv], type="YJ", robust=T,
                          prestandardize=F,
                          checkPars=list(silent=T))
      lambda    = tout$lambdahat[1]
      muhat     = tout$muhat
      sigmahat  = tout$sigmahat
      figparams = list(medOD = medOD,
                       medscores = medscores,
                       madscores = madscores,
                       medSD = medSD,
                       farscales = farscales,
                       mfar = mfar,
                       sfar = sfar,
                       lambda = lambda,
                       muhat = muhat,
                       sigmahat = sigmahat)
      # end of computing initial fig on training data
    } else { # test data, so we use the PCA fits of the training data
      #
      if(is.null(PCAfits)){
        stop(paste0("\nFor test data, the farness computation ",
                    "requires the trained object","\nto contain ",
                    "the value $PCAfits. Rerun the training ",
                    "with keepPCA = T"))
      }
      if(is.null(figparams)) stop("Argument figparams is missing")
      # if(is.null(farscales)) stop("Argument farscales is missing")
      medOD     = figparams$medOD
      madscores = figparams$madscores
      medscores = figparams$medscores
      medSD     = figparams$medSD
      farscales = figparams$farscales
      #
      # Computations for new data apply the existing figparams:
      #
      for(j in seq_len(nlab)){ # Loop over the classes # j=1
        loadings  = PCAfits[[j]]$lj
        cntr      = PCAfits[[j]]$cj
        sdev      = PCAfits[[j]]$sj
        tokeep    = length(which(sdev > 1e-8))
        # tokeep can still be zero when ncomps is 1
        ncomps = ncol(loadings)
        # wnq(paste0(" Using ",ncomps,
        #            " principal components in class ",j))
        clinds = which(yint == j) # class indices
        Xj = X[clinds,,drop=F] # class j can contain only 1 case
        dim(Xj); length(cntr); dim(loadings)
        pX = sweep(X,2,cntr)%*%loadings # projects all cases of X
        scores = pX[clinds,] # projected, ONLY for cases in class j
        dim(pX); dim(scores)
        Xdiffs = sweep(X,2,cntr) - pX%*%t(loadings) # for ALL cases
        # differences in full-dimensional space
        ODs[,j] = sqrt(rowSums(Xdiffs^2))
        ODs[,j] = ODs[,j]/medOD[j] # for all cases, but =zero in class j
        scscores = scale(pX,center=medscores[[j]],scale=madscores[[j]])
        SDs[,j] = sqrt(rowSums(scscores^2))
        SDs[,j] = SDs[,j]/medSD[j] # again for all cases
        if(tokeep == 0) SDs[,j] = rep(1,length(SDs[,j]))
        fig[,j] = sqrt(SDs[,j]^2 + ODs[,j]^2)
      } # ends loop over classes j
    } # end of computing initial fig on test data
    #
    # The following is for both training data and test data:
    #
    fig = scale(fig,center=F,scale=farscales)
    #
    mfar      = figparams$mfar
    sfar      = figparams$sfar
    lambda    = figparams$lambda
    muhat     = figparams$muhat
    sigmahat  = figparams$sigmahat
    # Now convert these distances to probabilities:
    for(g in seq_len(nlab)){ # cycle over columns of fig[,]
      fig[,g] = as.vector(scale(fig[,g],center=mfar,scale=sfar))
      # transform by YJ using lambda:
      fig[,g] = transformYJ(fig[,g],lambda)
      fig[,g] = as.vector(scale(fig[,g],center=muhat,scale=sigmahat))
      fig[,g] = pnorm(fig[,g])
    }
    for(g in ayint){ # g=1
      clinds = which(yint == g)
      farness[clinds] = fig[clinds,g]
      # so it remains NA outside of indsv
    }
  } # end of type == "pca"
  #
  ofarness = apply(fig,1,min) # ofarness even exists for cases
  # with missing yint, for which farness cannot be defined.
  if(!keepPCA) PCAfits = NULL
  return(list(fig = fig,
              farness = farness,
              ofarness = ofarness,
              figparams = figparams,
              PCAfits = PCAfits,
              medSD = medSD,
              medOD = medOD,
              SDs = SDs,
              ODs = ODs))
}
