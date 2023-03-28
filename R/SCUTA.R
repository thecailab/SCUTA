
#' fit linear mixed model for Single Cell RNA-seq Time Series data with multilevel structure
#' @importFrom lmerControl  lme4
#' @importFrom nextElem     iterators
#' @importFrom pbkrtest     get_SigmaG
#' @import     foreach
#' @importFrom variancePartition     colinearityScore
SCUTA<-function (exprObj, formula, data, L, ddf = c("Satterthwaite",
                                             "Kenward-Roger"), useWeights = TRUE, weightsMatrix = NULL,
          control = lme4::lmerControl(calc.derivs = FALSE, check.rankX = "stop.deficient"),
          suppressWarnings = FALSE, quiet = FALSE, BPPARAM = bpparam(),
          computeResiduals = FALSE, REML = TRUE, ...)
{
  exprObjInit = exprObj
  exprObjMat = as.matrix(exprObj)
  formula = stats::as.formula(formula)
  ddf = match.arg(ddf)
  colinearityCutoff = 0.999
  if (ncol(exprObj) != nrow(data)) {
    stop("the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)")
  }
  if (is(exprObj, "EList")) {
    if (useWeights) {
      weightsMatrix = exprObj$weights
    }
    .checkNA(exprObj$E)
  }else {
    .checkNA(exprObj)
  }
  if (!(ddf %in% c("Kenward-Roger", "Satterthwaite"))) {
    stop("Specify ddf correctly")
  }
  if (ddf == "Kenward-Roger" & !REML) {
    stop("Kenward-Roger must be used with REML")
  }
  if (useWeights && is.null(weightsMatrix)) {
    useWeights = FALSE
  }
  if (useWeights && !identical(dim(exprObj), dim(weightsMatrix))) {
    stop("exprObj and weightsMatrix must be the same dimensions")
  }
  if (!useWeights) {
    weightsMatrix = NULL
  }
  if (!identical(colnames(exprObj), rownames(data))) {
    warning("Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently.")
  }
  if (.isDisconnected()) {
    stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
  }

  univariateContrasts = FALSE
  if (missing(L)) {
    L = .getAllUniContrasts(exprObj, formula, data)
    univariateContrasts = TRUE
  } else {
    if (is(L, "numeric")) {
      L = as.matrix(L, ncol = 1)
    }
    else if (is(L, "data.frame")) {
      L = as.matrix(L)
    }
    if (is.null(colnames(L))) {
      colnames(L) = paste0("L", seq_len(ncol(L)))
    }
    tst = apply(L, 2, function(x) {
      length(x[x != 0]) == 1
    })
    if (any(tst)) {
      warning("Contrasts with only a single non-zero term are already evaluated by default.")
    }
    Luni = .getAllUniContrasts(exprObj, formula, data)
    L = cbind(L, Luni)
  }
  if (ncol(L) == 0) {
    stop("Must include fixed effect in the model for hypothesis testing")
  }
  if (length(unique(colnames(L))) != ncol(L)) {
    stop(paste("Contrast names must be unique: ", paste(colnames(L),
                                                        collapse = ", ")))
  }
  if (!.isMixedModelFormula(formula)) {
    if (!quiet) {
      message("Fixed effect model, using limma directly...")
      message("User can apply eBayes() afterwards...")
    }
    design = model.matrix(formula, data)
    ret = lmFit(exprObj, design, weights = weightsMatrix)
    if (computeResiduals) {
      ret$residuals = residuals(ret, exprObj)
    }
    if (!univariateContrasts) {
      ret = contrasts.fit(ret, L)
    }
  }else {
    form = paste("responsePlaceholder$E", paste(as.character(formula),
                                                collapse = ""))
    responsePlaceholder = nextElem(exprIter(exprObjMat, weightsMatrix,
                                            useWeights))
    timeStart = proc.time()
    fitInit <- lmerTest::lmer(eval(parse(text = form)), data = data,
                              ..., REML = REML, control = control)
    # fitInit <- lmerTest::lmer(eval(parse(text = form)), data = data,
    #                           REML = REML, control = control)
    sigGStruct = get_SigmaG(fitInit)$G
    if (!identical(rownames(L), names(fixef(fitInit)))) {
      stop("Names of entries in L must match fixed effects")
    }
    mod = .eval_lmm(fitInit, L, ddf)
    timediff = proc.time() - timeStart
    checkModelStatus(fitInit, showWarnings = !suppressWarnings,
                     dream = TRUE, colinearityCutoff = colinearityCutoff)
    a = names(fixef(fitInit))
    b = rownames(L)
    if (!identical(a, b)) {
      stop("Terms in contrast matrix L do not match model:\n  Model: ",
           paste(a, collapse = ","), "\n  L: ",
           paste(b, collapse = ","), "\nNote thhat order must be the same")
    }
    data2 = data.frame(data, expr = responsePlaceholder$E,
                       check.names = FALSE)
    form = paste("expr", paste(as.character(formula),
                               collapse = ""))
    # pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",
                           # total = nrow(exprObj), width = 60, clear = FALSE)
    # pids = .get_pids()
    timeStart = proc.time()
    .eval_models = function(responsePlaceholder, data2, form,
                            REML, theta, control, na.action = stats::na.exclude,
                            ...) {
      data2$expr = responsePlaceholder$E
      suppressWarnings({
        fit <- lmerTest::lmer(eval(parse(text = form)),
                              data = data2, REML = REML, ..., weights = responsePlaceholder$weights,
                              control = control, na.action = na.action)
      })
      mod = .eval_lmm(fit, L, ddf)
      res = NULL
      if (computeResiduals) {
        res = residuals(fit)
      }
      ret = list(coefficients = mod$beta, design = fit@pp$X,
                 df.residual = mod$df, Amean = mean(fit@frame[,
                                                              1]), method = "lmer", sigma = mod$sigma,
                 stdev.unscaled = mod$SE/mod$sigma, pValue = mod$pValue,
                 residuals = res)
      varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit,
                                                               "stddev")^2)
      varComp[["resid"]] = attr(lme4::VarCorr(fit),
                                "sc")^2
      if (univariateContrasts) {
        V = mod$vcov
      }
      else {
        V = crossprod(chol(mod$vcov) %*% L)
      }
      list(ret = new("MArrayLM", ret), varComp = varComp,
           edf = sum(hatvalues(fit)), vcov = V)
    }

    # $$$ important improvement:previous one only allow Elist. While now $$$
    # $$$ we improve the function to accept matrix as an object          $$$
    # mute the improvement now #
    if (is(exprObj, "EList")){
      .eval_master = function(obj, data2, form, REML, theta,
                              control, na.action = stats::na.exclude, ...) {
        lapply(seq_len(nrow(obj$E)), function(j) {
          .eval_models(list(E = obj$E[j, ], weights = obj$weights[j,
          ]), data2, form, REML, theta, control, na.action,
          ...)
        })
      }
    }
    # else{
    #   .eval_master = function(obj, data2, form, REML, theta,
    #                           control, na.action = stats::na.exclude, ...) {
    #     lapply(seq_len(nrow(obj)), function(j) {
    #       .eval_models(list(E = obj[j, ]),
    #                         data2, form, REML, theta, control, na.action,
    #       ...)
    #     })
    #   }
    # }
    # $$$ important improvement:previous one only allow Elist. While now $$$
    # $$$ we improve the function to accept matrix as an object          $$$

    resList<-.eval_master(obj=exprObj, data2=data2, form=form, REML=REML, theta=fitInit@theta,
                          control=control, na.action = stats::na.exclude, ...)
    # resList<-.eval_master(obj=exprObj, data2=data2, form=form, REML=REML, theta=fitInit@theta,
    #                       control=control, na.action = stats::na.exclude)
    # it = iterBatch(exprObjMat, weightsMatrix, useWeights,
    #                n_chunks = 100)
    # if (!quiet)
    #   message(paste0("Dividing work into ", attr(it,
    #                                              "n_chunks"), " chunks..."))
    # resList <- bpiterate(it, .eval_master, data2 = data2,
    #                      form = form, REML = REML, theta = fitInit@theta,
    #                      control = control, ..., REDUCE = c, reduce.in.order = TRUE,
    #                      BPPARAM = BPPARAM)
    names(resList) = seq_len(length(resList))
    if (!quiet)
      message("\nTotal:", paste(format((proc.time() -
                                          timeStart)[3], digits = 0), "s"))
    x = 1
    coefficients = foreach(x = resList, .combine = cbind) %do%{
        x$ret$coefficients
      }
    df.residual = foreach(x = resList, .combine = cbind) %do%{
        x$ret$df.residual
      }
    pValue = foreach(x = resList, .combine = cbind) %do%{
        x$ret$pValue
      }
    stdev.unscaled = foreach(x = resList, .combine = cbind) %do% {
        x$ret$stdev.unscaled
      }
    if (computeResiduals) {
      residuals = foreach(x = resList, .combine = cbind) %do% {
          x$ret$residuals
        }
    }
    coefficients = t(coefficients)
    df.residual = t(df.residual)
    pValue = t(pValue)
    stdev.unscaled = t(stdev.unscaled)
    if (computeResiduals) {
      residuals = t(residuals)
      rownames(residuals) = rownames(exprObj)
    }
    colnames(coefficients) = colnames(L)
    rownames(coefficients) = rownames(exprObj)
    design = resList[[1]]$ret$design
    colnames(df.residual) = colnames(L)
    rownames(df.residual) = rownames(exprObj)
    colnames(pValue) = colnames(L)
    rownames(pValue) = rownames(exprObj)
    Amean = sapply(resList, function(x) x$ret$Amean)
    names(Amean) = rownames(exprObj)
    method = "lmer"
    sigma = sapply(resList, function(x) x$ret$sigma)
    names(sigma) = rownames(exprObj)
    varComp = lapply(resList, function(x) as.data.frame(x$varComp))
    varComp = do.call("rbind", varComp)
    # rownames(varComp) = rownames(coefficients)
    edf = sapply(resList, function(x) x$edf)
    colnames(stdev.unscaled) = colnames(L)
    rownames(stdev.unscaled) = rownames(exprObj)
    ret = list(coefficients = coefficients, design = design,
               df.residual = df.residual, Amean = Amean, method = method,
               sigma = sigma, contrasts = L, stdev.unscaled = stdev.unscaled)
    if ("genes" %in% names(exprObjInit)) {
      ret$genes = exprObjInit$genes
    }
    ret = new("MArrayLM", ret)
    V = chol2inv(qr(ret$design)$qr)
    rownames(V) = colnames(ret$design)
    colnames(V) = colnames(ret$design)
    if (!univariateContrasts) {
      V = crossprod(chol(V) %*% L)
    }
    ret$cov.coefficients = V
    ret$randon.var = varComp
    ret$cov.coefficients.list = lapply(resList, function(x) as.matrix(x$vcov))
    if (computeResiduals) {
      ret$residuals = residuals
    }
    ret = as(ret, "MArrayLM2")
    attr(ret, "varComp") = varComp
    attr(ret, "sigGStruct") = sigGStruct
    attr(ret, "edf") = edf
    ret = .standard_transform(ret)
  }
  ret
}



# below are functions required in modDream function
exprIter = function( exprObj, weights, useWeights = TRUE, scale=TRUE, iterCount = "icount"){

  n_features = nrow(exprObj)

  if( iterCount == 'icount2'){
    xit <- icount2( n_features )
  }else{
    xit <- icount( n_features )
  }

  nextEl <- function() {
    j <- nextElem(xit)

    if( is.null(j) || j > n_features){
      res = NULL
    }else{
      if( useWeights && !is.null(weights) ){

        w = weights[j,]

        # scale weights to have mean of 1, otherwise it affects the residual variance too much
        # scale should be false when signa(fit) needs to be evaluted
        if(scale){
          w = w / mean(w)
        }
      }else{
        w = NULL
      }

      res = list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
    }
    res
  }
  it <- list(nextElem = nextEl)
  class(it) <- c("abstractiter", "iter")
  it
}

#' @import foreach
#' @importFrom lme4 fixef 
#' @importFrom stats sigma
.eval_lmm = function( fit, L, ddf ){

  j = 1
  # evaluate each contrast
  # cons = lmerTest::contest(fit, L, ddf=ddf)
  cons = foreach( j = 1:ncol(L), .combine=rbind) %do%
  {
    lmerTest::contest(fit, L[,j], ddf=ddf)
  }

  df = as.numeric(cons[,'DenDF'])

  if(ddf == "Kenward-Roger"){
    # KR
    V = pbkrtest::vcovAdj.lmerMod(fit, 0)

    # if matrix is not PSD
    if( min(diag(as.matrix(V))) < 0){
      warning("The adjusted Kenward-Roger covariance matrix is not positive definite.\nUsing Satterthwaite approximation instead")

      # Satterthwaite
      V = vcov(fit)
    }
    # df = pbkrtest::get_Lb_ddf(fit, L)
  }else{
    # Satterthwaite
    V = vcov(fit)
    # df = as.numeric(contest(fit, L, ddf="Sat")['DenDF'])
  }

  # sigma = attr(lme4::VarCorr(fit), "sc")

  # get contrasts
  # beta = as.matrix(sum(L * fixef(fit)), ncol=1)
  # colnames(beta) = "logFC"

  beta = foreach( j = 1:ncol(L), .combine=rbind) %do%
  {
    as.matrix(sum(L[,j] * fixef(fit)), ncol=1)
  }
  colnames(beta) = "logFC"
  rownames(beta) = colnames(L)

  # SE = as.matrix(sqrt(sum(L * (V %*% L))), ncol=1)
  # colnames(SE) = "logFC"
  SE = foreach( j = 1:ncol(L), .combine=rbind) %do%
  {
    as.matrix(sqrt(sum(L[,j] * (V %*% L[,j]))), ncol=1)
  }
  colnames(SE) = "logFC"
  rownames(SE) = colnames(L)

  # pValue = 2*pt(as.numeric(abs(beta / SE)), df, lower.tail=FALSE)
  pValue = as.numeric(cons[,'Pr(>F)'])

  list(	cons 	= cons,
        df		= df,
        sigma	= sigma(fit),
        beta	= beta,
        SE		= SE,
        pValue	= pValue,
        vcov 	= V )
}




setGeneric("checkModelStatus", signature="fit",
           function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
             standardGeneric("checkModelStatus")
)

setMethod("checkModelStatus", "lm",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
          {
            # if no intercept is specified, give warning
            if( showWarnings && length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
              warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
            }

            # if any coefficient is NA
            if( showWarnings && any(is.na(coef(fit))) ){
              stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
            }

            # check colinearity
            score = colinearityScore(fit)
            if( score > colinearityCutoff ){
              stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
            }
          }
)

setMethod("checkModelStatus", "lmerMod",
          function( fit, showWarnings=TRUE, dream=FALSE, colinearityCutoff=.999 )
          {
            # if no intercept is specified, give warning
            if( !dream && showWarnings && length(which(colnames(fit@pp$X) == "(Intercept)")) == 0 ){
              warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
            }

            # if any coefficient is NA
            if( ( showWarnings | dream) && any(is.na(coef(fit))) ){
              stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
            }

            # check colinearity
            ###################
            score = colinearityScore(fit)

            if( score > colinearityCutoff ){
              stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
            }

            # check that factors are random and continuous variables are fixed
            ###################################################################

            # remove backticks with gsub manually
            # solve issue that backticks are conserved is some but not all parts of lmer()

            # Simplified testing of random versus fixed effects
            # allows (A|B) only where A is continuous

            # variables fit by regression
            testVar = attr(attr(fit@frame, "terms"), "term.labels")
            testVar = gsub("`", "", testVar)

            # get type for each variable
            # keep only tested variables
            varType = attr(attr(fit@frame, "terms"), "dataClasses")[-1]
            varType = varType[testVar]

            # random effects
            randVar = names(fit@flist)

            # fixed effects
            # starting with all variables, remove random variables
            fixedVar = setdiff(testVar, randVar)

            for( i in 1:length(varType) ){

              # if factor is not random
              if( (showWarnings && ! dream) && varType[i] %in% c("factor", "character") && (! names(varType)[i] %in% randVar) ){
                stop(paste("Categorical variables modeled as fixed effect:", paste(names(varType)[i], collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))
              }

              # If numeric/double is not fixed
              if( (showWarnings && ! dream) && varType[i] %in% c("numeric", "double") && (!names(varType)[i] %in% fixedVar) ){
                stop(paste("Continuous variable cannot be modeled as a random effect:", names(varType)[i]))
              }
            }

            # show convergance message
            if( showWarnings && !is.null(fit@optinfo$conv$lme4$messages) && (fit@optinfo$conv$lme4$messages != "boundary (singular) fit: see ?isSingular")){
              stop(fit@optinfo$conv$lme4$messages)
            }
          }
)


iterBatch <- function(exprObj, weights, useWeights = TRUE, scale=TRUE, n_chunks = nrow(exprObj) / 500, min_chunk_size = 20, BPPARAM = NULL ) {
  # Adjust number of chunks upward to the next multiple of number of
  # workers in BPPARAM, if this can be determined. If any errors are
  # encountered, just continue without adjusting.
  tryCatch(
    if (is(BPPARAM, "BiocParallelParam")) {
      n_workers <- bpworkers(BPPARAM)
      if (!is.null(n_workers) && is.numeric(n_workers) && n_workers >= 1) {
        chunks_per_worker <- ceiling(n_chunks / n_workers)
        n_chunks <- chunks_per_worker * n_workers
      }
    },
    error = function(...) NULL
  )

  # Don't split into chunks smaller than min_chunk_size
  max_allowed_chunks <- floor(nrow(exprObj) / min_chunk_size)
  n_chunks = min(n_chunks, max_allowed_chunks)
  # Make sure we have at least 1 chunk (since we can get 0 if
  # min_chunk_size > nrow)
  n_chunks <- max(n_chunks, 1)

  # specify chunks
  idx <- parallel::splitIndices(nrow(exprObj), min(nrow(exprObj), n_chunks))
  i <- 0L

  f = function() {
    if (i == length(idx)){
      return(NULL)
    }
    i <<- i + 1L
    E = exprObj[ idx[[i]],, drop = FALSE ]

    if( useWeights && !is.null(weights) ){
      # scale weights to have mean of 1, otherwise it affects the residual variance too much
      if(scale){
        # w = weights[j,] /  mean(weights[j,])
        w = weights[idx[[i]],,drop=FALSE]
        # for each row, devide by row mean
        w = w / rowMeans(w)
      }else{
        # w = weights[j,]
        w = weights[idx[[i]],,drop=FALSE]
      }
    }else{
      w = matrix(1, nrow(E), ncol(E))
    }

    list(E = E, weights = w )
  }

  # get number of chunks
  attr( f, "n_chunks") = length(idx)
  f
}


.checkNA = function(exprObj){
  # check if values are NA
  countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
  if( countNA > 0 ){
    stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
  }

  # check if all genes have variance
  rv = apply( exprObj, 1, var)
  if( any( rv == 0) ){
    idx = which(rv == 0)
    stop(paste("Response variable", idx[1], 'has a variance of 0'))
  }
}


.isDisconnected = function(){
  i = NULL
  possibleError <- tryCatch( suppressWarnings(foreach(i = seq_len(2)) %dopar% {i}), error = function(e) e)
  return( isTRUE(inherits(possibleError, "error") && identical(possibleError$message, "invalid connection")) )
}

.getAllUniContrasts = function( exprObj, formula, data){

  Linit = .getContrastInit( exprObj, formula, data)

  Lall = lapply( seq_len(length(Linit)), function(i){
    Linit[i] = 1
    Linit
  })
  names(Lall) = names(Linit)
  Lall = do.call("rbind", Lall)

  # remove intercept contrasts
  # Lall[,-1,drop=FALSE]
  Lall
}


getContrast = function( exprObj, formula, data, coefficient){

  if( length(coefficient) > 2){
    stop("Length of coefficient array limited to 2")
  }

  L = .getContrastInit( exprObj, formula, data)

  # assign coefficient coding
  if( any(!coefficient %in% names(L)) ){
    stop("coefficient is not in the formula.  Valid coef are:\n", paste(names(L), collapse=', '))
  }
  L[coefficient[1]] = 1

  if( length(coefficient) == 2){
    L[coefficient[2]] = -1
  }

  L
}

# .getContrastInit( exprObj, formula, data)

.getContrastInit = function( exprObj, formula, data){

  exprObj = as.matrix( exprObj )
  formula = as.formula( formula )

  # only retain columns used in the formula
  data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]

  REML=TRUE
  useWeights=TRUE
  weightsMatrix=NULL
  showWarnings=FALSE
  dream=TRUE
  fxn=identity
  colinearityCutoff=.999
  control = lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

  # check dimensions of reponse and covariates
  if( ncol(exprObj) != nrow(data) ){
    stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e. rows)" )
  }

  # if weightsMatrix is not specified, set useWeights to FALSE
  if( useWeights && is.null(weightsMatrix) ){
    # warning("useWeights was ignored: no weightsMatrix was specified")
    useWeights = FALSE
  }

  # if useWeights, and (weights and expression are the same size)
  if( useWeights && !identical( dim(exprObj), dim(weightsMatrix)) ){
    stop( "exprObj and weightsMatrix must be the same dimensions" )
  }

  # If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
  if( ! identical(colnames(exprObj), rownames(data)) ){
    warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
  }
  form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

  # run lmer() to see if the model has random effects
  # if less run lmer() in the loop
  # else run lm()
  responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights, scale=FALSE))
  possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,control=control ), error = function(e) e)

  mesg <- "No random effects terms specified in formula"
  method = ''
  if( isTRUE(inherits(possibleError, "error") && identical(possibleError$message, mesg)) ){

    design = model.matrix( formula, data)

    L = rep(0, ncol(design))
    names(L) = colnames(design)

    # detect error when variable in formula does not exist
  }else if( isTRUE(inherits(possibleError, "error") && length( grep("object '.*' not found", possibleError$message)) > 0) ){
    stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
  }
  else{

    if( isTRUE(inherits(possibleError, "error") && grep('the fixed-effects model matrix is column rank deficient', possibleError$message) == 1) ){
      stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
    }

    fit = lmer( eval(parse(text=form)), data=data,control=control, REML=TRUE )

    L = rep(0, length(fixef(fit)))
    names(L) = names(fixef(fit))
  }
  L
}


.isMixedModelFormula = function(formula ){

  ! is.null( findbars( as.formula( formula ) ) )
}

.standard_transform = function(fit, sigma = fit$sigma){

  # If fit$df.prior is not defined, set df.prior to zero
  if( ! is.null(fit$df.prior) ){
    fit$df.total = fit$df.residual + fit$df.prior
  }else{
    fit$df.total = fit$df.residual
  }

  # t-test
  out = fit
  out$t <- fit$coefficients / fit$stdev.unscaled / sigma
  out$p.value <- 2*pt(-abs(out$t), df=fit$df.total )

  # F-test
  if(!is.null(out$design) && is.fullrank(out$design)) {

    # only evaluate F-stat on real coefficients, not contrasts
    realcoef = colnames(out)[colnames(out) %in% colnames(out$design)]
    realcoef = realcoef[realcoef!="(Intercept)"]

    if( is.null(realcoef) || (length(realcoef) == 0) ){

      # this happends when only the intercept term is included
      warning("No testable fixed effects were included in the model.\n  Running topTable() will fail.")
    }else{
      df = rowMeans(out[,realcoef]$df.total)

      F.stat <- classifyTestsF(out[,realcoef], df=df, fstat.only=TRUE)
      out$F <- as.vector(F.stat)
      df1 <- attr(F.stat,"df1")
      df2 <- attr(F.stat,"df2")
      if(df2[1] > 1e6){ # Work around bug in R 2.1
        out$F.p.value <- pchisq(df1*out$F,df1,lower.tail=FALSE)
      }else{
        out$F.p.value <- pf(out$F,df1,df2,lower.tail=FALSE)
      }
    }
  }

  # if fit$df.prior does not exist, then remove the df.total term
  if( is.null(fit$df.prior) ){
    out$df.total = NULL
  }

  out
}



