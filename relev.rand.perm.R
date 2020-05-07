#' @title Relevance by random permutations in the linear and additive models 
#' @name relev.rand.perm
#' 
#' @description  \code{\link{relev.rand.perm}} Computing the case-variable relevance case
#' \code{matrix A} using the random permutation strategy.
#' 
#' @return  
#' \code{relev.rand.perm} function that returns a list with following elements:
#' \itemize{
#' \item \code{A}{n x p case-variable relevance matrix}
#' \item \code{V}{relevance matrix}
#' \item \code{relev.rp}{diagonal of the relevance matrix}
#' \item \code{eig.V}{eiguen values of the relevance matrix}
#' \item \code{y.hat.test}{test set}
#' \item \code{diag.cov.X}{variances of the explanatory variables}
#' }
#' @param model ajusted model with training set. The variable relevance is computted for this model.
#' @param newdata another model in which random permutation are used in the relevance \code{matrix V}.
#' @param nperm random permutacions.
#' 
#' @seealso \code{\link{relev.ghost.var}} and \code{\link{relevance2cluster}}
#' @examples 
#' n1 <- 2000                                           # size of the training sample 
#' n2 <- 1000                                           # size of the test sample
#' sigma.1 <- 1                                         # sd for x1
#' sigma.2 <- 1                                         # sd for x2
#' sigma.3 <- 1                                         # sd for x3
#' sigma.eps <- 1                                       # residual sd for defining y
#' rho <- .95                                           # correlation between x2 and x3
#' beta1 <- 1                                           # coef. of y=x_1+...+x_{p1}
#' beta2 <- 1                                           # coef. of y=x_1+...+x_{p2}
#' beta3 <- 1                                           # coef. of y=x_1+...+x_{p2}
#' X1 <- sigma.1 * matrix(rnorm(n1+n2),ncol=1)          # Generating variables x2 and x3
#' Sigma.2 <- matrix(rho, nrow=2, ncol=2)               # Generating variables x2 and x3
#' diag(Sigma.2) <- 1
#' eig.Sigma.2 <- eigen(Sigma.2)
#' sqrt.Sigma.2 <- eig.Sigma.2$vectors %*% diag(eig.Sigma.2$values^.5) %*% t(eig.Sigma.2$vectors)
#' X23 <- matrix(rnorm((n1+n2)*2),ncol=2) %*% sqrt.Sigma.2 %*%diag(c(sigma.2,sigma.3))
#' X2<-X23[,1]
#' X3<-X23[,2]
#' y <- beta1*X1 + beta2*X2 + beta3*X3 + rnorm(n1+n2,sd=sigma.eps) # defining the response variable
#' X <- cbind(X1,X2,X3)
#' colnames(X) <- paste0("x",1:3)
#' yX <- as.data.frame(cbind(y,X))
#' colnames(yX) <- c("y",paste0("x",1:3))
#' tr.sample <- (1:n1) # Training sample:
#' lm.tr <- lm(y ~ ., data=yX, subset = tr.sample)                    # Fitting the linear model
#' (sum.lm.tr <- summary(lm.tr))
#' y.hat.ts <- as.numeric( predict(lm.tr,newdata = yX[-tr.sample,]) ) # Predicting in the test sample
#' relev.ghost.out <- relev.ghost.var(model=lm.tr, newdata = yX[-tr.sample,], func.model.ghost.var= lm)
#' relev.rand.out <- relev.rand.perm(model=lm.tr, newdata = yX[-tr.sample,], func.model.ghost.var= lm)
#' plot.relev.rand.perm(relev.rand.out, relev.ghost=relev.ghost.out$relev.ghost, sum.lm.tr=sum.lm.tr)
#' 
#' @rdname relev.rand.perm
#' @export

relev.rand.perm <- function(model,newdata=model$call$data,
  nperm=1, ...){
  #func.model <- eval(parse(text=class(model)[1])) # What kind of model has been fitted 
  #data.tr <- model$call$data[model$call$subset,] # data used for training the model
  #n <- dim(data.tr)[1]
  attr.model <- attributes(model$terms) #getting the varaible names in the model
  term.labels <- attr.model$term.labels #getting the varaible names in the model
  p <- length(term.labels)
  
  n2 <- dim(newdata)[1] # newdata is the test sample
  # Predicting in the test sample
  y.hat.ts <- as.numeric( predict(model,newdata = newdata) )
  
  A <- matrix(0,nrow=n2, ncol=p)
  colnames(A) <- term.labels
  diag.cov.X <- numeric(p)
  yX.ts <- newdata
  for (perm in 1:nperm){
    for (j in (1:p)){
      yX.ts.aux <- yX.ts
      xj.perm <- sample(yX.ts[,term.labels[j]])
      diag.cov.X[j] <- var(xj.perm)
      yX.ts.aux[,term.labels[j]] <- xj.perm
      y.hat.ts.j <- as.numeric( predict(model,newdata = yX.ts.aux) )
      A[,j] <- A[,j] + y.hat.ts - y.hat.ts.j
    }
  }
  A <- A/nperm
  V=(1/n2)*t(A)%*%A
  relev.rp <- diag(V)
  eig.V <- eigen(V)
  return(list(A=A, V=V, 
    relev.rp=relev.rp, 
    eig.V=eig.V,
    y.hat.test=y.hat.ts, 
    diag.cov.X=diag.cov.X)
  )
}
