#' @title Relevance by ghost variables in linear regression 
#' @name relev.ghost.var
#' 
#' @description  \code{\link{relev.ghost.var}} Computing the case-variable relevance case
#' \code{matrix A} by using the Ghost Variables strategy.
#' 
#' @return  
#' \code{relev.ghost.var} returns a list with following elements:
#' \itemize{
#' \item \code{A} { n x p case-variable relevance matrix}
#' \item \code{V} { relevance matrix}
#' \item \code{GhostX} { j-th ghost variable that it was replaced}
#' \item \code{relev.ghost} { diagonal of the relevance matrix}
#' \item \code{eig.V} { eigen values of the relevance matrix}
#' \item \code{y.hat.test} { test set}
#' }
#' 
#' @param model ajusted model with training set. The variable relevance is computted for this model.
#' @param newdata another model in which the variable of interest is substituted by its ghost variable, defined as the prediction of this variable by using the rest of explanatory variables.
#' @param func.model.ghost.var model we want to explain to build ghost variables
#' 
#' @seealso \code{\link{relev.rand.perm}} and \code{\link{relevance2cluster}}
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
#' plot.relev.ghost.var(relev.ghost.out, n1=n1, resid.var=sum.lm.tr$sigma^2, sum.lm.tr=sum.lm.tr)
#' 
#' @rdname relev.matrix.ghost
#' @export 

relev.ghost.var <- function(model,newdata=model$call$data,
  func.model.ghost.var=gam, ...){ #gam (Fits a generalized additive model)
  #func.model <- eval(parse(text=class(model)[1])) # What kind of model has been fitted 
  #data.tr <- model$call$data[model$call$subset,] # data used for training the model
  #n <- dim(data.tr)[1]
  ifelse( isS4(model),  #getting the varaible names in the model
          term.labels <- attr(model@terms,"term.labels"),
          term.labels <- attr(model$terms,"term.labels")
  )
  if (identical(mod.gh,"gam")|identical(mod.gh,gam)){
    s.term.labels <- paste0("s(",term.labels,")")
  }else{
    s.term.labels <- term.labels
  }
  p <- length(term.labels)
  
  n2 <- dim(newdata)[1] # newdata is the test sample (mostra de prova)
  # Predicting in the test sample
  y.hat.ts <- as.numeric( predict(model,newdata = newdata) )
  
  A <- matrix(0,nrow=n2, ncol=p)
  colnames(A) <- term.labels
  GhostX <- A
  yX.ts <- newdata
  for (j in (1:p)){
    yX.ts.aux <- yX.ts 
    form.j <- as.formula(paste0(term.labels[j],"~",paste(s.term.labels[-j],collapse="+")))
    xj.hat <- func.model.ghost.var(form.j, data=yX.ts[,term.labels],...)$fitted.values
    yX.ts.aux[,term.labels[j]] <- xj.hat
    y.hat.ts.j <- as.numeric( predict(model,newdata = yX.ts.aux) )
    A[,j] <- y.hat.ts - y.hat.ts.j
    GhostX[,j] <- xj.hat
  }
  V=(1/n2)*t(A)%*%A  #Relevance matrix
  relev.ghost <- diag(V)
  eig.V <- eigen(V)
  return(list(A=A, V=V, GhostX=GhostX, 
    relev.ghost=relev.ghost, 
    eig.V=eig.V,
    y.hat.test=y.hat.ts)
  )
}

