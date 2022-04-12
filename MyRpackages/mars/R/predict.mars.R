
#' Predict Method for mars object
#'
#' Predicted Values
#' @param object a mars object
#' @param newdata optional specification of new data on which prediction could be derived
#' @param ... other arguments
#' @return a vector of prediction values
#' @export
#'
#'@examples
#'m<- mars(wage~age+education,data=ISLR::Wage)
#'predict(m)
#'
predict.mars <- function(object,newdata,...) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

make_B <- function(X,Bfuncs) {
  B <- matrix(1,nrow=nrow(X),ncol=length(Bfuncs))

   for(i in 2:length(Bfuncs)) {
     for (j in 1:nrow(Bfuncs[[i]])){
       B <- B*h(X[,Bfuncs[[i]][j,"v"]], Bfuncs[[i]][j,"s"], Bfuncs[[i]][j,"t"])
       }

     }
 B
   }

