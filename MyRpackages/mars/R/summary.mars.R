#' Summary method for mars object
#'
#' @param object a mars object
#' @param ... other optional arguments to be passed to summary.mars
#'
#' @return the call to mars, summary of residuals of the fitted values, estimated basis function coefficients and residual standard error
#' @export
#' @import plotrix

#' @examples
#' m<- mars(wage~age+education,data=ISLR::Wage)
#' summary(m)
summary.mars<-function(object,...){
  cat("\n")

  cat("Call:","\n")

  print(object$call)

  cat("\n")
  cat("Residuals:","\n")
  print(summary(object$residuals))

  cat("\n")
  cat("Coefficients:","\n")
  print(object$coefficients)

  cat("\n")
  cat("Residual standard error:  "  )
  print(std.error(object$residuals))
}
