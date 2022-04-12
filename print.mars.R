#' @title  Print method for mars object
#'
#' @param object a mars object
#' @param ... other additional arguments to print.mars
#' @family methods
#' @return  returns the call to mars and  estimated basis function coefficients.
#' @export
#'
#' @examples
#' m<- mars(wage~age+education,data=ISLR::Wage)
#' print(m)
print.mars<-function(object,...){
  cat("\n")
  cat("Call:","\n")
  print(object$call)
  cat("\n")
  cat("Coefficients:","\n")
  print(object$coefficients)
}
