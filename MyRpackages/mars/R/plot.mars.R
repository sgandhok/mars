#' @title plot.mars method
#' @description plots the fitted basis that depends on one explanatory variable or two explanatory variables.
#'
#' @param x mars object x
#' @param ... other arguments like line size and color
#' @family methods
#'
#' @return returns a plot of fitted values
#' @export
#'
#' @examples m<- mars(wage~age+education,data=ISLR::Wage)
#' plot(m,col="red")
#' class(m)<-"lm";plot(m) #standard diagnostic plots
plot.mars<-function(x,...){
  data<-eval(x$call$data)
  tt<-terms(x$formula,data=data)
  tt<-delete.response(tt)
  mf<-model.frame(tt,data)
  mt<-attr(mf,"terms")
  X<-model.matrix(mt,mf)[,-1]
  Bf<-x$Bfuncs
  singleB<-which(sapply(Bf, function(x) NROW(x)==1))
  doubleB<-which(sapply(Bf, function(x) NROW(x)==2))
  nn<-ceiling(sqrt(length(singleB)+length(doubleB)))
  opar<-graphics::par(mfrow=c(nn,nn),mar=c(2,2,2,2))
  on.exit(graphics::par(opar))
  for(i in singleB){
    vv<-Bf[[i]][1,"v"];
    varname<-x$x_names[[vv]]
    xx<-seq(from=min(X[,vv]),to=max(X[,vv]),length=100)
    bb<-h(xx,Bf[[i]][1,"s"],Bf[[i]][1,"t"])
    plot(xx,bb,type="l",xlab=varname,main=varname,...)
  }

}
