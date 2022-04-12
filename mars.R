

#' @title  Multivariate Adaptive Regression Splines
#'@description Mars package run a step wise selection algorithm to look for clusters in a data and return a basis function. There are multiple methods like print, summary and plot to visualizing the mars object
#' @param formula an R formula to perform the regression on
#' @param data  data set containing the data to perform mars on
#' @param control object of class mars.control
#' @param ... additional argument to pass to mars()
#' @details call the mars function as follows:-

#' @return the mars() would return an object of class 'mars' that includes the basis function obtained through final regression. the mars object has plot,predict and summary methods for mars object.
#'
#'
#' @export
#' @examples
#'  mars(wage~age,data=ISLR::Wage)
#' @seealso [plot.mars] for plotting the results
#' @seealso [summary.mars] and [print.mars] for summarizing and printing mars object
#' @seealso [predict.mars] to obtain prediction values
#' @import stats
#' @import ISLR
#' @author Simarprit Singh
#' @references Jerome H Friedman. "Multivariate Adaptive Regression Splines"
#' Ann. Statist. 19 (1) 1-67,March, 1991. \url{https://doi.org/10.1214/aos/1176347963}
#' @references Brad McNeney STAT360 Advanced R for Data Science
#'
#'

mars <- function(formula,data,control=NULL,...) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  x_names <- colnames(x)
  if(is.null(control)) control <- mars.control()
  control<-validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}


# constructor for mars,control

new_mars.control <- function(control) {
  #stopifnot(is.list(data))
  structure(control,
            class="mars.control")

}

# validator for mars,control
validate_mars.control <- function(control) {

  if(control$Mmax<2) {
    warning("Input Mmax must be >= 2; setting to 2")
    Mmax <- 2}
  stopifnot(is.numeric(control$d))
  stopifnot(is.logical(control$trace))
  stopifnot(is.integer(control$Mmax))

  if(control$Mmax%%2 != 0){
    control$Mmax<-ceiling(control$Mmax)
  }
  control
}

#' @title helper for mars.control object
#' @param Mmax mmax
#' @param d d
#' @param trace trace
#'
#' @return returns an object of class mars.control
#' @export
#'
#' @examples
#' mc<-mars.control(Mmax=10,d=3,trace=FALSE)
#' mc
mars.control<- function(Mmax=2,d=3,trace=FALSE) {
  Mmax<-c(as.integer(Mmax))
  control<-list(Mmax=Mmax,d=d,trace=trace)
  d<-c(as.integer(d))
  trace<-c(as.logical(trace))
  data<-list(Mmax=Mmax,d=d,trace=trace)
  new_mars.control(control)
}


# constructor for mars,control



####################lab 5


##week 05 sol

fwd_stepwise <- function(y,x,control=mars.control()){
  #---------------------------------------------------
  # Error checking:
  #---------------------------------------------------
  # Initialize:
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- init_B(N,control$Mmax) # Exercise: write init_B()
  Bfuncs<- vector(mode="list",length=control$Mmax+1)

  #---------------------------------------------------
  # Looping for forward selection:
  for(i in 1:(control$Mmax/2)) { # contrast to indexing 2...Mmax in Friedman
    M<-2*i-1
    lof_best <- Inf
    for(m in 1:M) {# choose a basis function to split
      svars<-setdiff(1:n,Bfuncs[[m]][,"v"])
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m]) # Exercise: write split_points()
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1=B[,m]*h(x[,v],+1,t),
                             Btem2=B[,m]*h(x[,v],-1,t))

          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat,control) #  Use your LOF() from week 4
          if(lof < lof_best) {

            lof_best <- lof
            split_best<-c(m=m,v=v,t=t)

          } # end if
        } # end loop over splits
      } # end loop over variables
    } # end loop over basis functions to split

    #cat("best split",split_best,"\n")
    m<-split_best["m"];v<-split_best["v"];t<-split_best["t"]

    Bfuncs[[M+1]]<- rbind(Bfuncs[[m]], c(s=-1,v,t))

    Bfuncs[[M+2]]<- rbind(Bfuncs[[m]], c(s=1,v,t))

    B[,M+1:2] <- cbind(B[,m]*h(x[,v],-1,t),B[,m]*h(x[,v],+1,t))

  } # end loop over i
  colnames(B) <- paste0("B",(0:(ncol(B)-1)))

  return(list(y=y,B=B,Bfuncs=Bfuncs))
}


init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}
#H <- function(x) {
#  return(as.numeric(x>=0))
#}
h <- function(x,s,t) {
  return(pmax(0,s*(x-t)))
}

bwd_stepwise <- function(fwd,control) {
  Mmax <- ncol(fwd$B)-1
  Jstar <- 2:(Mmax+1)
  Kstar <- Jstar
  dat <- data.frame(y=fwd$y,fwd$B)
  lofstar <- LOF(y~.-1,dat, control)
  for(M in (Mmax+1):2) {
    b <- Inf
    L <- Kstar
    for(m in L){
      K <- setdiff(L,m)
      dat <- data.frame(y=fwd$y,fwd$B[,K])
      lof <- LOF(y~.,dat,control)
      if(lof < b) {
        b <- lof
        Kstar <- K
      }
      if(lof < lofstar) {
        lofstar <- lof
        Jstar <- K
      }
    }
  }
  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}
LOF <- function(form,data,control) {
  ff <- lm(form,data)
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  M <- length(coef(ff))-1
  Ctilde <- sum(diag(hatvalues(ff))) + control$d*M
  return(RSS * N/(N-Ctilde)^2)
}
#-------------------------------------
#Solutions are from here down


init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  #names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}
split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])


}
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testmars.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/Testfiles/fwdtest.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/Testfiles/mctest.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testLOF.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testfwd_stepwise.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testbwd_stepwise.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testfwd_stepwise.RData")
# load("C:/Users/Simar/OneDrive/Desktop/SFUStat360/Exercises/ProjectTestfiles/testthat/testmc.RData")
#
#
# set.seed(123); n <- 10
#  data <- data.frame(x1=rnorm(n),x2=rnorm(n),y=rnorm(n))
#  formula <- formula(y ~.)
#
#  mc <- mars.control()
#  all.equal(mc,mctest)
#
#  mf <- model.frame(formula,data)
#  y <- model.response(mf)
#  mt <- attr(mf, "terms")
#  x <- model.matrix(mt, mf)[,-1]
#  fwd <- fwd_stepwise(y,x,mc)
#
#  # fwd
# all.equal(fwd,fwdtest)
#
# dat <- data.frame(y=testfwd$y,testfwd$B)
#  ff <- lm(y~.,dat)
#  lof <- LOF(y~.-1,dat,testmc)
#  all.equal(lof,testLOF)
#
#  bwd <- bwd_stepwise(testfwd,testmc)
#  all.equal(bwd,testbwd)
#
#  formula <- formula(y ~.)
#  tmars<-mars(formula,marstestdata,testmc)
