##' n1qn1 optimization
##'
##' This is an R port of the n1qn1 optimization procedure in scilab.
##'
##' @param call_eval Objective function
##' @param call_grad Gradient Function
##' @param vars Initial starting point for line search
##' @param environment Environment where call_eval/call_grad are
##'     evaluated.
##' @param ... Ignored additional parameters.
##' @param epsilon Precision of estimate
##' @param max_iterations Number of iterations
##' @param nsim Number of function evaluations
##' @param imp Verbosity of messages.
##' @param invisible boolean to control if the output of the minimizer
##'     is suppressed.
##' @param zm Prior Hessian (in compressed format; This format is
##'     output in \code{c.hess}).
##' @param restart Is this an estimation restart?
##' @param assign Assign hessian to c.hess in environment environment?
##'     (Default FALSE)
##' @param print.functions Boolean to control if the function value
##'     and parameter estimates are echoed every time a function is
##'     called.
##'
##' @return The return value is a list with the following elements:
##'     \itemize{
##'        \item \code{value} The value at the minimized function.
##'        \item \code{par} The parameter value that minimized the function.
##'        \item \code{H} The estimated Hessian at the final parameter estimate.
##'        \item \code{c.hess} Compressed Hessian for saving curvature.
##'        \item \code{n.fn} Number of function evaluations
##'        \item \code{n.gr} Number of gradient evaluations
##'      }
##' @author C. Lemarechal, Stephen L. Campbell, Jean-Philippe
##'     Chancelier, Ramine Nikoukhah, Wenping Wang & Matthew L. Fidler
##' @importFrom Rcpp evalCpp
##' @useDynLib n1qn1, .registration=TRUE
##' @examples
##'
##' ## Rosenbrock's banana function
##' n=3; p=100
##'
##' fr = function(x)
##' {
##'     f=1.0
##'     for(i in 2:n) {
##'         f=f+p*(x[i]-x[i-1]**2)**2+(1.0-x[i])**2
##'     }
##'     f
##' }
##'
##' grr = function(x)
##' {
##'     g = double(n)
##'     g[1]=-4.0*p*(x[2]-x[1]**2)*x[1]
##'     if(n>2) {
##'         for(i in 2:(n-1)) {
##'             g[i]=2.0*p*(x[i]-x[i-1]**2)-4.0*p*(x[i+1]-x[i]**2)*x[i]-2.0*(1.0-x[i])
##'         }
##'     }
##'     g[n]=2.0*p*(x[n]-x[n-1]**2)-2.0*(1.0-x[n])
##'     g
##' }
##'
##' x = c(1.02,1.02,1.02)
##' eps=1e-3
##' n=length(x); niter=100L; nsim=100L; imp=3L;
##' nzm=as.integer(n*(n+13L)/2L)
##' zm=double(nzm)
##'
##' (op1 <- n1qn1(fr, grr, x, imp=3))
##'
##' ## Note there are 40 function calls and 40 gradient calls in the above optimization
##'
##' ## Now assume we know something about the Hessian:
##' c.hess <- c(797.861115,
##'             -393.801473,
##'             -2.795134,
##'             991.271179,
##'             -395.382900,
##'             200.024349)
##' c.hess <- c(c.hess, rep(0, 24 - length(c.hess)))
##'
##' (op2 <- n1qn1(fr, grr, x,imp=3, zm=c.hess))
##'
##' ## Note with this knowledge, there were only 29 function/gradient calls
##'
##' (op3 <- n1qn1(fr, grr, x, imp=3, zm=op1$c.hess))
##'
##' ## The number of function evaluations is still reduced because the Hessian
##' ## is closer to what it should be than the initial guess.
##'
##' ## With certain optimization procedures this can be helpful in reducing the
##' ## Optimization time.
##'
##' @export
n1qn1 <- function(call_eval, call_grad, vars, environment=parent.frame(1), ...,
                  epsilon=.Machine$double.eps, max_iterations=100, nsim=100,
                  imp=0,
                  invisible=NULL,
                  zm=NULL, restart=FALSE,
                  assign=FALSE,
                  print.functions=FALSE){
    if (!is.null(invisible)){
        if (invisible == 1){
            imp <- 0;
            print.functions <- FALSE
        } else {
            print.functions <- TRUE
        }
    }
    if (!missing(max_iterations) && missing(nsim)){
        nsim <- max_iterations * 10;
    }
    n <- as.integer(length(vars));
    imp <- as.integer(imp);
    max_iterations <- as.integer(max_iterations)
    nsim <- as.integer(nsim)
    nzm <- as.integer(n * (n + 13) / 2);
    nsim <- as.integer(nsim);
    epsilon <- as.double(epsilon)
    if (is.null(zm)){
        mode <- 1L
        zm <- double(nzm);
    } else {
        mode <- 2L
        if (restart) model <- 3L
        if (length(zm) != nzm){
            stop(sprintf("Compressed Hessian not the right length for this problem.  It should be %d.", nzm))
        }
    }
    ret <- .Call(n1qn1_wrap, call_eval, call_grad, environment,
                 vars, epsilon, n, mode, max_iterations, nsim, imp, nzm, zm, as.integer(print.functions));
    if (assign) environment$c.hess <- ret$hess;
    return(ret)
}


##' qnbd optimization
##'
##' This is an R port of the qnbd which is a BFGS-B optimization
##' procedure in scilab. (R has L-BFGS-B).
##'
##' @param par Initial parameter estimate
##' @param fn Function
##' @param gr Gradient
##' @param lower Lower Bound for optimization
##' @param upper Upper Bound for optimization
##' @param environment
##' @param zero Tolerance for Zero
##' @param maxFn Maximum function evaluations
##' @param maxIt Maximum iterations
##' @param epsf Function eps for exiting
##' @param epsg Gradient eps for exiting
##' @param epsx Parameter eps for exiting
##' @param print.functions
##' @inheritParams n1qn1
##' @export
##' @examples
##'
##' ## Rosenbrock's banana function
##' n=3; p=100
##'
##' fr = function(x)
##' {
##'     f=1.0
##'     for(i in 2:n) {
##'         f=f+p*(x[i]-x[i-1]**2)**2+(1.0-x[i])**2
##'     }
##'     f
##' }
##'
##' grr = function(x)
##' {
##'     g = double(n)
##'     g[1]=-4.0*p*(x[2]-x[1]**2)*x[1]
##'     if(n>2) {
##'         for(i in 2:(n-1)) {
##'             g[i]=2.0*p*(x[i]-x[i-1]**2)-4.0*p*(x[i+1]-x[i]**2)*x[i]-2.0*(1.0-x[i])
##'         }
##'     }
##'     g[n]=2.0*p*(x[n]-x[n-1]**2)-2.0*(1.0-x[n])
##'     g
##' }
##'
##' x = c(1.02,1.02,1.02)
##'
##' op1 <- qnbd(x, fr, grr)
##'
##' @export
qnbd <- function(par, fn, gr, lower= -Inf, upper=Inf, environment=parent.frame(1),
                 zero=sqrt(.Machine$double.eps/7e-07), maxFn=10000L, maxIt=10000L,
                 epsf=sqrt(.Machine$double.eps), epsg=sqrt(.Machine$double.eps),
                 epsx=sqrt(.Machine$double.eps), print.functions=FALSE){
    n <- length(par);
    if (length(lower) == 1) lower <- rep(lower, n);
    if (length(upper) == 1) upper <- rep(upper, n);
    print(par);
    ret <- .Call(qnbd_wrap, fn, gr, environment, par, lower, upper, zero, as.integer(maxFn),
                           as.integer(maxIt), epsf, epsg, epsx, as.integer(n), as.integer(print.functions));
    return(ret);
}
