# n1qn1 optimization

This is an R port of the n1qn1 optimization procedure in scilab.

## Usage

``` r
n1qn1(
  call_eval,
  call_grad,
  vars,
  environment = parent.frame(1),
  ...,
  epsilon = .Machine$double.eps,
  max_iterations = 100,
  nsim = 100,
  imp = 0,
  invisible = NULL,
  zm = NULL,
  restart = FALSE,
  assign = FALSE,
  print.functions = FALSE
)
```

## Arguments

- call_eval:

  Objective function

- call_grad:

  Gradient Function

- vars:

  Initial starting point for line search

- environment:

  Environment where call_eval/call_grad are evaluated.

- ...:

  Ignored additional parameters.

- epsilon:

  Precision of estimate

- max_iterations:

  Number of iterations

- nsim:

  Number of function evaluations

- imp:

  Verbosity of messages.

- invisible:

  boolean to control if the output of the minimizer is suppressed.

- zm:

  Prior Hessian (in compressed format; This format is output in
  `c.hess`).

- restart:

  Is this an estimation restart?

- assign:

  Assign hessian to c.hess in environment environment? (Default FALSE)

- print.functions:

  Boolean to control if the function value and parameter estimates are
  echoed every time a function is called.

## Value

The return value is a list with the following elements:

- `value` The value at the minimized function.

- `par` The parameter value that minimized the function.

- `H` The estimated Hessian at the final parameter estimate.

- `c.hess` Compressed Hessian for saving curvature.

- `n.fn` Number of function evaluations

- `n.gr` Number of gradient evaluations

## Author

C. Lemarechal, Stephen L. Campbell, Jean-Philippe Chancelier, Ramine
Nikoukhah, Wenping Wang & Matthew L. Fidler

## Examples

``` r
## Rosenbrock's banana function
n=3; p=100

fr = function(x)
{
    f=1.0
    for(i in 2:n) {
        f=f+p*(x[i]-x[i-1]**2)**2+(1.0-x[i])**2
    }
    f
}

grr = function(x)
{
    g = double(n)
    g[1]=-4.0*p*(x[2]-x[1]**2)*x[1]
    if(n>2) {
        for(i in 2:(n-1)) {
            g[i]=2.0*p*(x[i]-x[i-1]**2)-4.0*p*(x[i+1]-x[i]**2)*x[i]-2.0*(1.0-x[i])
        }
    }
    g[n]=2.0*p*(x[n]-x[n-1]**2)-2.0*(1.0-x[n])
    g
}

x = c(1.02,1.02,1.02)
eps=1e-3
n=length(x); niter=100L; nsim=100L; imp=3L;
nzm=as.integer(n*(n+13)/2)
zm=double(nzm)

(op1 <- n1qn1(fr, grr, x, imp=3))
#> $value
#> [1] 1
#> 
#> $par
#> [1] 1 1 1
#> 
#> $H
#>              [,1]      [,2]         [,3]
#> [1,]  799.9959094 -399.6141   -0.1962135
#> [2,] -399.6140752 1002.5950 -400.3193166
#> [3,]   -0.1962135 -400.3193  202.1709062
#> 
#> $c.hess
#>  [1]  799.9959094 -399.6140752   -0.1962135 1002.5949739 -400.3193166
#>  [6]  202.1709062    0.0000000    0.0000000    0.0000000    0.0000000
#> [11]    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000
#> [16]    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000
#> [21]    0.0000000    0.0000000    0.0000000    0.0000000
#> 
#> $n.fn
#> [1] 40
#> 
#> $n.gr
#> [1] 40
#> 

## Note there are 40 function calls and 40 gradient calls in the above optimization

## Now assume we know something about the Hessian:
c.hess <- c(797.861115,
            -393.801473,
            -2.795134,
            991.271179,
            -395.382900,
            200.024349)
c.hess <- c(c.hess, rep(0, 24 - length(c.hess)))

(op2 <- n1qn1(fr, grr, x,imp=3, zm=c.hess))
#> $value
#> [1] 1
#> 
#> $par
#> [1] 1 1 1
#> 
#> $H
#>               [,1]      [,2]          [,3]
#> [1,]  800.03080771 -399.8782   -0.05266924
#> [2,] -399.87816045 1001.8405 -399.89053754
#> [3,]   -0.05266924 -399.8905  201.93261767
#> 
#> $c.hess
#>  [1]  800.03080771 -399.87816045   -0.05266924 1001.84045503 -399.89053754
#>  [6]  201.93261767    0.00000000    0.00000000    0.00000000    0.00000000
#> [11]    0.00000000    0.00000000    0.00000000    0.00000000    0.00000000
#> [16]    0.00000000    0.00000000    0.00000000    0.00000000    0.00000000
#> [21]    0.00000000    0.00000000    0.00000000    0.00000000
#> 
#> $n.fn
#> [1] 29
#> 
#> $n.gr
#> [1] 29
#> 

## Note with this knowledge, there were only 29 function/gradient calls

(op3 <- n1qn1(fr, grr, x, imp=3, zm=op1$c.hess))
#> $value
#> [1] 1
#> 
#> $par
#> [1] 1 1 1
#> 
#> $H
#>              [,1]      [,2]         [,3]
#> [1,]  795.2593270 -397.1187   -0.1459179
#> [2,] -397.1186752 1000.8727 -400.2294075
#> [3,]   -0.1459179 -400.2294  202.1576894
#> 
#> $c.hess
#>  [1]  795.2593270 -397.1186752   -0.1459179 1000.8726637 -400.2294075
#>  [6]  202.1576894    0.0000000    0.0000000    0.0000000    0.0000000
#> [11]    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000
#> [16]    0.0000000    0.0000000    0.0000000    0.0000000    0.0000000
#> [21]    0.0000000    0.0000000    0.0000000    0.0000000
#> 
#> $n.fn
#> [1] 33
#> 
#> $n.gr
#> [1] 33
#> 

## The number of function evaluations is still reduced because the Hessian
## is closer to what it should be than the initial guess.

## With certain optimization procedures this can be helpful in reducing the
## Optimization time.
```
