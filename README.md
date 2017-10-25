[![Travis-CI Build Status](https://travis-ci.org/nlmixrdevelopment/n1qn1.svg?branch=master)](https://travis-ci.org/nlmixrdevelopment/n1qn1) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mattfidler/n1qn1?branch=master&svg=true)](https://ci.appveyor.com/project/mattfidler/n1qn1)
[![Coverage Status](https://img.shields.io/codecov/c/github/nlmixrdevelopment/n1qn1/master.svg)](https://codecov.io/github/nlmixrdevelopment/n1qn1?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/n1qn1)](https://cran.r-project.org/package=n1qn1)

*Ported By:* Matthew Fidler, Wenping Wang

*Algorithm Authors:* C. Lemarechal, Stephen L. Campbell, Jean-Philippe
Chancelier, Ramine Nikoukhah

R port of the Scilab n1qn1 module.  This package provides n1qn1, or
Quasi-Newton BFGS "qn" without constraints.  This takes more memory
than traditional L-BFGS.  This particual routine is useful since it
allows prespecification of a Hessian; If the Hessian is near the truth
in optimization it can speed up the optimization problem.
