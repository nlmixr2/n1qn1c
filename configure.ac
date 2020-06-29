AC_INIT([n1qn1], 6.0.1.8) n1qn1 optimization
AC_PROG_F77

N1QN1_FFLAGS=""
if test "$ac_cv_f77_compiler_gnu" == yes; then
   N1QN1_FFLAGS="-O0"
fi

AC_SUBST(N1QN1_FFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
   
echo --------[[begin src/Makevars]]--------
cat src/Makevars
echo --------[[end src/Makevars]]--------
