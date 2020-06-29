/* zqnbd.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
#include <float.h>
#include <math.h>

#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

typedef /* Subroutine */ int (*S_fp)();
typedef /* Subroutine */ int (*U_fp)();

/*     Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/*     Copyright (C) INRIA */

/*     Copyright (C) 2012 - 2016 - Scilab Enterprises */

/* This file is hereby licensed under the terms of the GNU GPL v2.0, */
/* pursuant to article 5.3.4 of the CeCILL v.2.1. */
/* This file was originally licensed under the terms of the CeCILL v2.1, */
/* and continues to be available under such terms. */
/* For more information, see the COPYING file which you should have received */
/* along with this program. */

/* Subroutine */ int zqnbd_(int *indqn, S_fp simul, double *dh, 
	int *n, double *binf, double *bsup, double *x, 
	double *f, double *g, double *zero, int *napmax, 
	int *itmax, int *indic, int *izig, int *nfac, 
	double *epsx, double *epsf, double *epsg, double *x1, 
	double *x2, double *g1, double *dir, double *df0, 
	int *ig, int *in, int *irel, int *izag, int *iact,
	 double *epsrel, int *ieps1, int *izs, float *rzs, 
	double *dzs)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int i__, j, k;
    static double t, v, y, d1, d2;
    static int i1, n1, n3;
    static double t1, aa, dd, bi;
    static int ic, ii, ij;
    static double fn, bs, ep;
    static int mk, ip;
    static double gr;
    static int ir, np, nm1;
    static double amd, amf;
    static int ndh, nap, ifp;
    static double sig, fpn;
    static int nip;
    static double cof1, cof2, sig1, eps0, eps1;
    static int ifac;
    static double diff, difg, scal;
    extern /* Subroutine */ int rlbd_(int *, int *, S_fp, double *
	    , double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, int *,
	     int *, double *, int *, float *, double *);
    static int mode, napm;
    static double teta;
    static int iter, irit;
    extern /* Subroutine */ int proj_(int *, double *, double *, 
	    double *);
    static double tmax;
    static int nfac1;
    static double difg0, difg1;
    static int n2fac;
    static double scal1;
    static int napm1;
    static double teta1, zsig1;
    static int idfac, nnfac;
    static double epsmc;
    static int indrl, iconv;
    extern /* Subroutine */ int ajour_(int *, int *, int *, 
	    int *, double *, double *, int *);
    static double tiers, tproj, cscal1;
    extern /* Subroutine */ int calmaj_(double *, int *, double *,
	     double *, double *, int *, int *, double *, 
	    int *);
    static int indsim;

    /* Parameter adjustments */
    --dh;
    --dir;
    --g1;
    --x2;
    --x1;
    --epsx;
    --izig;
    --indic;
    --g;
    --x;
    --bsup;
    --binf;
    --izs;
    --rzs;
    --dzs;

    /* Function Body */
    difg0 = 1.;
    difg1 = 0.;

/*     eps0 sert a partitionner les variables */
    eps0 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	izig[i__] = 0;
	eps0 += epsx[i__];
    }
    eps0 = eps0 * 10.f / *n;

/*     section 1  mise en forme de dh */
/*     si indqn=1 on init dh a ident puis scal a it=2 */

    proj_(n, &binf[1], &bsup[1], &x[1]);
    ndh = *n * (*n + 1) / 2;
    if (*indqn == 1) {
	goto L10;
    }
    if (*indqn == 2) {
	goto L30;
    }
    return 0;
L10:
/*     on initialise dh a l identite puis a l iteration 2 */
/*     on met a l echelle */
    *nfac = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indic[i__] = i__;
    }
    i__1 = ndh;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dh[i__] = 0.;
    }
L30:

/*     section 2  mise a jour dh */

/*     iter nombre d iterations de descente */
    iter = 0;
    scal = 1.;
    nap = 1;
    indsim = 4;
    if (*indqn == 1) {
	(*simul)(&indsim, n, &x[1], f, &g[1], &izs[1], &rzs[1], &dzs[1]);
    }
    if (indsim <= 0) {
	*indqn = -1;
	if (indsim == 0) {
	    *indqn = 0;
	}
	return 0;
    }
    if (*indqn != 1) {
	goto L200;
    }
/*     mise a echelle dh */
/*     df0 decroissance prevue . si mod quad df0=((dh)-1g,g)/2 */
/*     et on cherche dh diag de la forme cst/(dx)**2 */
/*     d ou cst=som((y(i)*(dx))**2))/(2*df0) */
    cof1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = g[i__] * epsx[i__];
	cof1 += d__1 * d__1;
    }
    cof1 /= *df0 * 2.;
    i1 = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = i1 + *n + 2 - i__;
/* Computing 2nd power */
	d__1 = epsx[i__];
	dh[i1] = (cof1 + *zero) / (d__1 * d__1 + *zero);
    }
    iconv = 0;
L200:
    ++iter;
    if (iter <= *itmax) {
	goto L202;
    }
    *indqn = 5;
    return 0;
L202:
    if (iter == 1) {
	goto L300;
    }
    cof1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__] = x[i__] - x1[i__];
	g1[i__] = g[i__] - g1[i__];
	cof1 += x1[i__] * g1[i__];
    }
    if (cof1 <= *zero) {
	goto L250;
    }
    if (iter > 2 || *indqn != 1) {
	goto L250;
    }
    cof2 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = g1[i__];
	cof2 += d__1 * d__1;
    }
    cof2 /= cof1;
    dh[1] = cof2;
    i1 = 1;
    i__1 = *nfac;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = i1 + *n + 1 - i__;
	dh[i1] = cof2;
    }

/*     scal= (y,s)/(y,y) */
/*     scal sert de coeff a g dans le calcul de dir pour i dans ib */
    scal = 1. / cof2;
L250:

/*     mise a jour dh par methode bfgs (majour) si iter ge 2 */
/*     dh1=dh +y*yt/(y,s) - dh*s*st*dh/(s,dh*s) */
/*     exprimons ds=x1 et y=g1 dans les nouv variables soit x2 et g1 */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = indic[i__];
	x2[i1] = g1[i__];
	dir[i1] = x1[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g1[i__] = x2[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = indic[i__];
	x2[i1] = x1[i__];
    }
/*     on stocke d abord dh*s dans x2 */
/*     calcul des nfac premieres variables,en deux fois */
    if (*nfac == 0) {
	goto L2312;
    }
    if (*nfac > 1) {
	goto L2300;
    }
    dir[1] *= dh[1];
    goto L2312;
L2300:
    np = *nfac + 1;
    ii = 1;
    n1 = *nfac - 1;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y = dir[i__];
	if (dh[ii] == 0.) {
	    goto L2302;
	}
	ij = ii;
	ip = i__ + 1;
	i__2 = *nfac;
	for (j = ip; j <= i__2; ++j) {
	    ++ij;
	    y += dir[j] * dh[ij];
	}
L2302:
	dir[i__] = y * dh[ii];
	ii = ii + np - i__;
    }
    dir[*nfac] *= dh[ii];
    i__1 = n1;
    for (k = 1; k <= i__1; ++k) {
	i__ = *nfac - k;
	ii = ii - np + i__;
	if (dir[i__] == 0.) {
	    goto L2311;
	}
	ip = i__ + 1;
	ij = ii;
	y = dir[i__];
	i__2 = *nfac;
	for (j = ip; j <= i__2; ++j) {
	    ++ij;
	    dir[j] += dh[ij] * dir[i__];
	}
    }
L2311:
L2312:
    nfac1 = *nfac + 1;
    n2fac = *nfac * nfac1 / 2;
    nnfac = *n - *nfac;
    k = n2fac;
    if (*nfac == *n) {
	goto L268;
    }
    i__1 = *n;
    for (i__ = nfac1; i__ <= i__1; ++i__) {
	dir[i__] = 0.;
    }
    if (*nfac == 0) {
	goto L265;
    }
    i__1 = *nfac;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = nfac1; j <= i__2; ++j) {
	    ++k;
	    if (x2[j] == 0.f) {
		goto L260;
	    }
	    dir[i__] += dh[k] * x2[j];
	}
    }
L260:
/*     calcul autres comp de dh*s=d en deux fois */
    k = n2fac;
    i__1 = *nfac;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = nfac1; i__ <= i__2; ++i__) {
	    ++k;
	    dir[i__] += dh[k] * x2[j];
	}
    }
L265:
    k = n2fac + *nfac * nnfac;
    i__1 = *n;
    for (j = nfac1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    ++k;
	    if (x2[j] == 0.f) {
		goto L266;
	    }
	    dir[i__] += dh[k] * x2[j];
	}
    }
L266:
    if (*nfac == *n - 1) {
	goto L268;
    }
    nm1 = *n - 1;
    k = n2fac + *nfac * nnfac;
    i__1 = nm1;
    for (i__ = nfac1; i__ <= i__1; ++i__) {
	++k;
	i1 = i__ + 1;
	i__2 = *n;
	for (j = i1; j <= i__2; ++j) {
	    ++k;
	    if (x2[j] == 0.f) {
		goto L267;
	    }
	    dir[i__] += dh[k] * x2[j];
	}
L267:
	;
    }
/*     calcul de dh*s fini */
/*     calcul sig1 pour 2eme mise a jour */
L268:
    sig1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sig1 += dir[i__] * x2[i__];
    }
    if (sig1 > 0.) {
	goto L272;
    }

/*     ****************************************************** */
    *indqn = 8;
    if (iter == 1) {
	*indqn = -5;
    }
L272:
    sig1 = -1. / sig1;
/*     truc powell si (y,s) negatif */
    if (cof1 > *zero) {
	goto L277;
    }
    teta = -1. / sig1;
    teta = teta * .8f / (teta - cof1);
    teta1 = 1. - teta;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g1[i__] = teta * g1[i__] + teta1 * dir[i__];
    }
    cof1 = -.2f / sig1;
L277:

/*      premiere mise a jour de dh */
    sig = 1. / cof1;
    zsig1 = 1. / sig1;
    mk = 0;
    ir = *nfac;
    epsmc = 1e-9;
    calmaj_(&dh[1], n, &g1[1], &sig, &x2[1], &ir, &mk, &epsmc, nfac);
    if (ir != *nfac) {
	goto L280;
    }
    calmaj_(&dh[1], n, &dir[1], &sig1, &x2[1], &ir, &mk, &epsmc, nfac);
    if (ir != *nfac) {
	goto L280;
    }
    goto L300;
L280:
    *indqn = 8;
    if (iter == 1) {
	*indqn = -5;
    }
    return 0;
L300:

/*     section 3 determination des variables libres et bloquees */

/*     calcul eps1 */

    scal1 = scal;
    if (*ieps1 == 1) {
	scal1 = 0.;
    }
    if (*ieps1 == 2) {
	scal1 = scal * cscal1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x1[i__] = x[i__] - scal1 * (d__1 = g[i__], fabs(d__1)) * g[i__];
    }
    proj_(n, &binf[1], &bsup[1], &x1[1]);
    eps1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eps1 += (d__1 = x1[i__] - x[i__], fabs(d__1));
    }
    eps1 = min(eps0,eps1);
    if (*ieps1 == 1) {
	eps1 = 0.;
    }
    if (*ieps1 == 2) {
	eps1 *= 1e4;
    }
    ifac = 0;
    idfac = 0;
    k = 0;


    gr = 0.;
    if (*ig == 1) {
	gr = difg * .2f / *n;
    }
    n3 = *n;
    if (*in == 1) {
	n3 = *n / 10;
    }
/*     si irit=1 on peut relacher des variables */
    irit = 0;
    if (difg1 <= *epsrel * difg0) {
	irit = 1;
    }
    if (*irel == 0 || iter == 1) {
	irit = 1;
    }

    tiers = .33333333333333331;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	--izig[k];
	if (izig[k] <= 0) {
	    izig[k] = 0;
	}
	bi = binf[k];
	bs = bsup[k];
	ic = indic[k];
	d1 = x[k] - bi;
	d2 = bs - x[k];
	dd = (bs - bi) * tiers;
	ep = min(eps1,dd);
	if (d1 > ep) {
	    goto L324;
	}
	if (g[k] > 0.f) {
	    goto L330;
	}
	goto L335;
L324:
	if (d2 > ep) {
	    goto L335;
	}
	if (g[k] > 0.f) {
	    goto L335;
	}
	goto L330;
/*     on defactorise si necessaire */
L330:
	if (ic > *nfac) {
	    goto L340;
	}
	++idfac;
	mode = -1;
	izig[k] += *izag;
	ajour_(&mode, n, &k, nfac, &dh[1], &x2[1], &indic[1]);
	if (mode == 0) {
	    goto L340;
	}
	*indqn = 8;
	if (iter == 1) {
	    *indqn = -5;
	}
	return 0;
/*     on factorise */
L335:
	if (irit == 0) {
	    goto L340;
	}
	if (ic <= *nfac) {
	    goto L340;
	}
	if (izig[k] >= 1) {
	    goto L340;
	}
	mode = 1;
	if (ifac >= n3 && iter > 1) {
	    goto L340;
	}
	if ((d__1 = g[k], fabs(d__1)) <= gr) {
	    goto L340;
	}
	++ifac;
	ajour_(&mode, n, &k, nfac, &dh[1], &x2[1], &indic[1]);
	if (mode == 0) {
	    goto L340;
	}
	*indqn = 8;
	if (iter == 1) {
	    *indqn = -5;
	}
	return 0;
    }
L340:

/*     *********************************************** a voir */
    if (iconv == 1) {
	return 0;
    }

/*     section 6 resolution systeme lineaire et expression de dir */
/*     on inverse le syst correspondant aux nl premieres composantes */
/*     dans le nouveau syst d indices */
    ir = *nfac;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = indic[i__];
	x2[i1] = g[i__];
    }
    if (ir < *nfac) {
	goto L412;
    }
    if (*nfac > 1) {
	goto L400;
    }
    x2[1] /= dh[1];
    goto L412;
L400:
    i__1 = *nfac;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ij = i__;
	i1 = i__ - 1;
	v = x2[i__];
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    v -= dh[ij] * x2[j];
	    ij = ij + *nfac - j;
	}
	x2[i__] = v;
	x2[i__] = v;
    }
    x2[*nfac] /= dh[ij];
    np = *nfac + 1;
    i__1 = *nfac;
    for (nip = 2; nip <= i__1; ++nip) {
	i__ = np - nip;
	ii = ij - nip;
	v = x2[i__] / dh[ii];
	ip = i__ + 1;
	ij = ii;
	i__2 = *nfac;
	for (j = ip; j <= i__2; ++j) {
	    ++ii;
	    v -= dh[ii] * x2[j];
	}
	x2[i__] = v;
    }
L412:
    if (ir == *nfac) {
	goto L660;
    }
    *indqn = 7;
    if (iter == 1) {
	*indqn = -6;
    }
    return 0;
L660:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = indic[i__];
	dir[i__] = -g[i__] * scal;
	if (i1 <= *nfac) {
	    dir[i__] = -x2[i1];
	}
    }

/*     gestion contraintes actives (si iact=1) */
    if (*iact != 1) {
	goto L675;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (izig[i__] > 0) {
	    dir[i__] = 0.f;
	}
	if (indic[i__] > *nfac) {
	    dir[i__] = 0.;
	}
    }
L675:

/*     recherche lineaire */
/*     conservation de x et g . calcul de dir+ et fpn */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g1[i__] = g[i__];
	x1[i__] = x[i__];
    }
/*     ifp =1 si fpn trop petit. on prend alors d=-g */
    ifp = 0;
    fn = *f;
L709:
    fpn = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] - binf[i__] <= epsx[i__] && dir[i__] < 0.f) {
	    dir[i__] = 0.;
	}
	if (bsup[i__] - x[i__] <= epsx[i__] && dir[i__] > 0.f) {
	    dir[i__] = 0.;
	}
	fpn += g[i__] * dir[i__];
    }
    if (fpn > 0.) {
	if (ifp == 1) {
	    *indqn = 6;
	    if (iter == 1) {
		*indqn = -3;
	    }
	    return 0;
	} else {
	    ifp = 1;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (izig[i__] > 0) {
		    dir[i__] = -scal * g[i__];
		}
	    }
	    irit = 1;
	    goto L709;
	}
    }
/*     calcul du t initial suivant une idee de fletcher */
    t1 = t;
    if (iter == 1) {
	diff = *df0;
    }
    t = diff * -2. / fpn;
    if (t > .3 && t < 3.) {
	t = 1.;
    }
    if (eps1 < eps0) {
	t = 1.;
    }
    if (*indqn == 2) {
	t = 1.;
    }
    if (iter > 1 && t1 > .01 && t1 < 100.) {
	t = 1.;
    }
    tmax = 1e10;
    t = min(t,tmax);
/* Computing MAX */
    d__1 = t, d__2 = *zero * 10.f;
    t = max(d__1,d__2);
/*     amd,amf tests sur h'(t) et diff */
    amd = .7f;
    amf = .1f;
    napm = 15;
    napm1 = nap + napm;
    if (napm1 > *napmax) {
	napm1 = *napmax;
    }
    rlbd_(&indrl, n, (S_fp)simul, &x[1], &binf[1], &bsup[1], &fn, &fpn, &t, &
	    tmax, &dir[1], &g[1], &tproj, &amd, &amf, zero, &nap, &napm1, &x2[
	    1], &izs[1], &rzs[1], &dzs[1]);
    if (indrl >= 10) {
	indsim = 4;
	++nap;
	(*simul)(&indsim, n, &x[1], f, &g[1], &izs[1], &rzs[1], &dzs[1]);
	if (indsim <= 0) {
	    *indqn = -3;
	    if (indsim == 0) {
		*indqn = 0;
	    }
	    return 0;
	}
    }
    if (indrl <= 0) {
	*indqn = 10;
	if (indrl == 0) {
	    *indqn = 0;
	}
	if (indrl == -3) {
	    *indqn = 13;
	}
	if (indrl == -4) {
	    *indqn = 12;
	}
	if (indrl <= -1000) {
	    *indqn = 11;
	}
	return 0;
    }

    if (nap < *napmax) {
	goto L758;
    }
    *f = fn;
    *indqn = 4;
    return 0;
L758:
/*     section 8 test de convergence */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = x[i__] - x1[i__], fabs(d__1)) > epsx[i__]) {
	    goto L806;
	}
    }
    *f = fn;
    *indqn = 3;
    return 0;
L806:
    difg = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	aa = g[i__];
	if (x[i__] - binf[i__] <= epsx[i__]) {
	    aa = min(0.,aa);
	}
	if (bsup[i__] - x[i__] <= epsx[i__]) {
	    aa = max(0.,aa);
	}
/* Computing 2nd power */
	d__1 = aa;
	difg += d__1 * d__1;
    }
    difg1 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (indic[i__] > *nfac) {
	    goto L820;
	}
	aa = g[i__];
	if (x[i__] - binf[i__] <= epsx[i__]) {
	    aa = min(0.,aa);
	}
	if (bsup[i__] - x[i__] <= epsx[i__]) {
	    aa = max(0.,aa);
	}
/* Computing 2nd power */
	d__1 = aa;
	difg1 += d__1 * d__1;
    }
L820:
    difg1 = sqrt(difg1);
    difg = sqrt(difg);
    difg /= sqrt((float) (*n));
    diff = (d__1 = *f - fn, fabs(d__1));
    *df0 = -diff;
    if (irit == 1) {
	difg0 = difg1;
    }
    *f = fn;
    if (diff < *epsf) {
	*indqn = 2;
	return 0;
    }
    if (difg > *epsg) {
	goto L200;
    }
    *indqn = 1;
    return 0;
} /* zqnbd_ */

