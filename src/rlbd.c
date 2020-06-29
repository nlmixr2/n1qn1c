/* rlbd.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/*     Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/* Copyright (C) 1988 - INRIA - F. BONNANS */

/* Copyright (C) 2012 - 2016 - Scilab Enterprises */

/* This file is hereby licensed under the terms of the GNU GPL v2.0, */
/* pursuant to article 5.3.4 of the CeCILL v.2.1. */
/* This file was originally licensed under the terms of the CeCILL v2.1, */
/* and continues to be available under such terms. */
/* For more information, see the COPYING file which you should have received */
/* along with this program. */

#include <stdio.h>
#include <float.h>
#include <math.h>

#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

typedef /* Subroutine */ int (*S_fp)();
typedef /* Subroutine */ int (*U_fp)();

/* Subroutine */ int rlbd_(int *indrl, int *n, S_fp simul, double 
	*x, double *binf, double *bsup, double *f, double *hp,
	 double *t, double *tmax, double *d__, double *gn, 
	double *tproj, double *amd, double *amf, double *zero,
	 int *nap, int *napmax, double *xn, int *izs, float *
	rzs, double *dzs)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, short, short);
    double sqrt(double);

    /* Local variables */
    static double a, b, e;
    static int i__, k;
    static double p, r__, a1, f0, h1, t1, t2, fa, f11, di, fn, ta, td, tg,
	     fa1, ta1, hpa, hpd, ftd, hpg, ftg, div, hpn, eps, tmi, xni;
    static int ico1;
    static double eps1;
    static char var2[3];
    static int icoi, icop;
    static double cres;
    static int icos, imax;
    static double hptd, hptg;
    extern /* Subroutine */ int proj_(int *, double *, double *, 
	    double *);
    static double epst, text, topt, hpta1;
    static int indic;
    static double difhp, extra;
    static int iproj;
    static double tmaxp, ttmin;
    extern /* Subroutine */ int satur_(int *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, int *,
	     int *, int *);
    static double extrp, ttsup, cofder;


/* !but */
/*     subroutine de recherche lineaire pour des problemes avec */
/*       contraintes de borne (traitees par projection) */
/*       le critere de retour est une extension de celui de wolfe */
/* !methode */
/*     pour chaque valeur du parametre t , sont calcules le critere */
/*     et son gradient. */
/*     une phase d extrapolation permet d obtenir un encadrement. */
/*     l intervalle est ensuite reduit suivant les cas par une methode */
/*     de dichotomie, d interpolation lineaire sur les derivees ou */
/*     d interpolation cubique. */

/* !impressions */
/*      si iprint > 2 , rlbd fournit les impressions suivantes : */

/*      la premiere ligne indique : */
/*      t      premiere valeur de t fournie en liste d' appel */
/*      tproj  plus petit t > 0 pour lequel on bute sur une borne */
/*      dh/dt  derivee en zero de h(t)=f(x+t*d)-f(x) */
/*      tmax   valeur maximale de t fournie en liste d' appel */

/*      lignes suivantes : */
/*      chaine de caracteres en debut de ligne : indique comment sera calcule */
/*      le pas de la ligne suivante ; */
/*      ic : interpolation cubique */
/*      s  : saturation d une variable sur une borne */
/*      id : interpolation lineaire sur la derivee */
/*      e  : extrapolation */
/*      d  :interpolation cubique ayant echouee t est calcule par dichotomie */
/*      b  :sauvegarde de convergence active */

/* !subroutines utilisees */
/*     proj et satur (bibl. modulopt) */
/* !liste d appel */

/*      subroutine rlbd(indrl,n,simul,proj,x,binf,bsup,f,hp,t,tmax,d,gn, */
/*     &  tproj,amd,amf,iprint,io,zero,nap,napmax,xn,izs,rzs,dzs) */

/*      e;s;e,s:parametres initialises en entree,en sortie,en entree et */
/*              en sortie */
/*      indrl<0:la recherche lineaire n a pas trouve de meilleur pas(e,s) */
/*           =0:arret demande par l'utilisateur dans simul */
/*           >0:meilleur pas fourni avec f et g */
/*           >9:meilleur pas fourni avec f et sans g */
/*           =14:deltat trop petit */
/*           =13:nap=napmax */
/*           =8:toutes les variables sont saturees */
/*           =4:deltat trop petit */
/*           =3:nap=napmax */
/*           =2:t=tmax */
/*           =1:descente serieuse avec t<tmax */
/*           =0:arret demande par l'utilisateur */
/*           =-3:nap=napmax */
/*           =-4:deltat trop petit */
/*           =-1000+indic:nap=napmax et indic<0 */
/*      n:dimension de x                                    (e,s) */
/*      simul: subroutine fournissant le critere et le gradient (e) */
/*      x:valeur initiale de la variable a optimiser en entree;valeur a */
/*        l optimum en sortie.                                     (e,s) */
/*      binf,bsup:bornes inf et sup de dimension n                 (e,s) */
/*      f:valeur du critere en x                            (e,s) */
/*      hp:derivee de f(x+t*d) par rapport a t en 0         (e) */
/*      t:pas                                               (e) */
/*      tmax:valeur maximal du pas                          (e,s) */
/*      d:direction de descente                             (e) */
/*      gn: gradient de f en xn                             (e,s) */
/*      tproj:plus petit pas saturant une nouvelle contrainte(e,s) */
/*      amf,amd:constantes du test de wolfe                 (e) */
/*      iprint<=2:pas d'impression                             (e) */
/*         >=3:une impression par calcul de simul           (e) */
/*      io:numero du fichier resultat                       (e) */
/*      zero:proche du zero machine                         (e) */
/*      nap:nombre d'appel a simul                          (e) */
/*      napmax:nombre maximum d'appel a simul               (e) */
/*      xn:tableau de travail de dimension n (=x+t*d) */
/*      izs,rzs,dzs:cf norme modulopt                       (e,s) */
/* ! */
/*       parametres de l algorithme */
/*      eps1:sauvegarde de conv.dans l interpolation lineaire sur la derivee */
/*      eps:sauvegarde de conv.dans la l interpolation par saturation */
/*          d une contrainte. */
/*      epst:sauvegerde de conv.dans l interpolation cubique */
/*      extra,extrp:servent a calculer la limite sur la variation relative */
/*      de t dans l extrapolation et l interpolation lineaire sur la derivee */
/*      cofder: intervient dans le choix entre les methodes d' interpolation */

/*        variables de travail */
/*      fn:valeur du critere en xn */
/*      hpn:derivee de f(x+t*d) par rapport a t */
/*      hpd:valeur de hpn a droite */
/*      hpg:valeur de hpn a gauche */
/*      td:pas trop grand */
/*      tg:pas trop petit */
/*      tproj:plus petit pas saturant une contrainte */
/*      tmaxp:plus grand pas saturant une contraite */
/*      ftd:valeur de f en x+td*d */
/*      ftg:valeur de f en x+tg*d */
/*      hptd:valeur de hpn en x+td*d */
/*      hptg:valeur de hpn en x+tg*d */
/*      imax=1:tmax a ete atteint */
/*          =0:tmax n a pas ete atteint */
/*      icos:indice de la variable saturee par la borne superieure */
/*      icoi:indice de la variable saturee par la borne inferieure */
/*      ico1:indice de la variable saturee en tmaxp */
/*      icop:indice de la variable saturee en tproj */


    /* Parameter adjustments */
    --xn;
    --gn;
    --d__;
    --bsup;
    --binf;
    --x;
    --izs;
    --rzs;
    --dzs;

    /* Function Body */
    *indrl = 1;
    eps1 = .9;
    eps = .1;
    epst = .1;
    extrp = 100.;
    extra = 10.;
    cofder = 100.f;
    s_copy(var2, "   ", (short)3, (short)3);

    ta1 = 0.;
    f0 = *f;
    fa1 = *f;
    hpta1 = *hp;
    imax = 0;
    hptg = *hp;
    ftg = *f;
    tg = 0.;
    td = 0.;
    icos = 0;
    icoi = 0;
    icop = 0;

/*     calcul de tproj:plus petit point de discontinuite de h'(t) */
    *tproj = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cres = d__[i__];
	if (cres < 0.) {
	    goto L4;
	} else if (cres == 0.) {
	    goto L7;
	} else {
	    goto L5;
	}
L4:
	t2 = (binf[i__] - x[i__]) / d__[i__];
	goto L6;
L5:
	t2 = (bsup[i__] - x[i__]) / d__[i__];
L6:
	if (t2 <= 0.) {
	    goto L7;
	}
	if (*tproj == 0.) {
	    *tproj = t2;
	}
	if (t2 > *tproj) {
	    goto L7;
	}
	*tproj = t2;
	icop = i__;
    }
L7:


/*              boucle */

/*     calcul de xn,de fn et de gn */
L200:
    if (*nap >= *napmax) {
	k = 3;
	goto L1000;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xn[i__] = x[i__] + *t * d__[i__];
    }
    proj_(n, &binf[1], &bsup[1], &xn[1]);
    if (icos > 0) {
	xn[icos] = bsup[icos];
    }
    if (icoi > 0) {
	xn[icoi] = binf[icoi];
    }
    indic = 4;
    (*simul)(&indic, n, &xn[1], &fn, &gn[1], &izs[1], &rzs[1], &dzs[1]);
    ++(*nap);
    if (indic < 0) {
/* $$$         if (iprint.ge.3) then */
/* $$$           write (bufstr,16000) indic,t */
/* $$$           call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$           endif */
	if (*nap >= *napmax) {
	    goto L1000;
	}
	*t = tg + (*t - tg) / 4.;
	*tmax = *t;
	imax = 1;
	icoi = 0;
	icos = 0;
	s_copy(var2, "dd ", (short)3, (short)3);
	goto L800;
    }
    if (indic == 0) {
	*indrl = 0;
	goto L1010;
    }

/*      calcul de hpg et hpd */
    hpg = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xn[i__] = x[i__] + *t * d__[i__];
    }
    if (icoi > 0) {
	xn[icoi] = bsup[icoi];
    }
    if (icos > 0) {
	xn[icos] = bsup[icos];
    }
    proj_(n, &binf[1], &bsup[1], &xn[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xni = xn[i__];
	if (binf[i__] < xni && xni < bsup[i__]) {
	    hpg += d__[i__] * gn[i__];
	}
    }
    hpd = hpg;
    if (icoi > 0) {
	hpg += d__[icoi] * gn[icoi];
    }
    if (icos > 0) {
	hpg += d__[icos] * gn[icos];
    }

    icoi = 0;
    icos = 0;
    if (hpd != 0. || hpg != 0.) {
	goto L360;
    }

/*      la derivee de h est nulle */
/*      calcul du pas saturant toutes les bornes:tmaxp */
    tmaxp = 0.;
    ico1 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cres = d__[i__];
	if (cres < 0.) {
	    goto L310;
	} else if (cres == 0.) {
	    goto L350;
	} else {
	    goto L320;
	}
L310:
	t2 = (binf[i__] - x[i__]) / d__[i__];
	goto L330;
L320:
	t2 = (bsup[i__] - x[i__]) / d__[i__];
L330:
	if (t2 <= 0.) {
	    goto L350;
	}
	if (tmaxp == 0.) {
	    tmaxp = t2;
	}
	if (tmaxp > t2) {
	    goto L350;
	}
	tmaxp = t2;
	ico1 = i__;
    }
L350:
    if (*t < tmaxp) {
	if (fn <= *f + *amf * *hp * *t) {
	    goto L1010;
	}
	*t /= 10.;
	s_copy(var2, "d  ", (short)3, (short)3);
	goto L800;
    }
    icos = ico1;
    icoi = 0;
    if (d__[ico1] < 0.) {
	icoi = ico1;
	icos = 0;
    }

/*     toutes les variables sont saturees */
/* $$$      if (iprint.ge.3) then */
/* $$$        write (bufstr,3330) tmaxp */
/* $$$        call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$        endif */
/* $$$3330   format ('toutes les variables sont saturees:tmaxp= ',e11.4) */
    *t = tmaxp;
    if (fn < *f + *amf * *hp * tmaxp) {
	*indrl = 8;
	goto L1010;
    }
    hpg = d__[ico1] * gn[ico1];
    if (fn < *f && hpg < 0.) {
	*indrl = 8;
	goto L1010;
    }
L360:

/*       test de wolfe */

    a = *f + *amf * *hp * *t;
    if (fn > a) {
/*      le pas est trop grand */
/*       (dans le cas quadratique changer eps1 et extra si td<tproj) */
	td = *t;
	t1 = *t - ta1;
	h1 = (fn - fa1) / t1;
	ftd = fn;
	hptd = hpg;
	ta = tg;
	hpn = hptd;
	hpa = hptg;
	fa = ftg;
    } else {
	if (hpd >= *amd * *hp) {
	    goto L1010;
	}
/*      le pas est trop petit */
	tg = *t;
	t1 = *t - ta1;
	h1 = (fn - fa1) / t1;
	ftg = fn;
	hptg = hpd;
	ta = td;
	hpn = hptg;
	hpa = hptd;
	fa = ftd;
	if (td == 0.) {
	    goto L700;
	}
	a1 = (d__1 = hptd / *hp, fabs(d__1));
	if (a1 > cofder && ftd > *f && hptg > *hp * .99f) {
	    hpta1 = *hp;
	    fa1 = *f;
	    ta1 = 0.;
	    goto L700;
	}
    }
    a1 = (d__1 = hpn / *hp, fabs(d__1));
    if (tg != 0. || fn <= *f || a1 <= cofder || hpn < 0.) {
	if (td <= *tproj) {
	    goto L600;
	}
	goto L500;
    }

/*       calcul du nouveau t */

/*      par interpolation lineaire sur la derivee */

    ta1 = *t;
    fa1 = fn;
    div = *hp - hptd;
    text = *t / 10.;
    if (fabs(div) > *zero) {
	text = *t * (*hp / div);
    }
    if (text > *tproj) {
	text = *t / 10.;
    }
/* Computing MAX */
    d__1 = text, d__2 = *t / (extrp * extra);
    text = max(d__1,d__2);
/* Computing MIN */
    d__1 = text, d__2 = *t * eps1;
    *t = min(d__1,d__2);
    ttsup = *t * 1.5;
    extrp = 10.f;
    if (*tproj > ta1) {
	s_copy(var2, "id ", (short)3, (short)3);
	goto L800;
    }
    ttmin = *t * .7;
    tmi = *t;
    topt = 0.;
    iproj = 0;
    satur_(n, &x[1], &binf[1], &bsup[1], &d__[1], &ttmin, &ttsup, &topt, &tg, 
	    &td, &tmi, &icoi, &icos, &iproj);
    s_copy(var2, "id ", (short)3, (short)3);
    if (topt != 0.) {
	*t = topt;
	s_copy(var2, "ids", (short)3, (short)3);
    }
    goto L800;

/*      interpolation par saturation d une contrainte */

L500:
    if (td <= *tproj) {
	goto L600;
    }
    topt = 0.;
    iproj = 1;
    ta1 = *t;
    fa1 = fn;
    ttmin = tg + eps * (td - tg);
    ttsup = td - eps * (td - tg);
    tmi = (td + tg) / 2.;
    satur_(n, &x[1], &binf[1], &bsup[1], &d__[1], &ttmin, &ttsup, &topt, &tg, 
	    &td, &tmi, &icoi, &icos, &iproj);
    if (topt == 0.) {
	goto L600;
    }
    *t = topt;
    s_copy(var2, "s  ", (short)3, (short)3);
    if (*t == ttsup || *t == ttmin) {
	s_copy(var2, "sb ", (short)3, (short)3);
    }
    goto L800;

/*      interpolation cubique */

/*      test de securite */
L600:
    if (td - tg < *zero * 100.) {
	k = 4;
	goto L1000;
    }

/*      calcul du minimum */
    b = 1.;
    p = hpn + hpa - (fn - fa) * 3. / (*t - ta);
    di = p * p - hpn * hpa;
    if (di < 0.) {
	goto L690;
    }
    if (*t - ta < 0.) {
	b = -1.;
    }
    div = hpn + p + b * sqrt(di);
    if (fabs(div) <= *zero) {
	goto L690;
    }
    r__ = hpn / div;
    topt = *t - r__ * (*t - ta);
    if (topt < tg || topt > td) {
	goto L690;
    }

/*      sauvegarde de convergence */
    e = epst * (td - tg);
    s_copy(var2, "ic ", (short)3, (short)3);
    if (topt > td - e) {
	topt = td - e;
	s_copy(var2, "icb", (short)3, (short)3);
    }
    if (topt < tg + e) {
	topt = tg + e;
	s_copy(var2, "icb", (short)3, (short)3);
    }
    ta1 = *t;
    fa1 = fn;
    *t = topt;
    goto L800;
L690:
    ta1 = *t;
    fa1 = fn;
    *t = (tg + td) * .5;
    s_copy(var2, "d  ", (short)3, (short)3);
    goto L800;

/*      extrapolation */

L700:
    if (imax >= 1) {
	k = 2;
	goto L1000;
    }
    text = *t * 10.;
    difhp = hptg - hpta1;
    if (difhp > *zero) {
	text = (*amd * *hp / 3. - hptg) * ((tg - ta1) / difhp) + tg;
	if (td != 0. && text >= td) {
	    goto L600;
	}
/*        dans le cas quadratique prendre extrp plus grand */
/* Computing MIN */
	d__1 = text, d__2 = extra * extrp * *t;
	text = min(d__1,d__2);
/* Computing MAX */
	d__1 = text, d__2 = *t * 2.5;
	text = max(d__1,d__2);
    } else {
	text = extra * extrp * *t;
    }
    ta1 = *t;
    fa1 = fn;
    hpta1 = hpn;
    extrp = 10.f;
    if (text >= *tmax / 2.) {
	text = *tmax;
	imax = 1;
    }
    if (*t < *tproj && text > *tproj) {
/* Computing MAX */
	d__1 = *tproj, d__2 = *t * 2.5;
	*t = max(d__1,d__2);
	icoi = 0;
	icos = icop;
	if (d__[icop] < 0.) {
	    icoi = icop;
	    icos = 0;
	}
	s_copy(var2, "es ", (short)3, (short)3);
	goto L800;
    }
/* Computing MIN */
    d__1 = text * 1.5;
    ttsup = min(d__1,*tmax);
    if (ttsup < *tproj) {
	goto L785;
    }
    ttmin = *t * 2;
    iproj = 0;
    tmi = text;
    topt = 0.;
    satur_(n, &x[1], &binf[1], &bsup[1], &d__[1], &ttmin, &ttsup, &topt, &tg, 
	    &td, &tmi, &icoi, &icos, &iproj);
    if (topt > 0.) {
	*t = topt;
	s_copy(var2, "es ", (short)3, (short)3);
	goto L800;
    }
L785:
    *t = text;
    s_copy(var2, "e  ", (short)3, (short)3);
L800:
    f11 = fn - *f;
/* $$$      if (iprint.ge.3.and.indic.gt.0) then */
/* $$$        write (bufstr,15000)var2,ta1,f11,hpn,h1,t1 */
/* $$$        call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$        endif */

/*      test sur deltat */
    if ((d__1 = ta1 - *t, fabs(d__1)) >= *zero * 100.) {
	goto L200;
    }
    k = 4;
/*      calcul de indrl */
L1000:
    if (indic < 0) {
	*indrl = 13;
	if (tg == 0.) {
	    *indrl = indic - 1000;
	}
	fn = ftg;
	hpn = hptg;
	*t = tg;
	goto L1010;
    }
    if (fn <= ftg) {
	*indrl = k;
	*t = tg;
	goto L1010;
    }
    if (tg == 0.) {
	*indrl = -k;
	goto L1010;
    }
    *indrl = k + 10;
    *t = tg;
    fn = ftg;
    hpn = hptg;

/*      fin du programme */
L1010:
    *f = fn;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += *t * d__[i__];
    }
    proj_(n, &binf[1], &bsup[1], &x[1]);
    if (icos > 0) {
	x[icos] = bsup[icos];
    }
    if (icoi > 0) {
	x[icoi] = binf[icoi];
    }

    if (*indrl < 0) {
	++(*nap);
	indic = 4;
	(*simul)(&indic, n, &x[1], f, &gn[1], &izs[1], &rzs[1], &dzs[1]);
    }

    t1 = *t - ta1;
    if (t1 == 0.) {
	t1 = 1.f;
    }
    h1 = (fn - fa1) / t1;
    *hp = hpd;
    f0 = *f - f0;
/* $$$      if (iprint.ge.3) then */
/* $$$        write (bufstr,15020)t,f0,hpd,h1,t1 */
/* $$$        call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$        endif */
    return 0;
} /* rlbd_ */

