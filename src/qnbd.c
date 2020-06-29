/* qnbd.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/*     removed io and comment out any printing from fortran. */
/*     Scilab ( http://www.scilab.org/ ) - This file is part of Scilab */
/*     Copyright (C) 1986 - INRIA - F. BONNANS */

/* Copyright (C) 2012 - 2016 - Scilab Enterprises */

/* This file is hereby licensed under the terms of the GNU GPL v2.0, */
/* pursuant to article 5.3.4 of the CeCILL v.2.1. */
/* This file was originally licensed under the terms of the CeCILL v2.1, */
/* and continues to be available under such terms. */
/* For more information, see the COPYING file which you should have received */
/* along with this program. */

// From f2c header
typedef int /* Unknown procedure type */ (*U_fp)();

/* Subroutine */ int qnbd_(int *indqn, U_fp simul, int *n, double 
	*x, double *f, double *g, double *zero, int *napmax, 
	int *itmax, double *epsf, double *epsg, double *epsx, 
	double *df0, double *binf, double *bsup, int *nfac, 
	double *trav, int *ntrav, int *itrav, int *nitrav, 
	int *izs, float *rzs, double *dzs)
{
    static int n1, n2, n3, n4, n5, ig, in, ni1, ni2, iact, izag, irel, 
	    ieps1;
    extern /* Subroutine */ int zqnbd_(int *, U_fp, double *, int 
	    *, double *, double *, double *, double *, 
	    double *, double *, int *, int *, int *, 
	    int *, int *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, int *, int *, int *, int *, int 
	    *, double *, int *, int *, float *, double *);
    static double epsrel;

/* !but */
/*     code de minimisation d une fonction reguliere sous contraintes */
/*     de bornes , aux normes modulopt */
/* !methode */
/*     principe de l algorithme : quasi-newton + projection */
/*     details dans le rapport inria n. 242,1983 */
/*     cette version permet de tester plusieurs variantes de l algorithme */
/*     en modifiant certains parametres internes (cf comment dans le code) */
/*     taille memoire necessaire de l ordre de n**2/2 */
/*     pour les problemes de grande taille le code gcbd est mieux adapte */

/* !sous programmes appeles */
/*     zqnbd   optimiseur effectif */
/*     proj    projection */
/*     calmaj  mise a jour du hessien */
/*     ajour mise a jour des facteurs de choleski */
/*     rlbd,satur   recherche lineaire avec bornes */

/* !liste d'appel */

/*     subroutine qnbd(indqn,simul,n,x,f,g,iprint,io,zero, */
/*    & napmax,itmax,epsf,epsg,epsx,df0,binf,bsup,nfac, */
/*    & trav,ntrav,itrav,nitrav,izs,rzs,dzs) */

/*     indqn   indicateur de qnbd                                  es */
/*       en entree =1 standard */
/*                 =2 dh et indic initialises au debut de trav et itrav */
/*                    ifac,f,g initialises */
/*       en sortie */
/*        si < 0 incapacite de calculer un point meilleur que le point initial */
/*        si = 0 arret demande par l utilisateur */
/*        si > 0 on fournit un point meilleur que le point de depart */
/*        < -10 parametres d entree non convenables */
/*        = -6  arret lors du calcul de la direction de descente et iter=1 */
/*        = -5  arret lors du calcul de l approximation du hessien  iter=1 */
/*        = -3  anomalie de simul : indic negatif en un point ou */
/*              f et g ont ete precedemment calcules */
/*        = -2  echec de la recherche lineaire a la premiere iteration */
/*        = -1  f non definie au point initial */
/*        =  1  arret sur epsg */
/*        =  2            epsf */
/*        =  3            epsx */
/*        =  4            napmax */
/*        =  5            itmax */
/*        =  6  pente dans la direction opposee au gradient trop petite */
/*        =  7  arret lors du calcul de la direction de descente */
/*        =  8  arret lors du calcul de l approximation du hessien */
/*        = 10  arret par echec de la recherche lineaire , cause non precisee */
/*        = 11  idem avec indsim < 0 */
/*        = 12            un pas trop petit proche d un pas trop grand */
/*                        ceci peut resulter d une erreur dans le gradient */
/*        = 13            trop grand nombre d appels dans une recherche lineaire */
/*     simul voir les normes modulopt */
/*     n dim de x                                                 e */
/*     binf,bsup bornes inf,sup,de dim n                          e */
/*     x variables a optimiser (controle)                          es */
/*     f valeur du critere                                         s */
/*     g gradient de f                                             s */
/*     zero  proche zero machine                                             e */
/*     napmax nombre maximum d appels de simul                               e */
/*     itmax nombre maximum d iterations de descente               e */
/*     itrav vect travail dim nitrav=2n , se decompose en indic et izig */
/*     nfac nombre de variables factorisees                  (e si indqn=2)  s */
/*     iprint facteur d impression                                              e */
/*     varie de 0 (pas d impressions) a 3 (nombreuses impressions) */
/*     io numero du fichier de resultats                                     e */
/*     epsx vect dim n precision sur x                                       e */
/*     epsf critere arret sur f                                              e */
/*     epsg arret si sup a norm2(g+)/n                                       e */
/*     trav vect travail dim ntrav */
/*     il faut ntrav > n(n+1)/2 +6n */
/*     df0>0 decroissance f prevue     (prendre 1. par defaut)               e */
/*     izs,rzs,dzs : cf normes modulopt                                     es */
/* ! */
/*     indications sur les variables internes a qnbd et zqnbd */
/*     izig  sert a la memorisation des contraintes (actif si izag>1) */
/*     si i ne change pas d ens on enleve 1 a izig  (positif) */
/*     sinon on ajoute izag */
/*     factorisation seulement si izig est nul */
/*     dh estimation hessien dim n(n+1)/2 rangee en troismorceaux es */
/*     indic(i) nouvel indice de l indice i */
/*     indic vect dim n ordre de rangement des indices                       es */
/*     pas necessaire de l initialiser si indqn=1 */

/*     parametres de la recherche lineaire */
/*     amd,amf param. du test de wolfe .    (.7,.1) */
/*     napm nombre max d appels dans la rl  (=15) */


/* ---- initial printing */
/* $$$      if(iprint.ge.1) then */
/* $$$         call basout(io_out, io, '') */
/* $$$         write(bufstr,1010) */
/* $$$         call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$         write(bufstr,750) n,epsg,iprint */
/* $$$         call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$         write(bufstr,751) itmax */
/* $$$         call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$         write(bufstr,752) napmax */
/* $$$         call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$         call basout(io_out ,io , */
/* $$$     $    '------------------------------------------------') */
/* $$$1010    format(' *********** qnbd (with bound cstr) ****************') */
/* $$$750     format('dimension=',i10,', epsq=',e24.16, */
/* $$$     $ ', verbosity level: iprint=',i10) */
/* $$$751     format('max number of iterations allowed: iter=',i10) */
/* $$$752     format('max number of calls to costf allowed: nap=',i10) */
/* $$$      endif */


/*     parametres caracteristiques de l algorithme */
/*     si les parametres sont nuls l algorithme est celui du rr 242 */
/*     ig=1 test sur grad(i) pour relach var */
/*     in=1 limite le nombre de factorisations par iter a n/10 */
/*     irel=1 test sur decroissance grad pour rel a iter courante */
/*     epsrel taux de decroissance permettant le relachement (cf irit) */
/*     iact blocage variables dans ib (gestion contraintes actives) */
/*     ieps1 =1 eps1 egal a zero */
/*     note eps1 correspond a eps(xk) */
    /* Parameter adjustments */
    --bsup;
    --binf;
    --epsx;
    --g;
    --x;
    --trav;
    --itrav;
    --izs;
    --rzs;
    --dzs;

    /* Function Body */
    ig = 0;
    in = 0;
    irel = 1;
    epsrel = .5;
    izag = 0;
    iact = 1;
    ieps1 = 0;

/*     decoupage du vecteur trav */
    n1 = *n * (*n + 1) / 2 + 1;
    n2 = n1 + *n;
    n3 = n2 + *n;
    n4 = n3 + *n;
    n5 = n4 + *n - 1;
    if (*ntrav < n5) {
/* $$$         if(iprint.gt.0) then */
/* $$$           write(bufstr,110)ntrav,n5 */
/* $$$           call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$           endif */
/* $$$110      format(' qnbd : ntrav=',i8,' devrait valoir ',i8) */
	*indqn = -11;
	return 0;
    }
    ni1 = *n + 1;
    if (*nitrav < *n << 1) {
	ni2 = *n << 1;
/* $$$         if(iprint.gt.0) then */
/* $$$           write(bufstr,111)nitrav,ni2 */
/* $$$           call basout(io_out ,io ,bufstr(1:lnblnk(bufstr))) */
/* $$$           endif */
/* $$$111      format(' qnbd : nitrav=',i8,'devrait valoir',i8) */
	*indqn = -12;
	return 0;
    }
    zqnbd_(indqn, (U_fp)simul, &trav[1], n, &binf[1], &bsup[1], &x[1], f, &g[
	    1], zero, napmax, itmax, &itrav[1], &itrav[ni1], nfac, &epsx[1], 
	    epsf, epsg, &trav[n1], &trav[n2], &trav[n3], &trav[n4], df0, &ig, 
	    &in, &irel, &izag, &iact, &epsrel, &ieps1, &izs[1], &rzs[1], &dzs[
	    1]);
    return 0;
} /* qnbd_ */

