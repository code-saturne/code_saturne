!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

!                              vortex.h
!===============================================================================

!  METHODE DES VORTEX POUR CONDITIONS AUX LIMITES D'ENTREE EN L.E.S.

! --------------
! PARAMETRES MAX
! --------------
integer    nentmx, ndatmx
parameter (nentmx = 10)
parameter (ndatmx = 10000)

! NENTMX    : NOMBRE D'ENTREE MAX
! NDATMX    : NOMBRE DE POINTS MAX POUR LE FICHIER DES DONNEES

! ----------
! DIMENSIONS
! ----------
integer         icvor(nentmx)   , icvor2(nentmx)  ,               &
                icvmax , nvomax

common /idimvo/ icvor  , icvor2 , icvmax , nvomax

! ICVOR  : NOMBRE DE FACES (GLOBAL) UTILISANT DES VORTEX
!          POUR CHAQUE ENTREE
! ICVOR2 : COMPTEUR DU NOMBRE LOCAL DE FACES UTILISANT DES VORTEX
! ICVMAX : NOMBRE MAX DE FACES UTILISANT DES VORTEX (SUR TOUTES ENTREES
!          CONFONDUES)
! NVOMAX : NOMBRE MAX DE VORTEX UTILISE (TOUTES ENTREES CONFONDUES)

! ---------
! POINTEURS
! ---------

integer         iirepv , iifagl , iivrce ,                        &
                ixyzv  , ivisv  ,                                 &
                iyzcel , iuvort , ivvort , iwvort ,               &
                iyzvor , iyzvoa , isignv , ixsigm ,               &
                ixgamm , ixtmp  , ixtmpl ,                        &
                iw1x   , iw1y   , iw1z   , iw1v   ,               &
                iw2x   , iw2y   , iw2z   , iw2v

common /iivort/ iirepv , iifagl , iivrce ,                        &
                ixyzv  , ivisv  ,                                 &
                iyzcel , iuvort , ivvort , iwvort ,               &
                iyzvor , iyzvoa , isignv , ixsigm ,               &
                ixgamm , ixtmp  , ixtmpl ,                        &
                iw1x   , iw1y   , iw1z   , iw1v   ,               &
                iw2x   , iw2y   , iw2z   , iw2v

! IIREPV    : DEBUT DU TABLEAU ASSOCIANT AUX FACES DE BORD
!             LE NUMERO D'UNE ENTREE
! IIFAGL    : DEBUT DU TABLEAU DE CONNECTIVITE
! IIVRCE    : DEBUT DU TABLEAU REPERANT LA CELLULE LA PLUS VOISINE
!             DE CHAQUE VORTEX
! IXYZV     : DEBUT DU TABLEAUX CONTENANT LES COORDONNEES DE
!             TOUTES LES FACES D'ENTREE
! IVISV     : DEBUT DU TABLEAU CONTENANT LA VISCOSITE SUR
!             TOUTES LES FACES D'ENTREE
! IYZCEL    : DEBUT DU TABLEAU CONTENANT LES COORDONNEES DES
!             FACES D'ENTREE DANS LE REPERE LOCAL
! IUVORT,...: DEBUTS DES TABLEAUX CONTENANT LES COMPOSANTES DE VITESSE
! IYZVOR    : DEBUT DU TABLEAU CONTENANT LA POSITION DES VORTEX
!             DANS LE REPERE LOCAL
! IYZVOA    : DEBUT DU TABLEAU CONTENANT LA POSITION DES VORTEX
!             DANS LE REPERE LOCAL AU PAS DE TEMPS PRECEDENT
! ISIGNV    : DEBUT DU TABLEAU CONTENANT LE SENS DE ROTATION DES
!             VORTEX
! IXSIGM    : DEBUT DU TABLEAU CONTENANT LA TAILLE DES VORTEX
! IXGAMM    : DEBUT DU TABLEAU CONTENANT L'INTENSITE DES VORTEX
! IXTMP     : DEBUT DU TABLEAU CONTENANT LE TEMPS CUMULE
! IXTMPL    : DEBUT DU TABLEAU CONTENANT LE TEMPS DE VIE DES VORTEX
! IW1X,..  : DEBUT DES TABLEAUX DE TRAVAILS SERVANT A COMMUNIQUER
!             LES DONNEES AUX ENTREES A TOUS LES PROCESSEURS
!             (PLUS UTILISE APRES VORPRE)

! -----------------
! OPTIONS DE CALCUL
! -----------------

integer         nnent  , nvort(nentmx)   ,                        &
                initvo(nentmx)  ,                                 &
                icas(nentmx)    , itlivo(nentmx)  ,               &
                isgmvo(nentmx)  , idepvo(nentmx)  ,               &
                iclvor(4,nentmx), ndat(nentmx)

common /ioptvo/ nnent  , nvort  , initvo ,                        &
                icas   , itlivo , isgmvo , idepvo ,               &
                iclvor , ndat

! NNENT  : NOMBRE D ENTREES UTILISEES
! NVORT  : NOMBRE DE VORTEX
! INITVO : INDICATEUR DE REINITIALISATION
! ICAS   : TYPE DE GEOMETRIE POUR L'ENTREE
! ITLIVO : TYPE DE MODELE POUR LA DUREE DE VIE
! ISGMVO : TYPE DE MODELE POUR LA TAILLE DES VORTEX
! IDEPVO : TYPE DE MODELE POUR LA MARCHE EN TEMPS
! ICLVOR : TYPE DE CONDITION AUX LIMITES
! NDAT   : NOMBRE DE LIGNES DU FICHIER DE DONNEES

! -------
! DONNEES
! -------

double precision tlimvo(nentmx), xsgmvo(nentmx), ud(nentmx),      &
                 xdat(ndatmx,nentmx),                             &
                 ydat(ndatmx,nentmx), zdat(ndatmx,nentmx),        &
                 udat(ndatmx,nentmx),                             &
                 vdat(ndatmx,nentmx), wdat(ndatmx,nentmx),        &
                 dudat(ndatmx,nentmx),                            &
                 kdat(ndatmx,nentmx), epsdat(ndatmx,nentmx),      &
                 udebit(nentmx), kdebit(nentmx), edebit(nentmx),  &
                 dir1(3,nentmx), dir2(3,nentmx), dir3(3,nentmx),  &
                 cen(3,nentmx) , surf(3,nentmx),                  &
                 ymax(nentmx)  , ymin(nentmx),                    &
                 zmax(nentmx)  , zmin(nentmx),                    &
                 xsurfv(nentmx), llz(nentmx),                     &
                 lly(nentmx)   , lld(nentmx)

common /roptvo/  tlimvo , xsgmvo , ud     ,                       &
                 xdat   , ydat   , zdat   ,                       &
                 udat   , vdat   , wdat   ,                       &
                 dudat  , kdat   ,                                &
                 epsdat , udebit , kdebit , edebit ,              &
                 dir1   , dir2   , dir3   , cen    , surf   ,     &
                 ymax   , ymin   , zmax   , zmin   ,              &
                 xsurfv , llz    , lly    , lld


! TLIMVO      : TEMPS DE VIE MAX DES VORTEX IMPOSE PAR L'UTILISATEUR
! XSGMVO      : DIAMETRE DES VORTEX IMPOSE PAR L'UTILISATEUR
! UD          : VITESSE DE DEPLACEMENT (MAX) IMPOSEE PAR L'UTILISATEUR
! XDAT, ...   : COORDONNEES DES POINTS OU SONT CONNUES LES DONNEES
! UDAT        : VITESSE MOYENNE PRINCIPALE (FICHIER DE DONNEES)
! VDAT,WDAT   : VITESSE MOYENNE TRANSVERSE (FICHIER DE DONNEES)
! DUDAT       : DERIVE NORMALE DE LA VITESSE PRINCIPALE (FICHIER D'ENTREE)
! KDAT        : EC MOYENNE (FICHIER D'ENTREE)
! EPSDAT      : DISSIPATION (FICHIER D'ENTREE)
! UDEBIT      : VITESSE MOYENNE IMPOSEE PAR L'UTILISATEUR EN ENTREE
! KDEBIT      : EC IMPOSEE PAR L'UTILISATEUR EN ENTREE
! EDEBIT      : DISSIPATION IMPOSEE PAR L'UTILISATEUR EN ENTREE
! DIR1,...    : VECTEURS DEFINISSANT LE REPERE LOCAL DANS LE PLAN D'ENTREE
! CEN         : COORDONNEES DU CENTRE DE L'ENTREE
! SURF        : VECTEUR SURFACE DU PLAN D'ENTREE (SUPPOSEE PLANE)
! XMAX,...    : DIMENSIONS MAX DE L'ENTREE DANS LE REPERE LOCAL
! LLZ,LLY,LLD : DIMENSIONS DE L'ENTREE DANS LE CALCUL

character*50     ficvor(nentmx)
common /coptvo/  ficvor

! FICVOR : NOM DU FICHIER DE DONNEE
