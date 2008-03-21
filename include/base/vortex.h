c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C                              vortex.h
C***********************************************************************
C
C  METHODE DES VORTEX POUR CONDITIONS AUX LIMITES D'ENTREE EN L.E.S.
C
C --------------
C PARAMETRES MAX
C --------------
      INTEGER    NENTMX, NDATMX
      PARAMETER (NENTMX = 10)
      PARAMETER (NDATMX = 10000)
C
C NENTMX    : NOMBRE D'ENTREE MAX
C NDATMX    : NOMBRE DE POINTS MAX POUR LE FICHIER DES DONNEES
C
C ----------
C DIMENSIONS
C ----------
      INTEGER         ICVOR(NENTMX)   , ICVOR2(NENTMX)  ,
     &                ICVMAX , NVOMAX
C
      COMMON /IDIMVO/ ICVOR  , ICVOR2 , ICVMAX , NVOMAX
C
C ICVOR  : NOMBRE DE FACES (GLOBAL) UTILISANT DES VORTEX
C          POUR CHAQUE ENTREE
C ICVOR2 : COMPTEUR DU NOMBRE LOCAL DE FACES UTILISANT DES VORTEX
C ICVMAX : NOMBRE MAX DE FACES UTILISANT DES VORTEX (SUR TOUTES ENTREES
C          CONFONDUES)
C NVOMAX : NOMBRE MAX DE VORTEX UTILISE (TOUTES ENTREES CONFONDUES)
C
C ---------
C POINTEURS
C ---------
C
      INTEGER         IIREPV , IIFAGL , IIVRCE ,
     &                IXYZV  , IVISV  ,
     &                IYZCEL , IUVORT , IVVORT , IWVORT ,
     &                IYZVOR , IYZVOA , ISIGNV , IXSIGM ,
     &                IXGAMM , IXTMP  , IXTMPL ,
     &                IW1X   , IW1Y   , IW1Z   , IW1V   ,
     &                IW2X   , IW2Y   , IW2Z   , IW2V
C
      COMMON /IIVORT/ IIREPV , IIFAGL , IIVRCE ,
     &                IXYZV  , IVISV  ,
     &                IYZCEL , IUVORT , IVVORT , IWVORT ,
     &                IYZVOR , IYZVOA , ISIGNV , IXSIGM ,
     &                IXGAMM , IXTMP  , IXTMPL ,
     &                IW1X   , IW1Y   , IW1Z   , IW1V   ,
     &                IW2X   , IW2Y   , IW2Z   , IW2V
C
C IIREPV    : DEBUT DU TABLEAU ASSOCIANT AUX FACES DE BORD
C             LE NUMERO D'UNE ENTREE
C IIFAGL    : DEBUT DU TABLEAU DE CONNECTIVITE
C IIVRCE    : DEBUT DU TABLEAU REPERANT LA CELLULE LA PLUS VOISINE
C             DE CHAQUE VORTEX
C IXYZV     : DEBUT DU TABLEAUX CONTENANT LES COORDONNEES DE
C             TOUTES LES FACES D'ENTREE
C IVISV     : DEBUT DU TABLEAU CONTENANT LA VISCOSITE SUR
C             TOUTES LES FACES D'ENTREE
C IYZCEL    : DEBUT DU TABLEAU CONTENANT LES COORDONNEES DES
C             FACES D'ENTREE DANS LE REPERE LOCAL
C IUVORT,...: DEBUTS DES TABLEAUX CONTENANT LES COMPOSANTES DE VITESSE
C IYZVOR    : DEBUT DU TABLEAU CONTENANT LA POSITION DES VORTEX
C             DANS LE REPERE LOCAL
C IYZVOA    : DEBUT DU TABLEAU CONTENANT LA POSITION DES VORTEX
C             DANS LE REPERE LOCAL AU PAS DE TEMPS PRECEDENT
C ISIGNV    : DEBUT DU TABLEAU CONTENANT LE SENS DE ROTATION DES
C             VORTEX
C IXSIGM    : DEBUT DU TABLEAU CONTENANT LA TAILLE DES VORTEX
C IXGAMM    : DEBUT DU TABLEAU CONTENANT L'INTENSITE DES VORTEX
C IXTMP     : DEBUT DU TABLEAU CONTENANT LE TEMPS CUMULE
C IXTMPL    : DEBUT DU TABLEAU CONTENANT LE TEMPS DE VIE DES VORTEX
C IW1X,..  : DEBUT DES TABLEAUX DE TRAVAILS SERVANT A COMMUNIQUER
C             LES DONNEES AUX ENTREES A TOUS LES PROCESSEURS
C             (PLUS UTILISE APRES VORPRE)
C
C -----------------
C OPTIONS DE CALCUL
C -----------------
C
      INTEGER         NNENT  , NVORT(NENTMX)   ,
     &                INITVO(NENTMX)  ,
     &                ICAS(NENTMX)    , ITLIVO(NENTMX)  ,
     &                ISGMVO(NENTMX)  , IDEPVO(NENTMX)  ,
     &                ICLVOR(4,NENTMX), NDAT(NENTMX)
C
      COMMON /IOPTVO/ NNENT  , NVORT  , INITVO ,
     &                ICAS   , ITLIVO , ISGMVO , IDEPVO ,
     &                ICLVOR , NDAT
C
C NNENT  : NOMBRE D ENTREES UTILISEES
C NVORT  : NOMBRE DE VORTEX
C INITVO : INDICATEUR DE REINITIALISATION
C ICAS   : TYPE DE GEOMETRIE POUR L'ENTREE
C ITLIVO : TYPE DE MODELE POUR LA DUREE DE VIE
C ISGMVO : TYPE DE MODELE POUR LA TAILLE DES VORTEX
C IDEPVO : TYPE DE MODELE POUR LA MARCHE EN TEMPS
C ICLVOR : TYPE DE CONDITION AUX LIMITES
C NDAT   : NOMBRE DE LIGNES DU FICHIER DE DONNEES
C
C -------
C DONNEES
C -------
C
      DOUBLE PRECISION TLIMVO(NENTMX), XSGMVO(NENTMX), UD(NENTMX),
     &                 XDAT(NDATMX,NENTMX),
     &                 YDAT(NDATMX,NENTMX), ZDAT(NDATMX,NENTMX),
     &                 UDAT(NDATMX,NENTMX),
     &                 VDAT(NDATMX,NENTMX), WDAT(NDATMX,NENTMX),
     &                 DUDAT(NDATMX,NENTMX),
     &                 KDAT(NDATMX,NENTMX), EPSDAT(NDATMX,NENTMX),
     &                 UDEBIT(NENTMX), KDEBIT(NENTMX), EDEBIT(NENTMX),
     &                 DIR1(3,NENTMX), DIR2(3,NENTMX), DIR3(3,NENTMX),
     &                 CEN(3,NENTMX) , SURF(3,NENTMX),
     &                 YMAX(NENTMX)  , YMIN(NENTMX),
     &                 ZMAX(NENTMX)  , ZMIN(NENTMX),
     &                 XSURFV(NENTMX), LLZ(NENTMX),
     &                 LLY(NENTMX)   , LLD(NENTMX)
C
      COMMON /ROPTVO/  TLIMVO , XSGMVO , UD     ,
     &                 XDAT   , YDAT   , ZDAT   ,
     &                 UDAT   , VDAT   , WDAT   ,
     &                 DUDAT  , KDAT   ,
     &                 EPSDAT , UDEBIT , KDEBIT , EDEBIT ,
     &                 DIR1   , DIR2   , DIR3   , CEN    , SURF   ,
     &                 YMAX   , YMIN   , ZMAX   , ZMIN   ,
     &                 XSURFV , LLZ    , LLY    , LLD
C
C
C TLIMVO      : TEMPS DE VIE MAX DES VORTEX IMPOSE PAR L'UTILISATEUR
C XSGMVO      : DIAMETRE DES VORTEX IMPOSE PAR L'UTILISATEUR
C UD          : VITESSE DE DEPLACEMENT (MAX) IMPOSEE PAR L'UTILISATEUR
C XDAT, ...   : COORDONNEES DES POINTS OU SONT CONNUES LES DONNEES
C UDAT        : VITESSE MOYENNE PRINCIPALE (FICHIER DE DONNEES)
C VDAT,WDAT   : VITESSE MOYENNE TRANSVERSE (FICHIER DE DONNEES)
C DUDAT       : DERIVE NORMALE DE LA VITESSE PRINCIPALE (FICHIER D'ENTREE)
C KDAT        : EC MOYENNE (FICHIER D'ENTREE)
C EPSDAT      : DISSIPATION (FICHIER D'ENTREE)
C UDEBIT      : VITESSE MOYENNE IMPOSEE PAR L'UTILISATEUR EN ENTREE
C KDEBIT      : EC IMPOSEE PAR L'UTILISATEUR EN ENTREE
C EDEBIT      : DISSIPATION IMPOSEE PAR L'UTILISATEUR EN ENTREE
C DIR1,...    : VECTEURS DEFINISSANT LE REPERE LOCAL DANS LE PLAN D'ENTREE
C CEN         : COORDONNEES DU CENTRE DE L'ENTREE
C SURF        : VECTEUR SURFACE DU PLAN D'ENTREE (SUPPOSEE PLANE)
C XMAX,...    : DIMENSIONS MAX DE L'ENTREE DANS LE REPERE LOCAL
C LLZ,LLY,LLD : DIMENSIONS DE L'ENTREE DANS LE CALCUL
C
      CHARACTER*50     FICVOR(NENTMX)
      COMMON /COPTVO/  FICVOR

C FICVOR : NOM DU FICHIER DE DONNEE
