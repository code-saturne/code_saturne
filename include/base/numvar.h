c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2007 EDF S.A., France
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
C                              numvar.h
C***********************************************************************
C
C POSITION DES VARIABLES
C  ( Dans RTP, RTPA )
C
C IPR                        Pression
C IU   IV   IW               Vitesse(x,y,z)
C IK                         Energie Turbulente en k-epsilon
C IR11, IR22, IR33, ...
C ... IR12, IR13, IR23       Tensions de Reynolds en Rij
C IEP                        Dissipation turbulente
C IPHI, IFB                  Variables phi et f_barre du v2f phi-model
C IOMG                       Variable Omega du k-omega SST
C ISCA(i)                    Scalaire numero i
C ISCAPP(i)                  No du scalaire physique particuliere i
C NSCAUS                     Nbre de scalaires utilisateur
C NSCAPP                     Nbre de scalaires physique particuliere
C IUMA, IVMA, IWMA           Vitesse de maillage en ALE
C
C
      INTEGER           IPR(NPHSMX) ,
     &                  IU(NPHSMX)  , IV(NPHSMX)    , IW(NPHSMX)  ,
     &                  IK(NPHSMX)  , IEP(NPHSMX)   ,
     &                  IR11(NPHSMX), IR22(NPHSMX)  , IR33(NPHSMX),
     &                  IR12(NPHSMX), IR13(NPHSMX)  , IR23(NPHSMX),
     &                  IPHI(NPHSMX), IFB(NPHSMX)   , IOMG(NPHSMX),
     &                  ISCA(NSCAMX), ISCAPP(NSCAMX),
     &                  NSCAUS      , NSCAPP        ,
     &                  IUMA        , IVMA          , IWMA
      COMMON / IPOSVR / IPR         ,
     &                  IU          , IV            , IW          ,
     &                  IK          , IEP           ,
     &                  IR11        , IR22          , IR33        ,
     &                  IR12        , IR13          , IR23        ,
     &                  IPHI        , IFB           , IOMG        ,
     &                  ISCA        , ISCAPP        ,
     &                  NSCAUS      , NSCAPP        ,
     &                  IUMA        , IVMA          , IWMA
C
C POSITION DES PROPRIETES (Physiques ou numeriques)
C  ( Dans PROPCE, PROPFA et PROPFB)
C    Le numero des proprietes est unique, quelle aue soit la
C      localisation de ces dernieres (cellule, face, face de bord)
C    Voir usclim pour quelques exemples
C
C IPPROC : Pointeurs dans PROPCE
C IPPROF : Pointeurs dans PROPFA
C IPPROB : Pointeurs dabs PROPFB
C
C IROM   : Masse volumique des phases
C IROMA  : Masse volumique des phases au pas de temps precedent
C IVISCL : Viscosite moleculaire dynamique en kg/(ms) des phases
C IVISCT : Viscosite turbulente des phases
C IVISLA : Viscosite moleculaire dynamique en kg/(ms) des phases au pas
C          de temps precedent
C IVISTA : Viscosite turbulente des phases au pas de temps precedent
C ICP    : Chaleur specifique des phases
C ICPA   : Chaleur specifique des phases au pas de temps precedent
C ITSNSA : Terme source Navier Stokes des phases au pas de temps precedent
C ITSTUA : Terme source des grandeurs turbulentes au pas de temps precedent
C ITSSCA : Terme source des scalaires au pas de temps precedent
C IESTIM : Estimateur d'erreur pour Navier-Stokes
C IFLUMA : Flux de masse associe aux variables
C IFLUAA : Flux de masse explicite (plus vu comme un tableau de travail)
C          associe aux variables
C ISMAGO : constante de Smagorinsky dynamique
C ICOUR  : Nombre de Courant des phases
C IFOUR  : Nombre de Fourier des phases
C IPRTOT : Pression totale au centre des cellules Ptot=P*+rho*g.(x-x0)
C                                                             -  - -
C IVISMA : Viscosite de maillage en ALE (eventuellement orthotrope)
C
      INTEGER           IPPROC(NPROMX), IPPROF(NPROMX), IPPROB(NPROMX),
     &                  IROM  (NPHSMX), IROMA (NPHSMX), IVISCL(NPHSMX),
     &                  IVISCT(NPHSMX), IVISLA(NPHSMX), IVISTA(NPHSMX),
     &                  ICP   (NPHSMX), ICPA  (NPHSMX), ITSNSA(NPHSMX),
     &                  ITSTUA(NPHSMX), ITSSCA(NSCAMX),
     &                  IESTIM(NESTMX,NPHSMX)         , IFLUMA(NVARMX),
     &                  IFLUAA(NVARMX), ISMAGO(NPHSMX), ICOUR (NPHSMX),
     &                  IFOUR (NPHSMX), IPRTOT(NPHSMX), IVISMA(3)
      COMMON / IPOSPP / IPPROC        , IPPROF        , IPPROB        ,
     &                  IROM          , IROMA         , IVISCL        ,
     &                  IVISCT        , IVISLA        , IVISTA        ,
     &                  ICP           , ICPA          , ITSNSA        ,
     &                  ITSTUA        , ITSSCA        ,
     &                  IESTIM                        , IFLUMA        ,
     &                  IFLUAA        , ISMAGO        , ICOUR         ,
     &                  IFOUR         , IPRTOT        , IVISMA
C
C
C POSITION DES CONDITIONS AUX LIMITES
C  (Position dans COEFA et COEFB des coef (coef. coef.f) relatifs a
C   une variable donnee)
C
C ICOEF   : coef numeros 1 (anciens coefa coefb)
C ICOEFF  : coef numeros 2 (anciens coefaf coefbf)
C ICLRTP  : pointeur dans COEFA et COEFB
C
      INTEGER           ICOEF , ICOEFF , ICLRTP(NVARMX,2)
      COMMON / IPOSCL / ICOEF , ICOEFF , ICLRTP
C
c@z
