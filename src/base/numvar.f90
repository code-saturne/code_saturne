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

!                              numvar.h
!===============================================================================

! POSITION DES VARIABLES
!  ( Dans RTP, RTPA )

! IPR                        Pression
! IU   IV   IW               Vitesse(x,y,z)
! IK                         Energie Turbulente en k-epsilon
! IR11, IR22, IR33, ...
! ... IR12, IR13, IR23       Tensions de Reynolds en Rij
! IEP                        Dissipation turbulente
! IPHI, IFB                  Variables phi et f_barre du v2f phi-model
! IOMG                       Variable Omega du k-omega SST
! ISCA(i)                    Scalaire numero i
! ISCAPP(i)                  No du scalaire physique particuliere i
! NSCAUS                     Nbre de scalaires utilisateur
! NSCAPP                     Nbre de scalaires physique particuliere
! IUMA, IVMA, IWMA           Vitesse de maillage en ALE


integer           ipr(nphsmx) ,                                   &
                  iu(nphsmx)  , iv(nphsmx)    , iw(nphsmx)  ,     &
                  ik(nphsmx)  , iep(nphsmx)   ,                   &
                  ir11(nphsmx), ir22(nphsmx)  , ir33(nphsmx),     &
                  ir12(nphsmx), ir13(nphsmx)  , ir23(nphsmx),     &
                  iphi(nphsmx), ifb(nphsmx)   , iomg(nphsmx),     &
                  isca(nscamx), iscapp(nscamx),                   &
                  nscaus      , nscapp        ,                   &
                  iuma        , ivma          , iwma
common / iposvr / ipr         ,                                   &
                  iu          , iv            , iw          ,     &
                  ik          , iep           ,                   &
                  ir11        , ir22          , ir33        ,     &
                  ir12        , ir13          , ir23        ,     &
                  iphi        , ifb           , iomg        ,     &
                  isca        , iscapp        ,                   &
                  nscaus      , nscapp        ,                   &
                  iuma        , ivma          , iwma

! POSITION DES PROPRIETES (Physiques ou numeriques)
!  ( Dans PROPCE, PROPFA et PROPFB)
!    Le numero des proprietes est unique, quelle aue soit la
!      localisation de ces dernieres (cellule, face, face de bord)
!    Voir usclim pour quelques exemples

! IPPROC : Pointeurs dans PROPCE
! IPPROF : Pointeurs dans PROPFA
! IPPROB : Pointeurs dabs PROPFB

! IROM   : Masse volumique des phases
! IROMA  : Masse volumique des phases au pas de temps precedent
! IVISCL : Viscosite moleculaire dynamique en kg/(ms) des phases
! IVISCT : Viscosite turbulente des phases
! IVISLA : Viscosite moleculaire dynamique en kg/(ms) des phases au pas
!          de temps precedent
! IVISTA : Viscosite turbulente des phases au pas de temps precedent
! ICP    : Chaleur specifique des phases
! ICPA   : Chaleur specifique des phases au pas de temps precedent
! ITSNSA : Terme source Navier Stokes des phases au pas de temps precedent
! ITSTUA : Terme source des grandeurs turbulentes au pas de temps precedent
! ITSSCA : Terme source des scalaires au pas de temps precedent
! IESTIM : Estimateur d'erreur pour Navier-Stokes
! IFLUMA : Flux de masse associe aux variables
! IFLUAA : Flux de masse explicite (plus vu comme un tableau de travail)
!          associe aux variables
! ISMAGO : constante de Smagorinsky dynamique
! ICOUR  : Nombre de Courant des phases
! IFOUR  : Nombre de Fourier des phases
! IPRTOT : Pression totale au centre des cellules Ptot=P*+rho*g.(x-x0)
!                                                             -  - -
! IVISMA : Viscosite de maillage en ALE (eventuellement orthotrope)

integer           ipproc(npromx), ipprof(npromx), ipprob(npromx), &
                  irom  (nphsmx), iroma (nphsmx), iviscl(nphsmx), &
                  ivisct(nphsmx), ivisla(nphsmx), ivista(nphsmx), &
                  icp   (nphsmx), icpa  (nphsmx), itsnsa(nphsmx), &
                  itstua(nphsmx), itssca(nscamx),                 &
                  iestim(nestmx,nphsmx)         , ifluma(nvarmx), &
                  ifluaa(nvarmx), ismago(nphsmx), icour (nphsmx), &
                  ifour (nphsmx), iprtot(nphsmx), ivisma(3)
common / ipospp / ipproc        , ipprof        , ipprob        , &
                  irom          , iroma         , iviscl        , &
                  ivisct        , ivisla        , ivista        , &
                  icp           , icpa          , itsnsa        , &
                  itstua        , itssca        ,                 &
                  iestim                        , ifluma        , &
                  ifluaa        , ismago        , icour         , &
                  ifour         , iprtot        , ivisma


! POSITION DES CONDITIONS AUX LIMITES
!  (Position dans COEFA et COEFB des coef (coef. coef.f) relatifs a
!   une variable donnee)

! ICOEF   : coef numeros 1 (anciens coefa coefb)
! ICOEFF  : coef numeros 2 (anciens coefaf coefbf)
! ICLRTP  : pointeur dans COEFA et COEFB

integer           icoef , icoeff , iclrtp(nvarmx,2)
common / iposcl / icoef , icoeff , iclrtp

