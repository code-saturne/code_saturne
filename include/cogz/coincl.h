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

!                             coincl.h

!===============================================================================

!            INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                           LA COMBUSTION GAZ

! Necessite ppppar.h

!-------------------------------------------------------------------------------

!--> MODELE FLAMME DE DIFFUSION (CHIMIE 3 POINTS)

! ---- Grandeurs fournies par l'utilisateur dans usd3pc.F

!       TINOXY       --> Temperature d'entree pour l'oxydant en K
!       TINFUE       --> Temperature d'entree pour le fuel en K
!       IENTOX       --> indicateur oxydant par type de facette d'entree
!       IENTFU       --> indicateur fuel    par type de facette d'entree

! ---- Grandeurs deduiites

!       HINOXY       --> Enthalpie massique d'entree pour l'oxydant
!       HINFUE       --> Enthalpie massique d'entree pour le fuel
!       HSTOEA       --> Temperature a la stoechiometrie adiabatique
!       NMAXF        --> Nb de points de tabulation en F
!       NMAXFM       --> Nb maximal de points de tabulation en F
!       NMAXH        --> Nb de points de tabulation en H
!       NMAXHM       --> Nb maximal de points de tabulation en H
!       HH           --> Enthalpie stoechiometrique tabulee
!       FF           --> Richesse tabulee
!       TFH(IF,IH)   --> Tabulation richesse - enthalpie stoechiometrique

integer    nmaxf, nmaxfm, nmaxh, nmaxhm
parameter( nmaxfm = 15 , nmaxhm = 15)
integer    ientox(nozppm), ientfu(nozppm)

double precision tinoxy, tinfue, hinfue, hinoxy, hstoea
double precision hh(nmaxhm), ff(nmaxfm), tfh(nmaxfm,nmaxhm)


!--> MODELE FLAMME DE PREMELANGE (MODELE EBU)

! ---- Grandeurs fournies par l'utilisateur dans usebuc.F

!       IENTGF       --> indicateur gaz frais  par type de facette d'entree
!       IENTGB       --> indicateur gaz brules par type de facette d'entree
!       QIMP         --> Debit impose en kg/s
!       FMENT        --> Taux de melange par type de facette d'entree
!       TKENT        --> Temperature en K par type de facette d'entree
!       FRMEL        --> Taux de melange constant pour modeles 0 et 1
!       TGF          --> Temperature gaz frais en K identique
!                        pour premelange frais et dilution
!       CEBU         --> Constante Eddy break-Up

! ---- Grandeurs deduites

!       HGF          --> Enthalpie massique gaz frais identique
!                        pour premelange frais et dilution
!       TGBAD        --> Temperature adiabatique gaz brules en K


integer          ientgf(nozppm), ientgb(nozppm)
double precision fment(nozppm), tkent(nozppm), qimp(nozppm)
double precision frmel, tgf, cebu, hgf, tgbad


!--> DEFINITION DES COMMONS

common / icocom / nmaxf , nmaxh , ientgf, ientgb, ientox, ientfu
common / rcocom / tinoxy, tinfue, hinfue, hinoxy, hstoea,         &
                  hh    , ff    , tfh   ,                         &
                  fment , tkent , qimp  ,                         &
                  frmel , tgf   , cebu  ,                         &
                  hgf   , tgbad


!--> MODELE DE FLAMME DE PREMELANGE LWC

!       NDRACM : nombre de pics de Dirac maximum
!       NDIRAC : nombre de Dirac (en fonction du modele)
  integer ndracm
  parameter (ndracm = 5)

  integer ndirac
  common / ilwcdi / ndirac

! --- Grandeurs fournies par l'utilisateur dans uslwc1.F

!       VREF : Vitesse de reference
!       LREF : Longueur de reference
!         TA : Temperature d'activation
!      TSTAR : Temperature de cross-over

integer irhol(ndracm), iteml(ndracm), ifmel(ndracm)
integer ifmal(ndracm), iampl(ndracm), itscl(ndracm)
integer imaml(ndracm), ihhhh(ndracm), imam

common / rlwcet / irhol, iteml, ifmel, ifmal, iampl,              &
                  itscl, imaml, ihhhh, imam

double precision vref, lref, ta, tstar
double precision fmin, fmax, hmin, hmax
double precision coeff1, coeff2, coeff3

common / rlwcst / vref, lref, ta, tstar,                          &
                  fmin, fmax, hmin, hmax,                         &
                  coeff1, coeff2, coeff3


