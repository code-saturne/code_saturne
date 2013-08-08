!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!> \file ppthch.f90
!> Module for specific physics thermophysical data

module ppthch

  !===========================================================================

  use cstphy

  implicit none

  !===========================================================================

  !--> CONSTANTES THERMOCHIMIE

  !       RR           --> Constante des gaz parfaits en J/mol/K
  !       TREFTH       --> Temperature de reference (K)
  !       VOLMOL       --> Volume molaire dans les conditions NTP
  !                        T = 0 C et P = 1 atm

  double precision rr
  double precision trefth, prefth, volmol
  parameter ( rr     = 8.31434d0      ,                             &
              trefth = 25.d0 + tkelvi ,                             &
              prefth = 1.01325d5      ,                             &
              volmol = 22.41d-3       )

  !--> DONNEES

  !       NRGAZ        --> Nb de reactions globales en phase gaz
  !       NRGAZM       --> Nb maximal de reactions globales en phase gaz
  !       NATO         --> Nb d especes atomiques (C,H,..)
  !       NATOM        --> Nb maximal d especes atomiques (C,H,..)
  !       NGAZE        --> Nb de constituants gazeux elementaires
  !       NGAZEM       --> Nb maximal de constituants gazeux elementaires
  !       NGAZG        --> Nb d especes globales (ex:Fuel,Oxyd,Prod1,Prod2)
  !       NGAZGM       --> Nb maximal d especes globales
  !       NPO          --> Nb de points de tabulation
  !       NPOT         --> Nb maximal de points de tabulation
  !       TH           --> Temperature en Kelvin
  !       EHGAZG(G,IT) --> Enthalpie massique (J/kg) de l espece globale
  !                        no G a la temperature T(IT)
  !       WMOLG(G)     --> Masse molaire de l espece globale
  !       EHGAZE(G)    --> Enthalpie massique (J/kg) constituant gazeux
  !                        elementaire no E a la temperature T(IT)
  !       WMOLE(G)     --> Masse molaire du constituant gazeux elementaire
  !       WMOLAT(E)    --> Masse molaire des atomes (C,H,..)
  !       IATC, IATH   --> Pointeur dans WMOLEL pour les ecpeces
  !       IATO, IATN, IATS       elementaires (C,H,..)
  !       FS(R)        --> Taux de melange pour la reaction gloable R
  !       STOEG(G,R)   --> Stoechio en especes globales des reactions
  !                        pour l espece no G et pour la reaction no R
  !       CKABSG(G)    --> Coefficient d'absorption des especes globales
  !       CKABS1       --> Coefficient d'absorption du melange gazeux
  !                        (en CP)
  !       DIFTL0       --> Diffusivite dynamique en kg/(m s)

  integer    ngazgm, ngazem, npot, natom, nrgazm
  parameter( ngazgm = 25 , ngazem = 20 ,                                     &
             npot  = 500 , natom  = 5   , nrgazm = 1 )
  integer    iatc, iath, iato, iatn , iats
  parameter( iatc = 1, iath = 2, iato = 3, iatn = 4 , iats = 5 )

  integer, save ::           npo, ngaze, ngazg, nato, nrgaz

  double precision, save ::  th(npot),                                       &
                             ehgaze(ngazem,npot), ehgazg(ngazgm,npot),       &
                             wmole(ngazem), wmolg(ngazgm), wmolat(natom),    &
                             stoeg(ngazgm,nrgazm), fs(nrgazm),               &
                             ckabsg(ngazgm), ckabs1,                         &
                             diftl0, xco2, xh2o

  !=============================================================================

end module ppthch


