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

!                             cpincl.h

!===============================================================================

!            INCLUDE POUR LA PHYSIQUE PARTICULIERE
!                        VARIABLE COMMUNE ENTRE
!                    COMBUSTION DU CHARBON PULVERISE
!                    COMBUSTION DU FIOUL LOURD

!        XSI         --> XSI = 3,76 pour de l'air

double precision  xsi
common / rcpfu1 / xsi

!   nb de moles de I dans J

double precision ao2f3,acof3,an2f3,ah2of3
double precision ao2f4,an2f4,ah2of4,aco2f4
double precision ah2of5
double precision ao2f6,an2f6,ah2of6,aco2f6
double precision ao2f7,an2f7,ah2of7,aco2f7

common / rcpfu2 / ao2f3,acof3,an2f3,ah2of3,                       &
                  ao2f4,an2f4,ah2of4,aco2f4,                      &
                  ah2of5,                                         &
                  ao2f6,an2f6,ah2of6,aco2f6,                      &
                  ao2f7,an2f7,ah2of7,aco2f7

! Equation sur YCO2

 integer         ieqco2 , iyco2
 common/equco2 / ieqco2 , iyco2

! Combustion heterogene avec le  CO2

 integer         ihtco2
 common/ehtco2 / ihtco2

! Equation sur NOX :
! ================

!   IEQNOX = 0 pas de NOx
!          = 1 calcul du NOx

 integer         ieqnox
 common/equnox / ieqnox

!   Scalaires supplementaires : fraction massique de HCN et NO
!                               temperature air

 integer         iyhcn , iyno , itaire
 common/equnox / iyhcn , iyno , itaire

!   Propce supplementaires :

!         Conversion HCN en NO       : EXP(-E1/RT)
!         Conversion HCN en NO       : EXP(-E2/RT)
!         NO thermique (Zel'dovitch) : EXP(-E3/RT)


 integer         ighcn1 , ighcn2 , ignoth
 common/pronox / ighcn1 , ighcn2 , ignoth

!   Temperature moyenne d'entree
!   Taux de vapeur moyen

 double precision taire
 common /noxdbl/  taire


! FIN
