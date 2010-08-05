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

!                              alstru.h
!===============================================================================

!  METHODE ALE - MOUVEMENT DE STRUCTURES EN COUPLAGE INTERNE

! NBSTRU : NOMBRE DE STRUCTURES MOBILES

! XMSTRU : MATRICE DE MASSE DES STRUCTURES (kg)
! XCSTRU : MATRICE DE FRICTION DES STRUCTURES (kg/s)
! XKSTRU : MATRICE DE RAIDEUR DES STRUCTURES (kg/s2 = N/m)
! XSTREQ : VECTEUR ECART DE LA POSITION DES STRUCTURE DANS LE MAILLAGE
!          INITIAL PAR RAPPORT A LEUR POSITION D'EQUILIBRE (m)
! XSTR   : VECTEUR DEPLACEMENT DES STRUCTURES PAR RAPPORT A LEUR POSITION
!          DANS LE MAILLAGE INITIAL (m)
! XPSTR  : VECTEUR VITESSE DES STRUCTURES (m/s)
! XPPSTR : VECTEUR ACCELERATION DES STRUCTURES (m/s2)
! FORSTR : VECTEUR FORCE EXERCE SUR LES STRUCTURES (N)
! XSTA   : VALEUR DE XSTR AU PAS DE TEMPS PRECEDENT
! XPSTA  : VALEUR DE XPSTR AU PAS DE TEMPS PRECEDENT
! XPPSTA : VALEUR DE XPPSTR AU PAS DE TEMPS PRECEDENT
! FORSTA : VALEUR DE FORSTR AU PAS DE TEMPS PRECEDENT
! XSTP   : VALEUR PREDITE DE XSTR
! FORSTP : VALEUR PREDITE DE FORSTR
! DTSTR  : PAS DE TEMPS ASSOCIE AU MOUVEMENT DES STRUCTURES
! AEXXST : COEFFICIENT DE PREDICTION DU DEPLACEMENT (SUR XPSTR)
! BEXXST : COEFFICIENT DE PREDICTION DU DEPLACEMENT (SUR XPSTR-XPSTA)
! CFOPRE : COEFFICIENT DE PREDICTION DES EFFORTS


integer          nbstru

common / istruc / nbstru

double precision xmstru(3,3,nstrmx)
double precision xcstru(3,3,nstrmx)
double precision xkstru(3,3,nstrmx)
double precision xstr(3,nstrmx)  ,xsta(3,nstrmx)
double precision xstp(3,nstrmx)  ,xstreq(3,nstrmx)
double precision xpstr(3,nstrmx) ,xpsta(3,nstrmx)
double precision xppstr(3,nstrmx),xppsta(3,nstrmx)
double precision forstr(3,nstrmx),forsta(3,nstrmx)
double precision forstp(3,nstrmx)
double precision dtstr(nstrmx)
double precision aexxst, bexxst, cfopre

common / rstruc / xmstru, xcstru, xkstru, xstr  , xsta  , xstp,   &
                  xstreq, xpstr , xpsta , xppstr, xppsta,         &
                  forstr, forsta, forstp, dtstr,                  &
                  aexxst, bexxst, cfopre

!  METHODE DE NEWMARK HHT

!  ALPNMK : COEFFICIENT ALPHA DE LA METHODE DE NEWMARK HHT
!  BETNMK : COEFFICIENT BETA  DE LA METHODE DE NEWMARK HHT
!  GAMNMK : COEFFICIENT GAMMA DE LA METHODE DE NEWMARK HHT
double precision alpnmk, gamnmk, betnmk
common / rnemrk / alpnmk, gamnmk, betnmk

! FIN

