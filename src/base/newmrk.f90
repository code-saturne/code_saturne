!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine newmrk &
!================

 (istr  , alpnmk, betnmk, gamnmk,                                 &
  xm    , xc    , xk    , xn0   , xn    , xpn   , xppn  , xnm1  , &
  xpnm1 , xppnm1, xfn   , xfnm1 , dt    )


!===============================================================================
! FONCTION :
! ----------

! RESOLUTION PAR LA METHODE DE NEWMARK HHT D'UNE EQUATION
!   DIFFERENTIELLE LINEAIRE DE SECOND ORDRE DU TYPE

!    M.X'' + C.X' + K.(X+X0) = F

! X EST UN CHAMP DE VECTEUR TRIDIMENSIONNEL ET M, C ET K
!   SONT DES MATRICES 3x3 QUELCONQUES

! LA RESOLUTION EST FAITE EN DEPLACEMENT

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! istr             ! e  ! <-- ! numero de structure                            !
! xm(3,3)          ! tr ! <-- ! matrice de masse du systeme                    !
! xc(3,3)          ! tr ! <-- ! matrice de friction du systeme                 !
! xk(3,3)          ! tr ! <-- ! matrice de raideur du systeme                  !
! xn0(3)           ! tr ! <-- ! champ de deplacement initial                   !
! xn(3)            ! tr ! --> ! champ de deplacement au temps n                !
! xnm1(3)          ! tr ! <-- ! champ de deplacement au temps n-1              !
! xpn(3)           ! tr ! --> ! champ de vitesse au temps n                    !
! xpnm1(3)         ! tr ! <-- ! champ de vitesse au temps n-1                  !
! xppn(3)          ! tr ! --> ! champ d'acceleration au temps n                !
! xppnm1(3)        ! tr ! <-- ! champ d'acceleration au temps n-1              !
! xftld(3)         ! tr ! <-- ! champ de force au temps intermediaire          !
! dt               ! r  ! <-- ! pas de temps                                   !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor

!===============================================================================

implicit none

! Arguments

integer          istr

double precision alpnmk, betnmk, gamnmk
double precision xm(3,3),xc(3,3),xk(3,3)
double precision xn0(3),xn(3),xpn(3),xppn(3)
double precision xnm1(3),xpnm1(3),xppnm1(3)
double precision xfn(3),xfnm1(3)
double precision dt

! Local variables

integer          ii, jj

double precision a(3,3), b1(3), b2(3), b(3)
double precision a0, a1, a2, a3, a4, a5, a6, a7
double precision det, det1, det2, det3, epsdet

!===============================================================================

!     Critere de nullite du determinant
epsdet = 1.d-12

!     Ceofficients des equations
a0 = 1.d0/betnmk/dt**2
a1 = (1.d0+alpnmk)*gamnmk/betnmk/dt
a2 = 1.d0/betnmk/dt
a3 = 1.d0/2.d0/betnmk - 1.d0
a4 = (1.d0+alpnmk)*gamnmk/betnmk - 1.d0
a5 = (1.d0+alpnmk)*dt*(gamnmk/2.d0/betnmk - 1.d0)
a6 = dt*(1.d0 - gamnmk)
a7 = gamnmk*dt

do ii = 1, 3
  do jj = 1, 3
    a(ii,jj) = (1.d0+alpnmk)*xk(ii,jj)                            &
         + a1*xc(ii,jj) + a0*xm(ii,jj)
  enddo
  b(ii)  = (1.d0+alpnmk)*xfn(ii) - alpnmk*xfnm1(ii)
  b1(ii) = a0*xnm1(ii) + a2*xpnm1(ii) + a3*xppnm1(ii)
  b2(ii) = a1*xnm1(ii) + a4*xpnm1(ii) + a5*xppnm1(ii)
enddo

do ii = 1, 3
  do jj = 1, 3
    b(ii) = b(ii)                                                 &
         + xm(ii,jj)*b1(jj) + xc(ii,jj)*b2(jj)                    &
         + xk(ii,jj)*( alpnmk*xnm1(jj) + xn0(jj) )
  enddo
enddo

det = a(1,1)*a(2,2)*a(3,3)                                        &
    + a(2,1)*a(3,2)*a(1,3)                                        &
    + a(3,1)*a(1,2)*a(2,3)                                        &
    - a(3,1)*a(2,2)*a(1,3)                                        &
    - a(2,1)*a(1,2)*a(3,3)                                        &
    - a(1,1)*a(3,2)*a(2,3)

if (dabs(det).lt.epsdet) then
  write(nfecra,1000) istr,dabs(det),epsdet
  ntmabs = ntcabs
endif

det1 = b(1)*a(2,2)*a(3,3)                                         &
     + b(2)*a(3,2)*a(1,3)                                         &
     + b(3)*a(1,2)*a(2,3)                                         &
     - b(3)*a(2,2)*a(1,3)                                         &
     - b(2)*a(1,2)*a(3,3)                                         &
     - b(1)*a(3,2)*a(2,3)

det2 = a(1,1)*b(2)*a(3,3)                                         &
     + a(2,1)*b(3)*a(1,3)                                         &
     + a(3,1)*b(1)*a(2,3)                                         &
     - a(3,1)*b(2)*a(1,3)                                         &
     - a(2,1)*b(1)*a(3,3)                                         &
     - a(1,1)*b(3)*a(2,3)

det3 = a(1,1)*a(2,2)*b(3)                                         &
     + a(2,1)*a(3,2)*b(1)                                         &
     + a(3,1)*a(1,2)*b(2)                                         &
     - a(3,1)*a(2,2)*b(1)                                         &
     - a(2,1)*a(1,2)*b(3)                                         &
     - a(1,1)*a(3,2)*b(2)

xn(1) = det1/det
xn(2) = det2/det
xn(3) = det3/det

do ii = 1, 3
  xppn(ii) = a0*(xn(ii)-xnm1(ii)) - a2*xpnm1(ii)  - a3*xppnm1(ii)
  xpn(ii)  = xpnm1(ii)            + a6*xppnm1(ii) + a7*xppn(ii)
enddo

!----
! Formats
!----

#if defined(_CS_LANG_FR)

 1000 format (                                                    &
'@                                                            ',/,&
'@ @@ ATTENTION : DEPLACEMENT DE STRUCTURES INTERNES ALE      ',/,&
'@    =========                                               ',/,&
'@  Structure : ',I10                                          ,/,&
'@  La valeur absolue du determinant de la matrice de         ',/,&
'@    deplacement vaut : ',E14.5                               ,/,&
'@  La matrice est consideree comme non inversible            ',/,&
'@    (valeur limite fixee a ',E14.5     ,')                  ',/,&
'@                                                            ',/,&
'@  Arret du calcul                                           ',/,&
'@                                                            '  )

#else

 1000 format (                                                    &
'@                                                            ',/,&
'@ @@ WARNING: ALE DISPLACEMENT OF INTERNAL STRUCTURES        ',/,&
'@    ========                                                ',/,&
'@  Structure: ',I10                                           ,/,&
'@  The absolute value of the discriminant of the             ',/,&
'@    displacement matrix is: ',E14.5                          ,/,&
'@  The matrix is considered to be not inversible             ',/,&
'@    (limit value fixed to ',E14.5     ,')                   ',/,&
'@                                                            ',/,&
'@  Calculation abort                                         ',/,&
'@                                                            '  )

#endif

!----
! End
!----

end subroutine
