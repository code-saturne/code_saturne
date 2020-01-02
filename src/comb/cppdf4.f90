!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine cppdf4 &
!================

 (  ncelet , ncel   ,                                             &
    f1m    , f2m    , f3m    , f4m    , f4p2m  ,                  &
    indpdf ,                                                      &
    si7    , si8    , sp2m   , f4i7   )

!===============================================================================
! FONCTION :
! --------

! CALCUL DU TYPE DE PDF
!  PDF CONJOINTE DEGENERRE EN UNE PDF 1D DE TYPE RECTANGLE - DIRAC
!  SUPPORT (I7,I8) passant par M

!  --> RECONSTITUTION DE 4 MOMENTS : I7 = I4 et I8 = I6
!  -->

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! f1m              ! tr ! <-- ! moyenne du traceur 1 (mv chz +co)              !
! f2m              ! tr ! <-- ! moyenne du traceur 2 (mv cxhy +co)             !
! f3m              ! tr ! <-- ! moyenne du traceur 3 (co c.het)                !
! f4m              ! tr ! <-- ! moyenne du traceur 4 (air)                     !
! f3p2m            ! tr ! <-- ! variance du traceur 3 (co c.het)               !
! indpdf           ! te !  <- ! passage par les pdf                            !
! si7              ! tr !  <- ! abscisse curviligne au point i7                !
! si8              ! tr !  <- ! abscisse curviligne au point i8                !
! sp2m             ! tr !  <- ! variance droite support                        !
! f4i7             ! tr !  <- ! f4 au point i7                                 !
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
use numvar
use optcal
use cstphy
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          indpdf(ncelet)

double precision f1m(ncelet), f2m(ncelet), f3m(ncelet)
double precision f4m(ncelet), f4p2m(ncelet)
double precision si7(ncelet), si8(ncelet), sp2m(ncelet)
double precision f4i7(ncelet)

! Local variables

integer          iel

double precision f1i7, f2i7, f3i7, f4i8, s2max, t1, t2
double precision xf3max

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! --- Parametres relatifs a la variance T1 et a la moyenne T2
t1 = 1.d-04
t2 = 5.d-03

! --- Initialisation des tableaux

do iel = 1, ncel
  f4i7(iel)   = 0.d0
  si7(iel)    = 0.d0
  si8(iel)    = 0.d0
  sp2m(iel)   = 0.d0
  indpdf(iel) = 0
enddo


!===============================================================================
! 2. TYPE DE PDF
!    INDPDF = 0  : PAS DE PDF
!    INDPDF = 3  : SUPPORT (I7,I8)= (I4,M)
!===============================================================================

do iel = 1, ncel
  if ( f4p2m(iel).gt.t1 .and.                                     &
       f4m(iel)  .ge.t2 .and. f4m(iel).le.(1.d0-t2) ) then
    indpdf(iel) = 3
  else
    indpdf(iel) = 0
  endif
enddo


!===============================================================================
! 3. CALCULS DES PARAMETRES DE LA PDF : SI7, SI8 et SP2M
!    CALCUL DE F4I7
!===============================================================================

xf3max =  2.d0*0.012d0 / (2.d0*0.028d0 + xsi*0.028d0)

do iel = 1, ncel

  if ( indpdf(iel).eq.3 ) then
    f4i7(iel) = 1.d0
    f1i7       = 0.d0
    f2i7       = 0.d0
    f3i7       = 0.d0
    si7(iel) = - sqrt (                                           &
      ( sqrt(6.d0)/2.d0*(f1m(iel)-f1i7) +                         &
        sqrt(6.d0)/4.d0*(f2m(iel)+f3m(iel)-f2i7-f3i7) )**2 +      &
      ( 3.d0*sqrt(2.d0)/4.d0*(f2m(iel)-f2i7) +                    &
        sqrt(2.d0)/4.d0*(f3m(iel)-f3i7) )**2 +                    &
       ( f3m(iel)-f3i7 )**2 )
    f4i8 = (1.d0-xf3max)*f3m(iel)                                 &
        / (f3m(iel)+xf3max*(1.d0-f3m(iel)-f4m(iel)))
    si8(iel)  = si7(iel)*(f4m(iel)-f4i8)/(f4m(iel)-f4i7(iel))
    sp2m(iel) = f4p2m(iel)/(f4m(iel)-f4i7(iel))**2*(si7(iel))**2
    s2max     = -si7(iel)*si8(iel)
    if (sp2m(iel).gt.s2max) indpdf(iel) = 0
  endif

enddo



!----
! FIN
!----

return
end subroutine
