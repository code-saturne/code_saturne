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

subroutine fupdfr &
!================

 ( ncelet , ncel   ,                                              &
   fvap   , fhtf   , f4p2m  ,                                     &
   indpdf ,                                                       &
   f4m1   , f4m2   , d4cl   , d4f4     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DU TYPE DE PDF

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! fvap             ! tr ! <-- ! moyenne du traceur 1 mvl [fovm+co]             !
! fhet             ! tr !  ->           ! moyenne du traceur 3 c héterogène    !
! f4p2m            ! tr ! <-- ! variance du traceur 4 (air)                    !
! indpdf           ! te ! <-- ! passage par les pdf                            !
! f4m1             ! tr ! <-- ! borne minimum                                  !
! f4m2             ! tr ! <-- ! borne max                                      !
! d4cl             ! tr ! <-- ! amplitude du pic de dirac en f4cl              !
! d4f4             ! tr ! <-- ! amplitude du pic de dirac en 1                 !
!                  !    !     ! (air pur)                                      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!==============================================================================
! Common blocks
!==============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "pointe.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          ncelet , ncel
integer          indpdf(ncelet)

double precision fvap(ncelet) , fhtf(ncelet) , f4p2m(ncelet)
double precision f4m1(ncelet) , f4m2(ncelet) , d4cl(ncelet)
double precision d4f4(ncelet)

! Local variables

integer          iel

!   Variables pour le support de la Pdf et sa description

double precision t1,t2
double precision f4m,f42m,f1cl,f3cl,f4cl,f1m,f3m

double precision d2pd1,d2md1

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

!     Bornes pour le passage par les pdf T1 sur var, T2 sur moy

t1 = 1.d-4
t2 = 5.d-3

! --- Initialisation des tableaux

do iel = 1, ncel
  f4m1(iel)   = 0.d0
  f4m2(iel)   = 0.d0
  d4cl(iel)   = 0.d0
  d4f4(iel)  = 0.d0

  indpdf(iel) = 0
enddo

!===============================================================================
! 2. TYPE DE PDF
!    INDPDF = 0  : PAS DE PDF
!    INDPDF = 1  : Passage par les PDF

!    CALCULS DES PARAMETRES DE LA PDF : F4M1,F4M2,D4CL,D4F4
!    CALCUL DE F4I7


!===============================================================================

do iel = 1, ncel

!       Traceur virtuel au point moyen
  f1m =  fvap(iel)
  f3m =  fhtf(iel)/ff3max
  f4m = 1.d0 - f1m - f3m

!       Calcul des caractéristiques du point correspondant au
!       combustible moyen
!       F3cl : fraction de masse provenant de F3max et non de F3

  f1cl = f1m*ff3max/(f3m+f1m*ff3max)
  f3cl = 1.d0-f1cl
  f4cl = (1.d0-ff3max)*f3cl

!       L'intervalle permis pour F4 est [F4cl , 1 ]
!       La valeur maximale de la variance est donc
!       F4p2max = (F4m-F4cl)*(1-F4m)
!       On ne passe par les PDF que si la variance représente
!       une fraction significative de son maximum

  if (  f4m .gt. ( f4cl + t2 ) .and.                              &
       f4m .lt. ( 1.d0 - t2 ) .and.                               &
       f4p2m(iel) .gt. (t1*(f4m-f4cl)*(1.d0-f4m)) ) then
    indpdf(iel) = 1
  else
    indpdf(iel) = 0
  endif

  if ( indpdf(iel) .eq.1  )then

!        Calcul préliminaire pour une pdf rectangle et pics de Dirac
!        sur une droite passant par l'entrée d'air (F4) et le
!        combsutible moyen local Fcl
!        Soit F4m1 la borne Min
!             F4m2 la borne Max
!             D4cl amplitude du pic de Dirac en F4cl
!             D4f4 amplitude du pic de Dirac en 1 (air pur)

    f42m = f4m**2 + f4p2m(iel)

!         rectangle seul

    f4m1(iel) = f4m - sqrt(3.d0*f4p2m(iel))
    f4m2(iel) = f4m + sqrt(3.d0*f4p2m(iel))
    d4cl(iel) = zero
    d4f4(iel) = zero

    if ( f4m1(iel).le.f4cl .or. f4m2(iel).ge. 1.d0) then

!         pic en Fcl

      f4m1(iel) = f4cl
      f4m2(iel) = (3.d0*f42m + f4cl**2 - 4.d0*f4cl*f4m)           &
                 /(2.d0*(f4m-f4cl))
      d4cl(iel) = (f4cl+f4m2(iel)-2.d0*f4m)                       &
                 /(f4m2(iel)-f4cl)

      if (f4m2(iel).ge.1.d0 .or. d4cl(iel).le.zero) then

!           pic en F4

        f4m2(iel) = 1.d0
        f4m1(iel) = (3.d0*f42m + 1.d0 - 4.d0*f4m)                 &
                   /(2.d0*(f4m-1.d0))
        d4f4(iel) = (2.d0*f4m-1.d0-f4m1(iel))                     &
                   /(1.d0-f4m1(iel))
        d4cl(iel) = 0.d0

        if (f4m1(iel).le.f4cl .or. d4f4(iel).le.zero) then

          f4m1(iel) = f4cl
          f4m2(iel) = 1.d0
          d2pd1 = (2.d0*f4m-1.d0-f4cl)/(1.d0-f4cl)
          d2md1 = ( 6.d0*(f42m-f4m*(1.d0+f4cl))                   &
                  + 1.d0+f4cl**2+4.d0*f4cl)                       &
                 / (1.d0-f4cl)**2
          d4f4(iel) = 0.5d0*(d2pd1+d2md1)
          d4cl(iel) = 0.5d0*(d2pd1-d2md1)
        endif

      endif

    endif

  endif

enddo

!----
! FIN
!----

return
end subroutine
