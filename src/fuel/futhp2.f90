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

subroutine futhp2 &
!================

 ( mode   , enthal , xsolid , temper )

!===============================================================================
!  FONCTION  :
!  --------
! CALCUL DE LA TEMPERATURE DES PARTICULES
!  EN FONCTION DE L'ENTHALPIE ET DES CONCENTRATIONS
!  SI IMODE = 1
! CALCUL DE L'ENTHALPIE DES PARTICULES
!  EN FONCTION DE LA TEMPERATURE ET DES CONCENTRATIONS
!  SI IMODE = -1

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mode             ! e  ! <-- !  -1 : t -> h  ;   1 : h -> t                   !
! enthal           ! r  ! <-- ! enthalpie massique j/kg                        !
! xsolid           ! tr ! <-- ! fraction massique des constituants             !
! temper           ! r  ! <-- ! temperature en kelvin                          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "pointe.f90"
include "entsor.f90"
include "cstnum.f90"
include "cstphy.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "fuincl.f90"
include "ppincl.f90"

!===============================================================================


! Arguments

integer          mode

double precision xsolid(2)
double precision temper , enthal

! Local variables

integer          it , isol , ihflt2

double precision eh1 , eh0

!===============================================================================
!===============================================================================
! 1. RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)
!===============================================================================

ihflt2 = 0

if ( ihflt2.eq.0 ) then

!===============================================================================
! 2. H2 FONCTION LINEAIRE T2
!===============================================================================


  if ( mode.eq.-1 ) then

! --> Loi temperature -> enthalpie (MODE = -1)

    enthal = h02fol + cp2fol*(temper-trefth)

  elseif ( mode.eq.1 ) then

! --> Loi enthalpie -> temperature (MODE = 1)

    temper =  (enthal-h02fol)/cp2fol + trefth

  else

    write(nfecra,1000) mode
    call csexit (1)
    !==========

  endif


else

!===============================================================================
! 3. H2 TABULE
!===============================================================================

  if ( mode.eq.-1 ) then

! --> Loi temperature -> enthalpie (MODE = -1)

    it = npo
    if ( temper.ge.th(it) ) then
      enthal = zero
      do isol = 1, 2
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif

    it = 1
    if ( temper.le.th(it) ) then
      enthal = zero
      do isol = 1, 2
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif
    it = 1
 10       continue

    it = it + 1
    if ( temper.le.th(it) ) then
      eh0 = zero
      eh1 = zero
      do isol = 1, 2
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it-1)
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it  )
      enddo
      enthal = eh0                                                &
             + (eh1-eh0)*(temper-th(it-1))/(th(it)-th(it-1))
      goto 11
    endif
    goto 10
 11       continue

  elseif ( mode.eq.1 ) then

! --> Loi enthalpie -> temperature (MODE = 1)

    it  = npo-1
    eh1 = zero
    do isol = 1, 2
      eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
    enddo
    if ( enthal.ge.eh1 ) temper = th(it+1)

    it  = 1
    eh0 = zero
    do isol = 1, 2
      eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
    enddo
    if ( enthal.le.eh0 ) temper = th(it)

    do it = 1, npo-1
      eh0 = zero
      eh1 = zero
      do isol = 1, 2
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
      enddo
      if ( enthal.ge.eh0 .and. enthal.le.eh1 )                    &
        temper = th(it)                                           &
               + (enthal-eh0)*(th(it+1)-th(it))/(eh1-eh0)
    enddo

  else

    write(nfecra,1000) mode
    call csexit (1)
    !==========

  endif

endif


!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS FUTHP2                          ',/,&
'@    =========                                               ',/,&
'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,&
'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return
end subroutine
