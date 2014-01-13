!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine cs_fuel_htconvers2 &
!============================
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
! icla             ! e  ! <-- ! numero de la classe                            !
! enthal           ! r  ! <-- ! enthalpie massique j/kg                        !
! xsolid           ! tr ! <-- ! fraction massique des constituants             !
! temper           ! r  ! <-- ! temperature en kelvin                          !
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
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

integer          mode

double precision xsolid(2)
double precision temper, enthal

! Local variables

integer          it, isol, ihflt2

double precision eh1, eh0

!===============================================================================
! Mode de conversion
ihflt2 = 0

if ( ihflt2.eq.0 ) then

!===============================================================================
! 2. H2 FONCTION LINEAIRE T2
!===============================================================================

  if ( mode.eq.-1 ) then

! --> Loi temperature -> enthalpie (MODE = -1)

    enthal =  h02fol + cp2fol*(temper-trefth)

  elseif ( mode.eq.1 ) then

! --> Loi enthalpie -> temperature (MODE = 1)

    temper =  (enthal-h02fol)/cp2fol + trefth
!
    if ( temper .le. th(1)   ) temper = th(1)
    if ( temper .gt. th(npo) ) temper = th(npo)

  else

    write(nfecra,1000) mode
    call csexit (1)
    !==========

  endif


elseif( ihflt2.ne.0 ) then

!===============================================================================
! 3. H2 TABULE
!===============================================================================

  if ( mode.eq.-1 ) then

! --> Loi temperature -> enthalpie (MODE = -1)

    it = npoc
    if ( temper.ge.thc(it) ) then
      enthal = zero
      do isol = 1, nsolid
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif

    it = 1
    if ( temper.le.thc(it) ) then
      enthal = zero
      do isol = 1, nsolid
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif
    it = 1
 10       continue

    it = it + 1
    if ( temper.le.thc(it) ) then
      eh0 = zero
      eh1 = zero
      do isol = 1, nsolid
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it-1)
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it  )
      enddo
      enthal = eh0                                                &
             + (eh1-eh0)*(temper-thc(it-1))/(thc(it)-thc(it-1))
      goto 11
    endif
    goto 10
 11       continue

  elseif ( mode.eq.1 ) then

! --> Loi enthalpie -> temperature (MODE = 1)

    it  = npoc-1
    eh1 = zero
    do isol = 1, nsolid
      eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
    enddo
    if ( enthal.ge.eh1 ) temper = thc(it+1)

    it  = 1
    eh0 = zero
    do isol = 1, nsolid
      eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
    enddo
    if ( enthal.le.eh0 ) temper = thc(it)

    do it = 1, npoc-1
      eh0 = zero
      eh1 = zero
      do isol = 1, nsolid
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
      enddo
      if ( enthal.ge.eh0 .and. enthal.le.eh1 )                      &
        temper = thc(it)+(enthal-eh0)*(thc(it+1)-thc(it))/(eh1-eh0)

    enddo

  else

    write(nfecra,1000) mode
    call csexit (1)
    !==========

  endif

endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS CS_FUEL_HTCONVERS2              ',/,&
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
! End
!----

return
end subroutine
