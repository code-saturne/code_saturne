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

subroutine elthht &
!================

 ( mode   , nesp , yesp , enthal , temper )

!===============================================================================
!  FONCTION  :
!  --------

! CALCULE L'ENTHALPIE OU LA TEMPERATURE  A PARTIR DE LA
!  COMPOSITION ET DE LA VALEUR DE LA TEMPERATURE OU DE
!                    ENTHALPIE
!        SPECIFIQUE AU MODULE ELECTRIQUE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! mode             ! e  ! <-- !  -1 : t -> h  ;   1 : h -> t                   !
! nesp             ! e  ! <-- ! nb de constituants                             !
! yesp             ! tr ! <-- ! fraction massique des constituants             !
! enthal           ! r  ! <-- ! enthalpie                                      !
! temper           ! r  ! <-- ! temperature                                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppthch.h"
include "elincl.h"

!===============================================================================


! Arguments

integer          mode , nesp

double precision enthal , temper , yesp(nesp)


! VARIABLES LOCALES

integer          it , iesp

double precision eh1 , eh0


!===============================================================================

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

if ( mode.eq.-1 ) then

  it=npo
  if (temper.ge.th(it)) then
    enthal = 0
    do iesp=1,nesp
      enthal = enthal + yesp(iesp)*ehgazg(iesp,it)
    enddo
    return
  endif

  it=1
  if (temper.le.th(it)) then
    enthal = 0
    do iesp=1,nesp
      enthal = enthal + yesp(iesp)*ehgazg(iesp,it)
    enddo
    return
  endif

 10     continue
    it=it+1
    if (temper.le.th(it)) then
      eh0 = 0
      eh1 = 0
      do iesp=1,nesp
        eh0 = eh0 + yesp(iesp)*ehgazg(iesp,it-1)
        eh1 = eh1 + yesp(iesp)*ehgazg(iesp,it)
      enddo
      enthal=eh0+(eh1-eh0)*(temper-th(it-1))                      &
                          /(th(it)-th(it-1))
      return
    endif
  goto 10


!===============================================================================
! 2. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
!===============================================================================

else if ( mode.eq.1 ) then

  it  = npo
  eh1 = 0
  do iesp=1,nesp
    eh1 = eh1+yesp(iesp)*ehgazg(iesp,it)
  enddo
  if ( enthal .ge. eh1 ) then
    temper = th (it)
    return
  endif

  it  = 1
  eh1 = 0
  do iesp=1,nesp
    eh1 = eh1+yesp(iesp)*ehgazg(iesp,it)
  enddo
  if ( enthal .le. eh1 ) then
     temper = th (it)
     return
  endif

 20     continue
    it  = it+1
    eh0 = eh1
    eh1 = 0
    do iesp=1,nesp
      eh1 = eh1+yesp(iesp)*ehgazg(iesp,it)
    enddo
    if ( enthal .le. eh1 ) then
      temper = th(it-1)+ (enthal-eh0)                             &
                        *(th(it)-th(it-1))/(eh1-eh0)
      return
    endif
  goto 20


else

  write(nfecra,1000) mode
  call csexit (1)
  !==========


endif


!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS ELTHHT                          ',/,&
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
