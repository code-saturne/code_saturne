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

subroutine memtsm &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   , nvar   , nphas  ,                            &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE POUR LES TERMES SOURCES DE MASSE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! nvar             ! e  ! <-- ! nombre de variables                            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nphas            ! i  ! <-- ! number of phases                               !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
!__________________.____._____.________________________________________________.

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
use numvar
use entsor
use pointe
use parall

!===============================================================================

implicit none

integer          idbia0 ,idbra0
integer          ncelet , ncel  ,nvar   , nphas
integer          ifinia , ifinra

integer          idebia, idebra, iok1, iphas
integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

!---> VERIFICATION DES DIMENSIONS

iok1 = 0
do iphas = 1, nphas
  if(ncetsm(iphas).gt.ncelet .or. ncetsm(iphas).lt.0) then
    write(nfecra,1000) iphas, ncetsm(iphas)
    iok1 = 1
  endif
enddo
if(iok1.ne.0) then
  call csexit (1)
endif

!---> CALCUL DU NOMBRE DE CELLULES AVEC TSM TOTAL

do iphas = 1, nphas
  nctsmt(iphas) = ncetsm(iphas)
enddo
if (irangp.ge.0) then
  call parism(nphas,nctsmt)
endif

!---> QUELQUES MESSAGES

do iphas = 1, nphas
  if(nctsmt(iphas).eq.0) then
    write(nfecra,2000) iphas, nctsmt(iphas)
    write(nfecra,3000)
  else
    write(nfecra,2001) iphas, nctsmt(iphas)
    write(nfecra,3000)
  endif
enddo

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia = idebia
ifinra = idebra

do iphas = 1, nphas

  iicesm(iphas) = ifinia
  iitpsm(iphas) = iicesm(iphas) + ncetsm(iphas)
  ifinia        = iitpsm(iphas) + ncetsm(iphas)*nvar

  ismace(iphas) = ifinra
  ifinra        = ismace(iphas) + ncetsm(iphas)*nvar

enddo

!---> VERIFICATION

call iasize('memtsm',ifinia)
!==========

call rasize('memtsm',ifinra)
!==========

!---> FORMATS

#if defined(_CS_LANG_FR)

 1000 format(/,' SORTIE DANS MEMTSM CAR LA DIMENSIONS DU TABLEAU ',/,   &
         '   RELATIF AUX SOURCES DE MASSE EST INCORRECTE ',/,     &
         '   PHASE ',I10,/,                                 &
         '   NCETSM = ',I10)

 2000 format(                                                           &
 /,'PHASE ',I6,' : TRAITEMENT DES SOURCES DE MASSE NON ACTIVE ',/,&
   '                 NCETSM = ',I10,/)
 2001 format(                                                           &
 /,'PHASE ',I6,' : TRAITEMENT DES SOURCES DE MASSE ACTIVE ',/,    &
   '                 SUR  UN TOTAL DE ',I10,' CELLULES')

 3000 format(                                                           &
'-------------------------------------------------------------',/)

#else

 1000 format(/,' ABORT IN MEMTSM BECAUSE THE DIMENSION OF THE ARRAY ',/,&
         '   RELATIVE TO THE MASS SOURCE TERMS IS INCORRECT ',/,  &
         '   PHASE ',I10,/,                                 &
         '   NCETSM = ',I10)

 2000 format(                                                           &
 /,'PHASE ',I6,' : MASS SOURCE TERMS TREATMENT NOT ACTIVATED ',/, &
   '                 NCETSM = ',I10,/)
 2001 format(                                                           &
 /,'PHASE ',I6,' : MASS SOURCE TERMS TREATMENT ACTIVATED ',/,     &
   '                 ON A TOTAL OF ',I10,' CELLS')

 3000 format(                                                           &
'-------------------------------------------------------------',/)

#endif

return
end subroutine
