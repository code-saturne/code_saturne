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

subroutine mempdc &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   , nphas  , ndim   , ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE TABLEAU PDC : CKUPDC

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0/idbra0    ! e  ! <-- ! pointeur de la premiere cas libre des          !
!                  !    !     !  tableaux ia/ra                                !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nphas            ! i  ! <-- ! number of phases                               !
! ndim             ! e  ! <-- ! dimension de l'espace (3)                      !
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
integer          ncelet , ncel , nphas, ndim
integer          ifinia , ifinra

integer          idebia, idebra, iok, iok1, iphas

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0


!---> VERIFICATION DES DIMENSIONS

iok1 = 0
do iphas = 1, nphas

  iok = 0
  if(ncepdc.gt.ncelet .or. ncepdc.lt.0) then
    iok = 1
  endif
  if(iok.ne.0) then
    write(nfecra,1000) ncepdc
    iok1 = 1
  endif
enddo

if(iok1.ne.0) then
  call csexit (1)
endif

!---> CALCUL DU NOMBRE DE CELLULES AVEC PDC TOTAL

do iphas = 1, nphas
  ncpdct = ncepdc
enddo
if (irangp.ge.0) then
  call parism(nphas,ncpdct)
endif

!---> QUELQUES MESSAGES

do iphas = 1, nphas
  if(ncpdct.eq.0) then
    write(nfecra,2000) ncpdct
    write(nfecra,3000)
  else
    write(nfecra,2001) ncpdct
    write(nfecra,3000)
  endif
enddo

!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

ifinia = idebia
ifinra = idebra

do iphas = 1, nphas

  iicepd = ifinia
  ifinia        = iicepd + ncepdc

  ickupd = ifinra
  ifinra        = ickupd + ncepdc*6

enddo

!     Si pour une des phases on a des pertes de charge
!       sur un des processeurs
!     et que les tableaux TPUCOU n'ont pas encore ete definis (IPUCOU=0)
!     il faut les dimensionner
if (ipucou.eq.0) then
  iok=0
  do iphas = 1, nphas
    if (ncpdct.gt.0) iok = 1
  enddo
  if (iok.eq.1) then
    itpuco = ifinra
    ifinra = itpuco + ncelet *ndim
  endif
endif



!---> VERIFICATION

call iasize('mempdc',ifinia)
!==========

call rasize('mempdc',ifinra)
!==========

!---> FORMATS

#if defined(_CS_LANG_FR)

 1000 format(/,' SORTIE DANS MEMPDC CAR LES DIMENSIONS DES TABLEAUX ',/,&
         '   RELATIFS AUX PERTES DE CHARGES SONT INCORRECTES ',/, &
         '     NCEPDC = ',I10)

 2000 format(/,                                                   &
 'TRAITEMENT DES PERTES DE CHARGES NON ACTIVE ',/, &
 '                 NCEPDC = ',I10,/)
 2001 format(                                                           &
 /,/,                                                       &
 'TRAITEMENT DES PERTES DE CHARGES ACTIVE ',/,     &
 '                 SUR  UN TOTAL DE NCEPDC = ',I10,' CELLULES',/)

 3000 format(                                                           &
'-------------------------------------------------------------',/)

#else

 1000 format(/,' ABORT IN MEMPDC BECAUSE THE DIMENSION OF THE ARRAYS',/,&
         '   RELATIVE TO THE HEAD LOSS IS INCORRECT ',/,    &
         '     NCEPDC = ',I10)

 2000 format(                                                           &
 /,'HEAD LOSS TREATMENT NOT ACTIVATED ',/,   &
   '                 NCEPDC = ',I10,/)
 2001 format(                                                           &
 /,'HEAD LOSS TERMS TREATMENT ACTIVATED ',/,     &
 '                 ON   A TOTAL OF NCEPDC = ',I10,' CELLS',/)


 3000 format(                                                           &
'-------------------------------------------------------------',/)

#endif

return
end subroutine
