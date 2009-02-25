!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine csinit &
!================

 ( argifo , irgpar , nrgpar , ilisr0, ilisrp )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DU LISTING ET DE PARAMETRES ASSOCIES AU PREPROCESSEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! argifo           ! e  ! <-- ! valeur de ifoenv ;                             !
!                  !    !     ! format de communication avec                   !
!                  !    !     ! le preprocesseur                               !
!                  !    !     !   0 : pas de communications                    !
!                  !    !     !   1 : communication par fichiers               !
! irgpar           ! e  ! <-- ! rang si parallele ; -1 si sequentiel           !
! nrgpar           ! e  ! <-- ! nombre de processus ; 1 si sequentiel          !
! ilisr0           ! e  ! <-- ! option de sortie du listing :                  !
!                  !    !     !   0 : rang 0 non redirige                      !
!                  !    !     !   1 : rang 0 dans fichier listing,             !
! ilisrp           ! e  ! <-- ! option de sortie du listing :                  !
!                  !    !     !   0 : rangs > 0 non rediriges                  !
!                  !    !     !       (pour debugger par exemple)              !
!                  !    !     !   1 : rangs > 0 rediriges dans                 !
!                  !    !     !       fichiers listing_n*                      !
!                  !    !     !   2 : rangs > 0 rediriges dans                 !
!                  !    !     !       /dev/null (suppression)                  !
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
include "optcal.h"
include "entsor.h"
include "parall.h"

!===============================================================================

integer          argifo, irgpar, nrgpar, ilisr0, ilisrp

character        name*300

!===============================================================================

!===============================================================================
! Initialisation du common IPARAL
!===============================================================================

irangp = irgpar
nrangp = nrgpar

!===============================================================================
! Initialisation des paramètres de lecture des données Préprocesseur
!===============================================================================

ifoenv = argifo

!-------------------------------------------
! --- stdout :
!     NFECRA = 6 par défaut
!   Les autres fichiers listing sont fermes dans csclli a la fin de cs_main.
!-------------------------------------------

nfecra = 6

if (irangp.le.0) then
  if (ilisr0.eq.1) then
    nfecra = 9
    NAME = 'listing'
  endif
else
  if (ilisrp.eq.1) then
    nfecra = 9
    WRITE (NAME,'(A9,I4.4)') 'listing_n', IRANGP + 1
  else if (ilisrp.eq.2) then
    nfecra = 9
    NAME = '/dev/null'
  endif
endif

if (nfecra.eq.9) then
   open (file=name, unit=nfecra,                                  &
         FORM='FORMATTED', STATUS='UNKNOWN', ERR=900)
endif

goto 950

 900  write (0, 999) name
call csexit (1)

 950  continue

#if defined(_CS_LANG_FR)

 999  format(/,                                                   &
'Code_Saturne : Erreur d''initialisation :',/,              &
'Impossible d''ouvrir le fichier : ',A,/)

#else

 999  format(/,                                                   &
'Code_Saturne: Initialization error:',/,                    &
'Impossible to open the file: ',A,/)

#endif

return
end
