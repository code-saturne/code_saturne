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

subroutine modpar &
!================

 ( ntcabs , ntmabs )

!===============================================================================

! FONCTION :
! --------
!          MODIFICATION DE NTMABS AU COURS DU CALCUL
!          POUR ARRET INTERACTIF

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ntcabs           ! e  ! <-- ! numero absolu du pas de temps courant          !
! ntmabs           ! e  ! <-- ! numero absolu du pas de temps final            !
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
use entsor
use parall

!===============================================================================

implicit none

! Arguments

integer          ntcabs , ntmabs

! Local variables

integer          irangs, lng, itmp(1)
logical exstp

!===============================================================================

! On passe ici en sequentiel, ou en parallele avec le proc 0.
if (irangp.le.0) then
!---> ARRET D'URGENCE

  inquire (file=ficstp,exist=exstp)

!       SI UN FICHIER FICSTP EXISTE

  if (exstp) then

!         LIRE LE NOMBRE MAX D'ITERATIONS (ABSOLU)

    open (file=ficstp,unit=impstp)
    read (impstp,*,err=5200,end=5200)
 5200     read (impstp,*,err=5100,end=5100)ntmabs
 5100     continue
    CLOSE (IMPSTP,STATUS='DELETE')

!         COMPARER LE TEMPS ECOULE ET LE TEMPS MAX
!          MODIFIER FICSTP SI BESOIN

    if(ntcabs.gt.ntmabs)then
      ntmabs = ntcabs
    endif

!         SORTIES

    write (nfecra,1000) ntcabs,ntmabs

    OPEN (FILE=FICSTP//'.mod',UNIT=IMPSTP)
    write (impstp,1000) ntcabs,ntmabs
    close (impstp)
  endif

endif

! En parallele, bcast.
if(irangp.ge.0) then
  irangs  = 0
  lng     = 1
  itmp(1) = ntmabs
  call parbci(irangs,lng,itmp)
  ntmabs = itmp(1)
endif

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'=============================================================',/,&
'            NTCABS COURANT  = ',I10                           ,/,&
'            NTMABS RESET TO = ',I10                           ,/,&
'=============================================================',/,&
                                                                /)

#else

 1000 format(/,                                                   &
'=============================================================',/,&
'            NTCABS CURRENT  = ',I10                           ,/,&
'            NTMABS RESET TO = ',I10                           ,/,&
'=============================================================',/,&
                                                                /)

#endif

end subroutine
