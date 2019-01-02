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

subroutine tdesi1 &
!================

 ( jj     , nfabor , nn     , ifabor , iclass )

!===============================================================================
!  FONCTION
!  --------

!              SOUS-PROGRAMME DE DESCENTE D'UN ARBRE BINAIRE POUR
!               L'ALGORITHME DE TRI PAR ARBRE (HEAPSORT)

!              VERSION 2D

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! jj               ! e  ! <-- ! niveau de l'arbre binaire traite               !
! nn               ! e  ! <-- ! nombre de faces    actives                     !
! nfabor           ! e  ! <-- ! taille du tableau a descendre                  !
! iclass           ! e  ! <-- ! pointeurs sur les numeros d'elements           !
!                  !    !     ! dans l'arbre binaire                           !
! ifabor           ! tr ! <-- ! voisins des faces                              !
!  (  nfabor)      !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

integer          jj, nn, nfabor

integer          iclass(nfabor)
integer          ifabor(nfabor)

! VARAIBLES LOCALES

integer      ii, ll, isvnum
logical      desc, isup

!===============================================================================

ii = jj
ll = 2*jj

desc = .true.

if ((ll.gt.nn).or.(.not. desc)) goto 20

 10   continue

   if (ll .lt. nn) then

!--->       ISUP = .TRUE. SI ICLASS (LL) > ICLASS (LL+1)

      isup = .true.
      if ( ifabor (iclass (ll+1)) .gt.                            &
                                     ifabor (iclass (ll)) ) then
         isup = .false.
      endif

      if (isup) ll = ll + 1

   endif

!--->    ISUP = .TRUE. SI ICLASS (II) > ICLASS (LL)

   isup = .true.
   if ( ifabor (iclass (ll)) .gt.                                 &
                                    ifabor (iclass (ii)) ) then
      isup = .false.
   endif

   if (isup) then
      isvnum = iclass (ii)
      iclass (ii) = iclass (ll)
      iclass (ll) = isvnum
      ii = ll
      ll = 2*ll
      desc = .true.
   else
      desc = .false.
   endif

   if ((ll .le. nn) .and. desc) goto 10

 20   continue

return
end subroutine
