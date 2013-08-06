!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine alemaj &
!================

 ( itrale ,                                                       &
   nvar   , nscal  ,                                              &
   impale ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , depale , xyzno0 )

!===============================================================================
! FONCTION :
! ----------

! MISE A JOUR DU MAILLAGE EN METHODE ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! impale(nnod)     ! te ! <-- ! indicateur de delacement impose                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! depale(3,nnod)   ! tr ! <-- ! deplacement aux noeuds                         !
! xyzno0(3,nnod)   ! tr ! <-- ! coordonnees noeuds maillage initial            !
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
use entsor
use cstphy
use cstnum
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          itrale
integer          nvar   , nscal

integer          impale(nnod)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision depale(3,nnod), xyzno0(3,nnod)

! Local variables

integer          inod
integer          iel
integer          icluma, iclvma, iclwma, idim

double precision, allocatable, dimension(:,:) :: dproj

!===============================================================================


!===============================================================================
! 1.  INITIALISATION
!===============================================================================

if(iwarni(iuma).ge.1) then
  write(nfecra,1000)
endif

icluma = iclrtp(iuma, icoef)
iclvma = iclrtp(ivma, icoef)
iclwma = iclrtp(iwma, icoef)

!===============================================================================
! 2.  MISE A JOUR DE LA GEOMETRIE
!===============================================================================
! (en utilisant le deplacement predit)

!     Projection du deplacement calcule sur les noeuds

! Allocate a temporary array
allocate(dproj(nnod,3))

call aldepl                                                       &
!==========
 ( ifacel , ifabor , ipnfac , nodfac , ipnfbr , nodfbr ,          &
   rtp(1,iuma), rtp(1,ivma), rtp(1,iwma),                         &
   coefa(1,icluma), coefa(1,iclvma), coefa(1,iclwma),             &
   coefb(1,icluma), coefb(1,iclvma), coefb(1,iclwma),             &
   dt, dproj  )

!     Mise a jour du deplacement sur les noeuds ou on ne l'a pas impose
!       (DEPALE a alors la valeur du deplacement au pas de temps precedent)

do inod = 1, nnod
  if (impale(inod).eq.0) then
    do idim = 1, 3
      depale(idim,inod) = depale(idim,inod) + dproj(inod,idim)
    enddo
  endif
enddo

! Free memory
deallocate(dproj)

!     Mise a jour de la geometrie

do inod = 1, nnod
  do idim = 1, ndim
    xyznod(idim,inod) = xyzno0(idim,inod) + depale(idim,inod)
  enddo
enddo

call algrma
!==========

! Abort at the end of the current time-step if there is a negative volume
if (volmin.le.0.d0) ntmabs = ntcabs


! Si on est a l'iteration d'initialisation, on remet les vitesses de maillage
!   a leur valeur initiale
if (itrale.eq.0) then
  do iel = 1, ncelet
    rtp(iel,iuma) = rtpa(iel,iuma)
    rtp(iel,ivma) = rtpa(iel,ivma)
    rtp(iel,iwma) = rtpa(iel,iwma)
  enddo
endif
!--------
! FORMATS
!--------

 1000 format(/,                                                   &
' ------------------------------------------------------------',/,&
                                                              /,/,&
'  MISE A JOUR DU MAILLAGE (ALE)                              ',/,&
'  =============================                              ',/)

!----
! FIN
!----

end subroutine
