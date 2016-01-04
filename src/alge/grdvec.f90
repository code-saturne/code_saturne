!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine grdvec &
!================

 ( ivar   , imrgra , inc    , nswrgp , imligp ,                   &
   iwarnp , epsrgp , climgp ,                                     &
   ilved  ,                                                       &
   pvar   , coefav , coefbv ,                                     &
   gradv )

!===============================================================================
! FONCTION :
! ----------

! APPEL DES DIFFERENTES ROUTINES DE CALCUL DE GRADIENT CELLULE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivar             ! e  ! <-- ! numero de la variable                          !
!                  !    !     !   destine a etre utilise pour la               !
!                  !    !     !   periodicite uniquement (pering)              !
!                  !    !     !   on pourra donner ivar=0 si la                !
!                  !    !     !   variable n'est ni une composante de          !
!                  !    !     !   la vitesse, ni une composante du             !
!                  !    !     !   tenseur des contraintes rij                  !
! imrgra           ! e  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindres carres                             !
!                  !    !     !  2 moindres carres support etendu              !
!                  !    !     !    complet                                     !
!                  !    !     !  3 moindres carres avec selection du           !
!                  !    !     !    support etendu                              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! pvar(3,ncelet)   ! tr ! <-- ! variable (vectorielle)                         !
! coefav,coefbv    ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (3,nfabor)     !    !     !  sur la normale a la face de bord              !
! gradv            ! tr ! --> ! gradient de la variable vectorielle            !
!   (3,3,ncelet)   !    !     !                                                !
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
use parall
use period
use mesh
use cstphy
use cstnum
use albase
use cplsat

!===============================================================================

implicit none

! Arguments

integer          ivar   , imrgra , inc    , nswrgp
integer          imligp , iwarnp
double precision epsrgp , climgp

double precision pvar(*)
double precision coefav(*), coefbv(*)
double precision gradv(3,3,ncelet)

logical ilved

! Local variables

integer          iel, isou

double precision, dimension(:,:), allocatable :: pvari

!===============================================================================

!===============================================================================
! 1. Computation of the gardient
!===============================================================================

! the velocity and the gradient fields are interleaved
if (ilved) then

  call cgdvec &
  !==========
 ( ivar   ,                                                       &
   imrgra , inc    , nswrgp , iwarnp , imligp , epsrgp , climgp , &
   coefav , coefbv , pvar   ,                                     &
   gradv  )

! We interleave the velocity
else

  !Allocation
  allocate(pvari(3,ncelet))

  do isou = 1, 3
    do iel = 1, ncelet
      pvari(isou,iel) = pvar(iel + (isou-1)*ncelet)
    enddo
  enddo

  call cgdvec &
  !==========
 ( ivar   ,                                                       &
   imrgra , inc    , nswrgp , iwarnp , imligp , epsrgp , climgp , &
   coefav , coefbv , pvari  ,                                     &
   gradv  )

  ! Free memory
  deallocate(pvari)

endif

return
end subroutine
