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

subroutine grdcel &
!================

 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   pvar   , coefap , coefbp ,                                     &
   grad   )

!===============================================================================
! Purpose:
! --------

! Call different cell gradient subroutines

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivar             ! i  ! <-- ! numero de la variable                          !
!                  !    !     !   destine a etre utilise pour la               !
!                  !    !     !   periodicite uniquement (pering)              !
!                  !    !     !   on pourra donner ivar=0 si la                !
!                  !    !     !   variable n'est ni une composante de          !
!                  !    !     !   la vitesse, ni une composante du             !
!                  !    !     !   tenseur des contraintes rij                  !
! imrgra           ! i  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindres carres                             !
!                  !    !     !  2 moindres carres support etendu              !
!                  !    !     !    complet                                     !
!                  !    !     !  3 moindres carres avec selection du           !
!                  !    !     !    support etendu                              !
! inc              ! i  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! iccocg           ! i  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! i  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! i  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! nfecra           ! i  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar  (ncelet    ! ra ! <-- ! variable (pression)                            !
! coefap,coefbp    ! ra ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! grad(ncelet,3)   ! ra ! --> ! gradient de pvar                               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use albase
use cplsat
use pointe
use parall
use period
use mesh
use numvar, only: ivarfl

!===============================================================================

implicit none

! Arguments

integer          ivar   , imrgra , inc    , iccocg , nswrgp
integer          imligp ,iwarnp  , nfecra
double precision epsrgp , climgp , extrap


double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision grad(ncelet,3)

! Local variables

integer          iphydp, ipond, ilved
integer          idimtr

double precision rvoid(1)
double precision climin

!===============================================================================

!===============================================================================
! 0. Preparation for periodicity of rotation
!===============================================================================

! By default, the gradient will be treated as a vector ...
!   (i.e. we assume it is the gradient of a scalar field)

! If there is no rotation, halos are synchronized by synvec/synvni
!   (periodicity is implicit)

! If rotational periodicities are present,
!   we determine if the variable is a vector (velocity) or a tensor
!   (Reynolds stresses) so as to apply the necessary treatment.
!   We set idimtr and we retrieve the matching gradient.
! Note that if halo gradients have not been saved before, they cannot be
!   retrieved here (...)
!   So this subroutine is called by phyvar (in perinr)
!   to compute gradients at the beginning of the time step and save them
!   in dudxyz et drdxyz

! It is necessary for idimtr to always be initialized, even with no
!   periodicity of rotation, so it's default value is set.

idimtr = 0

if (iperot.eq.1.and.ivar.gt.0) then
  call pering(ivarfl(ivar), idimtr, grad(1,1), grad(1,2), grad(1,3))
  !==========
endif

!===============================================================================
! 1. Compute gradient
!===============================================================================

! This subroutine is never used to compute the pressure gradient

ilved = 0
iphydp = 0
ipond  = 0

call cgdcel &
!==========
 ( ivar   , imrgra , ilved  , inc    , iccocg , nswrgp ,          &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , rvoid  , coefap , coefbp ,                            &
   pvar   , rvoid  , grad   )

return
end subroutine
