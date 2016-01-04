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

subroutine grdpre &
!================

 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   pvar   , ktvar  , coefap , coefbp ,                            &
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
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar  (ncelet    ! ra ! <-- ! variable (pression)                            !
! ktvar (ncelet    ! ra ! <-- ! variable (coefficient du gradient de pression) !
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

!===============================================================================

implicit none

! Arguments

integer          ivar   , imrgra , inc    , iccocg , nswrgp
integer          imligp , iwarnp
double precision epsrgp , climgp , extrap


double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision ktvar(ncelet)
double precision grad(ncelet,3)

! Local variables

integer          iphydp, ipond, ilved
integer          idimtr

double precision rvoid(1)

!===============================================================================

!===============================================================================
! 1. Compute gradient
!===============================================================================

ilved = 0

! The current variable is a scalar (and the gradient is a vector)
idimtr = 0

! In apriori the hydrostatic pressure gradient is computed without
! extern hydrostatic force

iphydp = 0

! the pressure gradient coefficient ponderation activated
ipond = 1

call cgdcel &
!==========
 ( ivar   , imrgra , ilved  , inc    , iccocg , nswrgp ,          &
   idimtr , iphydp , ipond  , iwarnp , imligp , epsrgp , extrap , &
   climgp , rvoid  , coefap , coefbp ,                            &
   pvar   , ktvar  , grad   )

return
end subroutine
