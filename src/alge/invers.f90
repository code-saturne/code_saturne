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

subroutine invers &
!================

 ( cnom   , isym   , ibsize , iesize , ipol   , ireslp , nitmap , &
   imgrp  , ncymxp , nitmfp ,                                     &
   iwarnp , niterf , icycle , iinvpe ,                            &
   epsilp , rnorm  , residu ,                                     &
   dam    , xam    , smbrp  , vx     )

!===============================================================================
! Purpose:
! -------

! Call linear system resolution:
! - multigrid + -conjugate gradient or Jacobi or Bi-CGstab)
! - Jacobi
! - Bi-CGstab
! we assume vx in initialized on input.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! cnom             ! a  ! <-- ! variable name                                  !
! isym             ! e  ! <-- ! flag = 1: symmetric matrix                     !
!                  !    !     !        2: non-symmetric matrix                 !
! ibsize           !    ! <-- ! flag = 1: diag block size 1                    !
!                  !    !     !      = 3: diag block size 3                    !
! iesize           !    ! <-- ! flag = 1: extra diag block size                !
!                  !    !     !      = 3: diag block size 3                    !
! ipol             ! e  ! <-- ! polynomial degree for preconditioning          !
!                  !    !     !   (0 <-- diagonal)                             !
! ireslp           ! e  ! <-- ! solver type: 0 conjugate gradient              !
!                  !    !     !              1 Jacobi                          !
!                  !    !     !              2 CG-stab                         !
! nitmap           ! e  ! <-- ! max number of iterations for soultion          !
! imgrp            ! e  ! <-- ! 1 for multigrid, 0 otherwise                   !
! ncymxp           ! e  ! <-- ! max. number of multigrid cycles                !
! nitmfp           ! e  ! <-- ! number of equivalent iterations on fine mesh   !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! niterf           ! e  ! --> ! number of iterations done (non-multigrid)      !
! icycle           ! e  ! --> ! number of multigrid cycles done                !
! iinvpe           ! e  ! <-- ! flag to cancel increments in rotational        !
!                  !    !     ! periodicity (=2) or to exchange them normally  !
!                  !    !     ! in a scalar fashion (=1)                       !
! epsilp           ! r  ! <-- ! precision for iterative resolution             !
! rnorm            ! r  ! <-- ! residue normalization                          !
! residu           ! r  ! --> ! final non-normalized residue                   !
! dam(ncelet       ! tr ! <-- ! diagonal (fine mesh if multigrid)              !
! xam(nfac,isym    ! tr ! <-- ! extradiagonal (fine mesh if multigrid)         !
! smbrp(ncelet     ! tr ! <-- ! right hand side (fine mesh if multigrid)       !
! vx(ncelet)       ! tr ! <-- ! system solution                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use mesh

!===============================================================================

implicit none

! Arguments

character*16     cnom
integer          isym   , ipol   , ireslp , nitmap , ibsize , iesize
integer          imgrp  , ncymxp , nitmfp
integer          iwarnp
integer          niterf , icycle , iinvpe
double precision epsilp , rnorm  , residu

double precision dam(*), xam(*)
double precision smbrp(*)
double precision vx(*)

! Local variables

integer          lnom
integer          iresds, iresas, nitmds, nitmas, ilved

!===============================================================================

! Initialization

lnom = len(cnom)

icycle = 0
niterf = 0
ilved = 2

! xam and dam are interleaved if ibsize is greater than 1
if (ibsize.gt.1) ilved = 1

! Resolution

if (imgrp.eq.1) then

  iresds = ireslp
  iresas = ireslp

  nitmds = nitmfp
  nitmas = nitmfp

  call resmgr                                                     &
  !==========
 ( cnom   , lnom   ,                                              &
   iresds , iresas , ireslp , ipol   ,                            &
   ncymxp , nitmds , nitmas , nitmap , iinvpe ,                   &
   iwarnp , icycle , niterf , epsilp , rnorm  , residu ,          &
   smbrp  , vx     )

elseif(imgrp.eq.0) then

  call reslin &
  !==========
 ( cnom   , lnom   , ncelet , ncel   , nfac   ,                            &
   isym   , ilved  , ibsize , iesize , ireslp , ipol   , nitmap , iinvpe , &
   iwarnp , niterf , epsilp , rnorm  , residu ,                            &
   !        ------                     ------
   ifacel , dam    , xam    , smbrp  , vx     )
   !                          -----

endif

!----
! End
!----

return

end subroutine
