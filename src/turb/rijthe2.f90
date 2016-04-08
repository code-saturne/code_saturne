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


!===============================================================================
! Function:
! ---------
!> \file rijthe2.f90
!> \brief Gravity terms
!>        For \f$R_{ij}\f$
!>        \f[ var = R_{11} \: R_{22} \: R_{33} \:R_{12} \:R_{23} \:R_{13}\f]

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     gradro        work array for \f$ \grad{rom} \f$
!> \param[in,out] smbr          work array for second member
!______________________________________________________________________________!

subroutine rijthe2 &
 ( nscal  ,                                                       &
   gradro , smbr   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nscal

double precision gradro(3,ncelet)
double precision smbr(6,ncelet)

! Local variables

integer          iel, dimrij, isou

double precision uns3, const, kseps, csttmp
double precision prdtur, r1t, r2t, r3t
double precision g11, g22, g33, g12, g13, g23, gkks3
double precision phit11, phit22, phit33, phit12, phit13, phit23

double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:,:), pointer :: cvara_rij

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! ebrsm
if (iturb.eq.32) then
  csttmp = cebmr6
else
  csttmp = crij3
endif

if(iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

const = -1.5d0*cmu/prdtur
uns3  = 1.d0/3.d0

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_v(ivarfl(irij), cvara_rij)

call field_get_dim(ivarfl(irij), dimrij)

!===============================================================================
! 2. Terms for Rij:
!      rom*volume*dRij/dt =
!                     ... + (Gij - CRIJ3*(Gij-Delta ij Gkk/3))*volume
!            With Gij = -(1.5 cmu/prdtur) (k/eps) (Rit Gj + Rjt Gi)
!                 Rit = Rik drom/dxk (sum on k)
!===============================================================================
do isou = 1, dimrij
  if     (isou.eq.1) then

    do iel = 1, ncel

      r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
          + cvara_rij(4,iel)*gradro(2,iel)                            &
          + cvara_rij(6,iel)*gradro(3,iel)
      r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
          + cvara_rij(2,iel)*gradro(2,iel)                            &
          + cvara_rij(5,iel)*gradro(3,iel)
      r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
          + cvara_rij(5,iel)*gradro(2,iel)                            &
          + cvara_rij(3,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g11 = const*kseps*2.d0*(r1t*gx       )
      g22 = const*kseps*2.d0*(r2t*gy       )
      g33 = const*kseps*2.d0*(r3t*gz       )
      gkks3 = uns3*(g11+g22+g33)

      phit11 = -csttmp*(g11-gkks3)

      smbr(isou,iel) = smbr(isou,iel) + (g11+phit11)*volume(iel)

    enddo

  elseif (isou.eq.2) then

    do iel = 1, ncel

      r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
          + cvara_rij(4,iel)*gradro(2,iel)                            &
          + cvara_rij(6,iel)*gradro(3,iel)
      r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
          + cvara_rij(2,iel)*gradro(2,iel)                            &
          + cvara_rij(5,iel)*gradro(3,iel)
      r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
          + cvara_rij(5,iel)*gradro(2,iel)                            &
          + cvara_rij(3,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g11 = const*kseps*2.d0*(r1t*gx       )
      g22 = const*kseps*2.d0*(r2t*gy       )
      g33 = const*kseps*2.d0*(r3t*gz       )
      gkks3 = uns3*(g11+g22+g33)

      phit22 = -csttmp*(g22-gkks3)

      smbr(isou,iel) = smbr(isou,iel) + (g22+phit22)*volume(iel)

    enddo

  elseif (isou.eq.3) then

    do iel = 1, ncel

      r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
          + cvara_rij(4,iel)*gradro(2,iel)                            &
          + cvara_rij(6,iel)*gradro(3,iel)
      r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
          + cvara_rij(2,iel)*gradro(2,iel)                            &
          + cvara_rij(5,iel)*gradro(3,iel)
      r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
          + cvara_rij(5,iel)*gradro(2,iel)                            &
          + cvara_rij(3,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g11 = const*kseps*2.d0*(r1t*gx       )
      g22 = const*kseps*2.d0*(r2t*gy       )
      g33 = const*kseps*2.d0*(r3t*gz       )
      gkks3 = uns3*(g11+g22+g33)

      phit33 = -csttmp*(g33-gkks3)

      smbr(isou,iel) = smbr(isou,iel) + (g33+phit33)*volume(iel)

    enddo

  elseif (isou.eq.4) then

    do iel = 1, ncel

      r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
          + cvara_rij(4,iel)*gradro(2,iel)                            &
          + cvara_rij(6,iel)*gradro(3,iel)
      r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
          + cvara_rij(2,iel)*gradro(2,iel)                            &
          + cvara_rij(5,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g12 = const*kseps*     (r1t*gy+r2t*gx)

      phit12 = -csttmp* g12

      smbr(isou,iel) = smbr(isou,iel) + (g12+phit12)*volume(iel)

    enddo

  elseif (isou.eq.6) then

    do iel = 1, ncel

      r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
          + cvara_rij(4,iel)*gradro(2,iel)                            &
          + cvara_rij(6,iel)*gradro(3,iel)
      r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
          + cvara_rij(5,iel)*gradro(2,iel)                            &
          + cvara_rij(3,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g13 = const*kseps*     (r1t*gz+r3t*gx)

      phit13 = -csttmp* g13

      smbr(isou,iel) = smbr(isou,iel) + (g13+phit13)*volume(iel)

    enddo

  elseif (isou.eq.5) then

    do iel = 1, ncel

      r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
          + cvara_rij(2,iel)*gradro(2,iel)                            &
          + cvara_rij(5,iel)*gradro(3,iel)
      r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
          + cvara_rij(5,iel)*gradro(2,iel)                            &
          + cvara_rij(3,iel)*gradro(3,iel)

      kseps = (cvara_rij(1,iel)+cvara_rij(2,iel)+cvara_rij(3,iel))  &
             /(2.d0*cvara_ep(iel))

      g23 = const*kseps*(r2t*gz+r3t*gy)

      phit23 = -csttmp* g23

      smbr(isou,iel) = smbr(isou,iel) + (g23+phit23)*volume(iel)

    enddo


  endif
enddo

return

end subroutine rijthe2
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


!===============================================================================
! Function:
! ---------
!> \brief Gravity terms for \f$\epsilon\f$
!>

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     gradro        work array for \f$ \grad{rom} \f$
!> \param[in,out] smbr          work array for second member
!______________________________________________________________________________!

subroutine rijtheps &
 ( nscal  ,                                                       &
   gradro , smbr   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nscal

double precision gradro(3,ncelet)
double precision smbr(ncelet)

! Local variables

integer          iel

double precision uns3, const, csttmp
double precision prdtur, r1t, r2t, r3t
double precision g11p, g22p, g33p
double precision aa, bb

double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:,:), pointer :: cvara_rij

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! ebrsm
if (iturb.eq.32) then
  csttmp = cebmr6
else
  csttmp = crij3
endif

if(iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

const = -1.5d0*cmu/prdtur
uns3  = 1.d0/3.d0

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_v(ivarfl(irij), cvara_rij)

!===============================================================================
! 2. Terms for epsilon:
!      rom*volumr*deps/dt =
!                     ... + CEPS1*(EPS/K)*Max(0,(Gkk/2))*volume
!            With Gij = -(1.5 cmu/prdtur) (k/eps) (Rit Gj + Rjt Gi)
!                 Rit = Rik drom/dxk (sum on k)
!            We simplify (eps/k) by noting
!                GijP = -(1.5 cmu/prdtur)         (Rit Gj + Rjt Gi)
!      rom*volume*deps/dt =
!                     ... + CEPS1*        Max(0,(GkkP/2))*volume
!===============================================================================


do iel = 1, ncel

  r1t = cvara_rij(1,iel)*gradro(1,iel)                            &
      + cvara_rij(4,iel)*gradro(2,iel)                            &
      + cvara_rij(6,iel)*gradro(3,iel)
  r2t = cvara_rij(4,iel)*gradro(1,iel)                            &
      + cvara_rij(2,iel)*gradro(2,iel)                            &
      + cvara_rij(5,iel)*gradro(3,iel)
  r3t = cvara_rij(6,iel)*gradro(1,iel)                            &
      + cvara_rij(5,iel)*gradro(2,iel)                            &
      + cvara_rij(3,iel)*gradro(3,iel)

  g11p = const*2.d0*(r1t*gx)
  g22p = const*2.d0*(r2t*gy)
  g33p = const*2.d0*(r3t*gz)

  aa = 0.d0
  bb = 0.5d0*(g11p+g22p+g33p)
  smbr(iel) = smbr(iel) + ce1*max(aa,bb)*volume(iel)

enddo

return

end subroutine rijtheps
