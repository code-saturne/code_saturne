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

!> \file csc2cl.f90
!> \brief Translation of the "itypfb(*, *) = icscpl" condition.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvcp
!> \param[in]     nvcpto
!> \param[in]     nfbcpl
!> \param[in]     nfbncp
!> \param[out]    icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in]     lfbcpl
!> \param[in]     lfbncp
!> \param[out]    itypfb        boundary face types
!> \param[in]     dt            time step (per cell)
!> \param[out]    rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!> \param[in]     rvcpfb
!> \param[in]     pndcpl
!> \param[in]     dofcpl
!______________________________________________________________________________

subroutine csc2cl &
 ( nvcp   , nvcpto , nfbcpl , nfbncp ,                            &
   icodcl , itypfb ,                                              &
   lfbcpl , lfbncp ,                                              &
   dt     ,                                                       &
   rcodcl ,                                                       &
   rvcpfb , pndcpl , dofcpl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use entsor
use parall
use period
use cplsat
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          nvcp   , nvcpto
integer          nfbcpl , nfbncp

integer          icodcl(nfabor,nvar)
integer          lfbcpl(nfbcpl)  , lfbncp(nfbncp)
integer          itypfb(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)
double precision rvcpfb(nfbcpl,nvcpto), pndcpl(nfbcpl)
double precision dofcpl(3,nfbcpl)

! Local variables

integer          ifac, iel, isou, icscp
integer          inc, iccocg, iprev
integer          ivar, iflmab
integer          ipt

double precision xip   , xiip  , yiip  , ziip
double precision xjp
double precision xipf, yipf, zipf, ipf

double precision xif, yif, zif, xopf, yopf, zopf
double precision gradi, pondj, flumab

double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: grad

double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_var

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

! Pointer to the boundary mass flux
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

if (ifaccp.eq.0) then
  icscp = icscpl
else
  icscp = icscpd
endif

!===============================================================================
! 1.  Translation of the coupling to boundary conditions
!===============================================================================

! Allocate a temporary array for gradient computation
allocate(grad(3,ncelet))
allocate(gradv(3,3,ncelet))

! Reminder: variables are received in the order of VARPOS;
! looping on variables is thus sufficient.

do ivar = 1, nvcp

  ! --- Compute gradient of variable if it is interpolated.
  !       Exchanges for parallelism and periodicity have already been
  !       done in CSCPFB. Non need to do them again.

  iprev  = 0
  inc    = 1
  iccocg = 1

  isou = 0 ! set for iu, iv, iw only

  if (ivar.ne.iu .and. ivar.ne.iv .and. ivar.ne.iw) then

    call field_get_val_s(ivarfl(ivar), cvar_var)
    call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,   &
                               iccocg,                             &
                               grad)

  else if (ivar.eq.iu) then

    call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,  &
                               gradv)

  endif

  ! For a specific face to face coupling, geometric assumptions are made

  if (ifaccp.eq.1) then


    do ipt = 1, nfbcpl

      ifac = lfbcpl(ipt)
      iel  = ifabor(ifac)

      ! Information for current instance interpolated at I'
      xiip = diipb(1,ifac)
      yiip = diipb(2,ifac)
      ziip = diipb(3,ifac)

      xif = cdgfbo(1,ifac) -xyzcen(1,iel)
      yif = cdgfbo(2,ifac) -xyzcen(2,iel)
      zif = cdgfbo(3,ifac) -xyzcen(3,iel)

      xipf = cdgfbo(1,ifac)-xiip - xyzcen(1,iel)
      yipf = cdgfbo(2,ifac)-yiip - xyzcen(2,iel)
      zipf = cdgfbo(3,ifac)-ziip - xyzcen(3,iel)

      ipf = sqrt(xipf**2+yipf**2+zipf**2)

      xopf = dofcpl(1,ipt)
      yopf = dofcpl(2,ipt)
      zopf = dofcpl(3,ipt)

      if (ivar.eq.ipr) then

        ! --- We want to prescribe a Dirichlet for pressure so as to conserve
        !     the pressure gradient through the coupling and remain consistent
        !     with the resolution of the pressure gradient on an orthogonal mesh.

        xip = cvar_var(iel) &
            + (grad(1,iel)*xiip + grad(2,iel)*yiip + grad(3,iel)*ziip)

      else if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw) then

        isou = ivar - iu + 1

        ! --- For all other variables, we want to prescribe a Dirichlet matching
        !     the convective fluxes at the center. We resrve a choice between
        !     UPWIND, SOLU, and CENTERED. Only the centered case respects the diffusion
        !     at the domain's interior faces. For UPWIND and SOLU, the decentering
        !     is done here and in bilsc2.f90 for coupled faces.

        ! -- UPWIND

        !        xip =  cvar_var(iel)

        ! -- SOLU

        !        xip =   cvar_var(iel)                &
        !              + (   gradv(1,isou,iel)*xif    &
        !                 +  gradv(2,isou,iel)*yif    &
        !                 +  gradv(3,isou,iel)*zif)

        ! -- CENTERED

        xip =   vel(isou,iel)                &
              + (  gradv(1,isou,iel)*xiip    &
                 + gradv(2,isou,iel)*yiip    &
                 + gradv(3,isou,iel)*ziip)

      else

        ! -- UPWIND

        !        xip =  cvar_var(iel)

        ! -- SOLU

        !        xip = cvar_var(iel) &
        !            + (grad(1,iel)*xif + grad(2,iel)*yif + grad(3,iel)*zif)

        ! -- CENTERED

        xip =  cvar_var(iel) &
            + (grad(1,iel)*xiip + grad(2,iel)*yiip + grad(3,iel)*ziip)

      endif

      ! -- We need alpha_ij for centered interpolation and flumab for decentering

      pondj = pndcpl(ipt)
      flumab = bmasfl(ifac)

      ! Information received from distant instance at J'/O'
      xjp = rvcpfb(ipt,ivar)

      itypfb(ifac) = icscp

      if (ivar.eq.ipr) then

        icodcl(ifac,ivar  ) = 1
        rcodcl(ifac,ivar,1) = (1.d0-pondj)*xjp + pondj*xip + p0

      else if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw) then

        icodcl(ifac,ivar  ) = 1

        ! -- DECENTERED (SOLU or UPWIND)

        !        if (flumab.ge.0.d0) then
        !          rcodcl(ifac,ivar,1) = xip
        !        else
        !          rcodcl(ifac,ivar,1) = xjp
        !        endif

        ! -- CENTERED

        rcodcl(ifac,ivar,1) = (1.d0-pondj)*xjp + pondj*xip

      else

        icodcl(ifac,ivar  ) = 1

        ! -- DECENTERED (SOLU or UPWIND)

        !        if(flumab.ge.0.d0) then
        !          rcodcl(ifac,ivar,1) = xip
        !        else
        !          rcodcl(ifac,ivar,1) = xjp
        !        endif

        ! -- CENTERED

        rcodcl(ifac,ivar,1) = (1.d0-pondj)*xjp + pondj*xip

      endif

    enddo

    ! For a generic coupling, no assumption can be made

  else


    ! --- Translation in terms of boundary conditions for located boundary faces
    !     --> Dirichlet BC type

    do ipt = 1, nfbcpl

      ifac = lfbcpl(ipt)
      iel  = ifabor(ifac)

      ! Information from local instance interpolated at I'
      xiip = diipb(1,ifac)
      yiip = diipb(2,ifac)
      ziip = diipb(3,ifac)

      xif = cdgfbo(1,ifac) -xyzcen(1,iel)
      yif = cdgfbo(2,ifac) -xyzcen(2,iel)
      zif = cdgfbo(3,ifac) -xyzcen(3,iel)

      xipf = cdgfbo(1,ifac)-xiip - xyzcen(1,iel)
      yipf = cdgfbo(2,ifac)-yiip - xyzcen(2,iel)
      zipf = cdgfbo(3,ifac)-ziip - xyzcen(3,iel)

      ipf = sqrt(xipf**2+yipf**2+zipf**2)

      xopf = dofcpl(1,ipt)
      yopf = dofcpl(2,ipt)
      zopf = dofcpl(3,ipt)

      ! Local information interpolated at I'/O'

      if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw) then

        isou = ivar - iu + 1

        xip =   vel(isou,iel)           &
              + gradv(1,isou,iel)*xiip  &
              + gradv(2,isou,iel)*yiip  &
              + gradv(3,isou,iel)*ziip

        gradi =  (  gradv(1,isou,iel)*xipf          &
                  + gradv(2,isou,iel)*yipf          &
                  + gradv(3,isou,iel)*zipf) / ipf

      else

        xip =   cvar_var(iel)            &
              + grad(1,iel)*(xiip+xopf)  &
              + grad(2,iel)*(yiip+yopf)  &
              + grad(3,iel)*(ziip+zopf)

        gradi = (grad(1,iel)*xipf+grad(2,iel)*yipf+grad(3,iel)*zipf)/ipf

      endif

      ! Information received from distant instance at J'/O'
      xjp = rvcpfb(ipt,ivar)

      itypfb(ifac) = icscp

      if (ivar.ne.ipr) then
        icodcl(ifac,ivar  ) = 1
        rcodcl(ifac,ivar,1) = 0.5d0*(xip+xjp)
      else
        icodcl(ifac,ivar  ) = 3
        rcodcl(ifac,ivar,3) = -0.5d0*dt(iel)*(gradi+xjp)
      endif


    enddo

  endif

  ! --- Non-located boundary faces
  !     --> Homogeneous Neuman BC type

  do ipt = 1, nfbncp

    ifac = lfbncp(ipt)

    itypfb(ifac) = icscp

    icodcl(ifac,ivar  ) = 3
    rcodcl(ifac,ivar,3) = 0.d0

  enddo

enddo

! Free memory
deallocate(gradv)
deallocate(grad)

!----
! Formats
!----

!----
! End
!----

return
end subroutine csc2cl
