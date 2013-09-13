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
!> \param[out]    icodcl        boundary condition code at boundary faces
!>                               - 1   -> dirichlet
!>                               - 3   -> flux density
!>                               - 4   -> sliding and u.n=0 (velocity)
!>                               - 5   -> friction and u.n=0 (velocity)
!>                               - 9   -> free inlet/outlet (inlet velocity
!>                              possibly fixed)
!>                               - 10  -> free inlet/outlet (possible inlet
!>                              volocity not fixed: prescribe a Dirichlet
!>                              value for scalars k, eps, scal in addition to
!>                              the usual Neumann
!> \param[in]     lfbcpl
!> \param[in]     lfbncp
!> \param[out]    itypfb        boundary face types
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                              (at current time step)
!> \param[in]     coefa         boundary conditions
!> \param[in]     coefb         boundary conditions
!> \param[out]    rcodcl        value of boundary conditions at boundary faces
!>                              rcodcl(1) = Dirichlet value
!>                              rcodcl(2) = ext. exchange coefficient value
!>                              (infinite if no exchange)
!>                              rcodcl(3) = value of the flux density
!>                              (negative if gain) in w/m2
!>                              for velocities:   (vistl+visct)*gradu
!>                              for pressure:                dt*gradp
!>                              for scalars:
!>                              cp*(viscls+visct/sigmas)*gradt
!> \param[in]     rvcpfb
!> \param[in]     pndcpl
!> \param[in]     dofcpl
!______________________________________________________________________________

subroutine csc2cl &
 ( nvcp   , nvcpto , nfbcpl , nfbncp ,                            &
   icodcl , itypfb ,                                              &
   lfbcpl , lfbncp ,                                              &
   dt     , rtp    ,                                              &
   coefa  , coefb  , rcodcl ,                                     &
   rvcpfb , pndcpl , dofcpl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvcp   , nvcpto
integer          nfbcpl , nfbncp

integer          icodcl(nfabor,nvarcl)
integer          lfbcpl(nfbcpl)  , lfbncp(nfbncp)
integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)
double precision rvcpfb(nfbcpl,nvcpto), pndcpl(nfbcpl)
double precision dofcpl(3,nfbcpl)

! Local variables

integer          ifac, iel
integer          inc, iccocg, iclvar, nswrgp, imligp
integer          iwarnp, ivar, iflmab
integer          ipt

double precision epsrgp, climgp, extrap
double precision xip   , xiip  , yiip  , ziip
double precision xjp
double precision xipf, yipf, zipf, ipf

double precision xif, yif, zif, xopf, yopf, zopf
double precision gradi, pondj, flumab

double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: bmasfl

!===============================================================================


!===============================================================================
! 1.  Translation of the coupling to boundary conditions
!===============================================================================

! Pointer to the boundary mass flux
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

! Allocate a temporary array for gradient computation
allocate(grad(ncelet,3))

! Reminder: variables are received in the order of VARPOS;
! loopin on variables is thus sufficient.

do ivar = 1, nvcp

  ! --- Compute gradient of variable if it is interpolated.
  !       Exchanges for parallelism and periodicity have already been
  !       done in CSCPFB. Non need to do them again.

  inc    = 1
  iccocg = 1

  iclvar = iclrtp(ivar,icoef)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                     &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
   grad   )


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

        ! --- We want to prescribe a Direchlet for pressure so as to conserve
        !     the pressure gradient through the coupling and remain consistent
        !     with the resolution of the pressure gradient on an orthogonal mesh.

        xip = rtp(iel,ivar) &
            + (grad(iel,1)*xiip + grad(iel,2)*yiip + grad(iel,3)*ziip)

      else if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw) then

        ! --- For all other variables, we want to prescribe a Dirichlet matching
        !     the convective fluxes at the center. We resrve a choice between
        !     UPWIND, SOLU, and CENTERED. Only the centered case respects the diffusion
        !     at the domain's interior faces. For UPWIND and SOLU, the decentering
        !     is done here and in bilsc2.f90 for coupled faces.

        ! -- UPWIND

        !        xip =  rtp(iel,ivar)

        ! -- SOLU

        !        xip =  rtp(iel,ivar) &
        !            + (grad(iel,1)*xif + grad(iel,2)*yif + grad(iel,3)*zif)

        ! -- CENTERED

        xip = rtp(iel,ivar) &
            + (grad(iel,1)*xiip + grad(iel,2)*yiip + grad(iel,3)*ziip)

      else

        ! -- UPWIND

        !        xip =  rtp(iel,ivar)

        ! -- SOLU

        !        xip = rtp(iel,ivar) &
        !            + (grad(iel,1)*xif + grad(iel,2)*yif + grad(iel,3)*zif)

        ! -- CENTERED

        xip =  rtp(iel,ivar) &
            + (grad(iel,1)*xiip + grad(iel,2)*yiip + grad(iel,3)*ziip)

      endif

      ! -- We need alpha_ij for centered interpolation and flumab for decentering

      pondj = pndcpl(ipt)
      flumab = bmasfl(ifac)

      ! Information received from distant instance at J'/O'
      xjp = rvcpfb(ipt,ivar)


      itypfb(ifac)  = icscpl

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

      xip =  rtp(iel,ivar) &
          + grad(iel,1)*(xiip+xopf) &
          + grad(iel,2)*(yiip+yopf) &
          + grad(iel,3)*(ziip+zopf)

      ! Information received from distant instance at J'/O'
      xjp = rvcpfb(ipt,ivar)


      gradi = (grad(iel,1)*xipf+grad(iel,2)*yipf+grad(iel,3)*zipf)/ipf

      itypfb(ifac)  = icscpl

      if(ivar.ne.ipr) then
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

    itypfb(ifac)  = icscpl

    icodcl(ifac,ivar  ) = 3
    rcodcl(ifac,ivar,3) = 0.d0

  enddo

enddo

! Free memory
deallocate(grad)

!----
! Formats
!----

!----
! End
!----

return
end subroutine
