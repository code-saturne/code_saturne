!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine diffst &
!================

 ( nscal  ,                                              &
   propce )

!===============================================================================
! Function :
! --------

! Weakly compressible algorithm (semi-analytic):
!  Computation of scalar diffusion terms

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use albase
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nscal

double precision propce(ncelet,*)

! Local variables

integer          ivar  , iel   , ifac  , iscal, ivar0
integer          nswrgp, imligp, iwarnp
integer          iccocg, inc
integer          iconvp, idiffp, ircflp
integer          ischcp, isstpp
integer          ifcvsl, iflmas, iflmab
integer          imucpp, idftnp
double precision epsrgp, climgp, extrap
double precision blencp, relaxp, thetex
double precision qimp , hint

integer          icvflb
integer          ivoid(1)

double precision rvoid(1)

double precision, allocatable, dimension(:) :: vistot, viscf, viscb
double precision, allocatable, dimension(:) :: whsad
double precision, allocatable, dimension(:) :: xcpp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: visct, cpro_cp, cpro_viscls
double precision, dimension(:), pointer :: cvar_scal

!===============================================================================

! Memory allocation
allocate(vistot(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(xcpp(ncelet))

if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

do iscal = 1, nscal

  ! Index for variable
  ivar = isca(iscal)

  call field_get_val_s(ivarfl(isca(iscal)), cvar_scal)

  imucpp = 0
  if (iscavr(iscal).gt.0) then
    if (abs(iscacp(iscavr(iscal))).eq.1) then
      imucpp = 1
    endif
  else
    if (abs(iscacp(iscal)).eq.1) then
      imucpp = 1
    endif
  endif

  if (imucpp.eq.0) then
    do iel = 1, ncel
      xcpp(iel) = 1.d0
    enddo
  elseif (imucpp.eq.1) then
    if (icp.gt.0) then
      do iel = 1, ncel
        xcpp(iel) = cpro_cp(iel)
      enddo
   else
      do iel = 1, ncel
        xcpp(iel) = cp0
      enddo
    endif
  endif

  ! Handle parallelism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(xcpp)
  endif

  ivar0  = 0
  iconvp = 0
  idiffp = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  inc    = 1
  iccocg = 1
  idftnp = 1 !idften(ivar)!FIXME when activating GGDH
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetex = 1.d0
  ! all boundary convective flux with upwind
  icvflb = 0

  ! Pointers to the mass fluxes
  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  ! Diffusion velocity

  ! Index for molecular diffusivity
  call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_viscls)
  endif

  ! Index for turbulent diffusivity
  call field_get_val_s(iprpfl(ivisct), visct)

  if (idiff(ivar).ge.1) then

    ! Only the positive part of mu_t is considered (MAX(mu_t,0)),
    ! Dynamic LES can cause negative mu_t (clipping on (mu+mu_t))
    ! The positive part of (K+K_t) would have been considered
    ! but should allow negative K_t that is considered non physical here

    if(ifcvsl.lt.0)then
      do iel = 1, ncel
        vistot(iel) = visls0(iscal)                                     &
           + idifft(ivar)*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
      enddo
    else
      do iel = 1, ncel
        vistot(iel) = cpro_viscls(iel)                                  &
           + idifft(ivar)*xcpp(iel)*max(visct(iel),zero)/sigmas(iscal)
      enddo
    endif

    call viscfa ( imvisf , vistot , viscf , viscb )
    !==========

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo
    do iel = 1, ncel
      vistot(iel) = 0.d0
    enddo

  endif

  ! If the combustion model is used, gap between enthalpy and
  ! adiabatic enthalpy has to be used. The last one is already considered
  ! through the mixture fraction contribution.

  if (ippmod(icod3p).eq.1.and.iscal.eq.iscalt) then

    ! Memory allocation
    allocate(coefap(nfabor), coefbp(nfabor))
    allocate(cofafp(nfabor), cofbfp(nfabor))
    allocate(whsad(ncelet))

    ! Hs is store in a local array and the source term is initialized
    do iel = 1, ncel
      whsad(iel) = propce(iel,ipproc(iustdy(iscalt)))
      propce(iel,ipproc(iustdy(iscalt))) = 0.d0
    enddo

    ! Parallel and periodic exchanges
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(whsad)
    endif

    ! Boundary condition on Hsad: Homogenous Neumann
    do ifac = 1, nfabor
      iel = ifabor(ifac)

      hint = vistot(iel)/distb(ifac)
      qimp = 0.d0

      call set_neumann_scalar &
           !=================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )

    enddo

    ! Diffusion term calculation
    call bilsca &
   !==========
  ( idtvar , ivar0  , iconvp , idiffp , nswrgp , imligp , ircflp , &
    ischcp , isstpp , inc    , imrgra , iccocg ,                   &
    iwarnp , imucpp , idftnp ,                                     &
    blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
    whsad  , whsad  ,                                              &
    coefap , coefbp ,                                              &
    cofafp , cofbfp ,                                              &
    imasfl , bmasfl ,                                              &
    viscf  , viscb  , rvoid  , rvoid  ,                            &
    rvoid  , rvoid  ,                                              &
    icvflb , ivoid  ,                                              &
    propce(1,ipproc(iustdy(iscal))) )

    ! Free memory
    deallocate(coefap, coefbp)
    deallocate(cofafp, cofbfp)
    deallocate(whsad)

  else

    ! Diffusion term calculation

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)
    call field_get_coefaf_s(ivarfl(ivar), cofafp)
    call field_get_coefbf_s(ivarfl(ivar), cofbfp)

    call bilsca &
    !==========
  ( idtvar , ivar0  , iconvp , idiffp , nswrgp , imligp , ircflp , &
    ischcp , isstpp , inc    , imrgra , iccocg ,                   &
    iwarnp , imucpp , idftnp ,                                     &
    blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
    cvar_scal       , cvar_scal       ,                            &
    coefap , coefbp , cofafp , cofbfp ,                            &
    imasfl , bmasfl ,                                              &
    viscf  , viscb  , rvoid  , xcpp   ,                            &
    rvoid  , rvoid  ,                                              &
    icvflb , ivoid  ,                                              &
    propce(1,ipproc(iustdy(iscal))) )

  endif

enddo

! Free memory
deallocate(viscf, viscb)
deallocate(vistot)
deallocate(xcpp)

!----
! Formats
!----

!----
! End
!----

return
end subroutine
