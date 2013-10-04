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

subroutine diffst &
!================

 ( nscal  ,                                              &
   rtp    , propce ,                                     &
   coefa  , coefb  )

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
! rtp,             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  at current time step                          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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

double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          ivar  , iel   , ifac  , iscal
integer          ipcvst
integer          nswrgp, imligp, iwarnp, iclvar, iclvaf
integer          iccocg, inc
integer          iconvp, idiffp, ircflp
integer          ischcp, isstpp, ippvar
integer          ipcvsl, iflmas, iflmab
integer          imucpp, idftnp
double precision epsrgp, climgp, extrap
double precision blencp, relaxp, thetex
double precision qimp , hint

integer          icvflb
integer          ivoid(1)

double precision rvoid(1)

double precision, allocatable, dimension(:) :: vistot, viscf, viscb
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp
double precision, allocatable, dimension(:) :: whsad
double precision, allocatable, dimension(:) :: xcpp
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================

! Memory allocation
allocate(vistot(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(xcpp(ncelet))

do iscal = 1, nscal

  ! Index for variable
  ivar = isca(iscal)

  if (iscalt.gt.0) then
    if (ivar.eq.isca(iscalt) .or. iscavr(iscal).eq.iscalt) then
      if (abs(iscsth(iscalt)).eq.1) then
        imucpp = 1
      else
        imucpp = 0
      endif
    else
      imucpp = 0
    endif
  else
    imucpp = 0
  endif

  if (imucpp.eq.0) then
    do iel = 1, ncel
      xcpp(iel) = 1.d0
    enddo
  elseif (imucpp.eq.1) then
    if (icp.gt.0) then
      do iel = 1, ncel
        xcpp(iel) = propce(iel,ipproc(icp))
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
  !
  ! Index for Boundary conditions
  iclvar = iclrtp(ivar,icoef)
  iclvaf = iclrtp(ivar,icoeff)

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
  ippvar = ipprtp(ivar)
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetex = 1.d0

  ! Pointers to the mass fluxes
  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  ! Diffusion velocity

  ! Index for molecular diffusivity
  if (ivisls(iscal).gt.0) then
    ipcvsl = ipproc(ivisls(iscal))
  else
    ipcvsl = 0
  endif

  ! Index for turbulent diffusivity
  ipcvst = ipproc(ivisct)

  if (idiff(ivar).ge.1) then

    ! Only the positive part of mu_t is considered (MAX(mu_t,0)),
    ! Dynamic LES can cause negative mu_t (clipping on (mu+mu_t))
    ! The positive part of (K+K_t) would have been considered
    ! but should allow negative K_t that is considered non physical here

    if(ipcvsl.eq.0)then
      do iel = 1, ncel
        vistot(iel) = visls0(iscal)                                     &
           + idifft(ivar)*xcpp(iel)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
      enddo
    else
      do iel = 1, ncel
        vistot(iel) = propce(iel,ipcvsl)                                &
           + idifft(ivar)*xcpp(iel)*max(propce(iel,ipcvst),zero)/sigmas(iscal)
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

  if (ippmod(icod3p).eq.1.and.iscal.eq.ihm) then
  ! if (nscapp.gt.0.and.iscal.eq.ihm) then

    ! Memory allocation
    allocate(coefap(nfabor), coefbp(nfabor))
    allocate(cofafp(nfabor), cofbfp(nfabor))
    allocate(whsad(ncelet))

    ! Hs is store in a local array and the source term is initialized
    do iel = 1, ncel
      whsad(iel) = propce(iel,ipproc(iustdy(ihm)))
      propce(iel,ipproc(iustdy(ihm))) = 0.d0
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

    ! all boundary convective flux with upwind
    icvflb = 0

    ! Diffusion term calculation
    call bilsca &
   !==========
  ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
    ischcp , isstpp , inc    , imrgra , iccocg ,                   &
    ippvar , iwarnp , imucpp , idftnp ,                            &
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
    call bilsca &
    !==========
  ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
    ischcp , isstpp , inc    , imrgra , iccocg ,                   &
    ippvar , iwarnp , imucpp , idftnp ,                            &
    blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
    rtp(1,ivar)     , rtp(1,ivar)     ,                            &
    coefa(1,iclvar) , coefb(1,iclvar) ,                            &
    coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
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
