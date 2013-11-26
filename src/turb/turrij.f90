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

!===============================================================================
! Function:
! ---------

!> \file turrij.f90
!>
!> \brief Solving the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
!>  slightly compressible flows for one time step.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of the ncepdp cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        mass source type for the variables (cf. ustsma)
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp           calculated variables at cell centers
!>                               (at the current time step)
!> \param[in]     rtpa          calculated variables at cell centers
!>                               (at the previous time step)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     tslagr        coupling term of the lagangian module
!> \param[in]     coefa, coefb  boundary conditions
!>
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        values of the variables associated to the
!>                               mass source
!>                               (for ivar=ipr, smacel is the mass flux)
!_______________________________________________________________________________

subroutine turrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr ,                                                       &
   coefa  , coefb  , ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use lagdim, only: ntersl
use numvar
use entsor
use cstphy
use optcal
use lagran
use pointe, only: coefau, coefbu
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp)

integer, dimension(ncesmp,nvar), target :: itypsm

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6)

double precision, dimension(ncesmp,nvar), target ::  smacel
double precision, dimension(ncelet,ntersl), target :: tslagr

! Local variables

integer          ifac  , iel   , ivar  , isou  , ii
integer          inc   , iccocg
integer          ipp   , iwarnp, iclip
integer          icliup, iclivp, icliwp
integer          nswrgp, imligp
integer          ipcrom, ipbrom, ipcroo, ipbroo, iivar
integer          iitsla
double precision epsrgp, climgp, extrap

logical          ilved

integer,          dimension(1), target :: ivoid(1)
double precision, dimension(1), target :: rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:,:,:) :: grdvel
double precision, allocatable, dimension(:,:) :: produc
double precision, allocatable, dimension(:,:) :: gradu, gradv, gradw, gradro

integer,          pointer, dimension(:) :: itpsmp
double precision, pointer, dimension(:) :: smcelp, gammap, tslage, tslagi

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

itpsmp => ivoid1

smcelp => rvoid1
gammap => rvoid1
tslage => rvoid1
tslagi => rvoid1

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(ncelet), rovsdt(ncelet))
allocate(grdvel(ncelet,3,3))

! Allocate other arrays, depending on user options
if (iturb.eq.30) then
  allocate(produc(6,ncelet))
endif

icliup = iclrtp(iu,icoef)
iclivp = iclrtp(iv,icoef)
icliwp = iclrtp(iw,icoef)

ipcrom = ipproc(irom  )
ipbrom = ipprob(irom  )

if(iwarni(iep).ge.1) then
  if (iturb.eq.30) then
    write(nfecra,1000)
  elseif (iturb.eq.31) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

!===============================================================================
! 2.1 Compute the velocity gradient
! WARNING: grdvel(iel, xyz, uvw)
!===============================================================================

iccocg = 1
inc    = 1

nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(iu)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

if (ivelco.eq.1) then

  ilved = .false.

  call grdvec &
  !==========
( iu     , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  ilved  ,                                                       &
  rtpa(1,iu) ,  coefau , coefbu,                                 &
  grdvel  )

else

  call grdvni &
  !==========
( iu  , imrgra , inc    , iccocg , nswrgp , imligp ,             &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
  rtpa(1,iu)   , coefa(1,icliup) , coefb(1,icliup) ,             &
  grdvel )

endif

!===============================================================================
! 2.2 Compute the production term for Rij LRR (iturb =30)
!===============================================================================

if (iturb.eq.30) then

  do ii = 1 , 6
    do iel = 1, ncel
      produc(ii,iel) = 0.0d0
    enddo
  enddo

  do iel = 1 , ncel

    ! grad u

    produc(1,iel) = produc(1,iel)                                &
                  - 2.0d0*(rtpa(iel,ir11)*grdvel(iel,1,1) +      &
                           rtpa(iel,ir12)*grdvel(iel,2,1) +      &
                           rtpa(iel,ir13)*grdvel(iel,3,1) )

    produc(4,iel) = produc(4,iel)                                 &
                  - (rtpa(iel,ir12)*grdvel(iel,1,1) +             &
                     rtpa(iel,ir22)*grdvel(iel,2,1) +             &
                     rtpa(iel,ir23)*grdvel(iel,3,1) )

    produc(5,iel) = produc(5,iel)                                 &
                  - (rtpa(iel,ir13)*grdvel(iel,1,1) +             &
                     rtpa(iel,ir23)*grdvel(iel,2,1) +             &
                     rtpa(iel,ir33)*grdvel(iel,3,1) )

    ! grad v

    produc(2,iel) = produc(2,iel)                                 &
                  - 2.0d0*(rtpa(iel,ir12)*grdvel(iel,1,2) +       &
                           rtpa(iel,ir22)*grdvel(iel,2,2) +       &
                           rtpa(iel,ir23)*grdvel(iel,3,2) )

    produc(4,iel) = produc(4,iel)                                 &
                  - (rtpa(iel,ir11)*grdvel(iel,1,2) +             &
                     rtpa(iel,ir12)*grdvel(iel,2,2) +             &
                     rtpa(iel,ir13)*grdvel(iel,3,2) )

    produc(6,iel) = produc(6,iel)                                 &
                  - (rtpa(iel,ir13)*grdvel(iel,1,2) +             &
                     rtpa(iel,ir23)*grdvel(iel,2,2) +             &
                     rtpa(iel,ir33)*grdvel(iel,3,2) )

    ! grad w

    produc(3,iel) = produc(3,iel)                                 &
                  - 2.0d0*(rtpa(iel,ir13)*grdvel(iel,1,3) +       &
                           rtpa(iel,ir23)*grdvel(iel,2,3) +       &
                           rtpa(iel,ir33)*grdvel(iel,3,3) )

    produc(5,iel) = produc(5,iel)                                 &
                  - (rtpa(iel,ir11)*grdvel(iel,1,3) +             &
                     rtpa(iel,ir12)*grdvel(iel,2,3) +             &
                     rtpa(iel,ir13)*grdvel(iel,3,3) )

    produc(6,iel) = produc(6,iel)                                 &
                  - (rtpa(iel,ir12)*grdvel(iel,1,3) +             &
                     rtpa(iel,ir22)*grdvel(iel,2,3) +             &
                     rtpa(iel,ir23)*grdvel(iel,3,3) )

  enddo

endif

!===============================================================================
! 3. Compute the density gradient for buoyant terms
!===============================================================================

if(igrari.eq.1) then

  ! Allocate a temporary array for the gradient calculation
  allocate(gradro(ncelet,3))

! Conditions aux limites : Dirichlet ROMB
!   On utilise VISCB pour stocker le coefb relatif a ROM
!   On impose en Dirichlet (COEFA) la valeur ROMB

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

! Le choix ci dessous a l'avantage d'etre simple

  nswrgp = nswrgr(ir11)
  imligp = imligr(ir11)
  iwarnp = iwarni(ir11)
  epsrgp = epsrgr(ir11)
  climgp = climgr(ir11)
  extrap = extrag(ir11)

  iivar = 0

!     Si on extrapole les termes sources et rho, on utilise cpdt rho^n
  ipcroo = ipcrom
  ipbroo = ipbrom
  if(isto2t.gt.0.and.iroext.gt.0) then
    ipcroo = ipproc(iroma)
    ipbroo = ipprob(iroma)
  endif

  call grdcel &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   propce(1,ipcroo), propfb(1,ipbroo), viscb           ,          &
   gradro )

endif

!===============================================================================
! 4.  Boucle sur les variables Rij (6 variables)
!     L'ordre est R11 R22 R33 R12 R13 R23 (La place de ces variables
!     est IR11.    ..
!     On resout les equation dans une routine semblable a covofi.f90
!===============================================================================

do isou = 1, 6
  if    (isou.eq.1) then
    ivar   = ir11
  elseif(isou.eq.2) then
    ivar   = ir22
  elseif(isou.eq.3) then
    ivar   = ir33
  elseif(isou.eq.4) then
    ivar   = ir12
  elseif(isou.eq.5) then
    ivar   = ir13
  elseif(isou.eq.6) then
    ivar   = ir23
  endif
  ipp    = ipprtp(ivar)

  if (iilagr.eq.2) then
    iitsla = itsr11 + (isou-1)
    tslage => tslagr(1:ncelet,iitsla)
    tslagi => tslagr(1:ncelet,itsli)
  endif

  if (ncesmp.gt.0) then
    itpsmp => itypsm(1:ncesmp,ivar)
    smcelp => smacel(1:ncesmp,ivar)
    gammap => smacel(1:ncesmp,ipr)
  endif

  ! Rij-epsilon standard (LRR)
  if (iturb.eq.30) then
    call resrij &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , produc , gradro ,                            &
   ckupdc , smcelp , gammap ,                                     &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbr   , rovsdt )

  ! Rij-epsilon SSG or EBRSM
  elseif (iturb.eq.31.or.iturb.eq.32) then

    call resssg &
    !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvel , gradro ,                            &
   ckupdc , smcelp , gammap ,                                     &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbr   , rovsdt )
  endif

enddo

!===============================================================================
! 5. Solve Epsilon
!===============================================================================

ivar   = iep
ipp    = ipprtp(ivar)
isou   = 7

if (ncesmp.gt.0) then
  itpsmp => itypsm(1:ncesmp,ivar)
  smcelp => smacel(1:ncesmp,ivar)
  gammap => smacel(1:ncesmp,ipr)
endif

call reseps &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvel , produc , gradro ,                   &
   ckupdc , smcelp , gammap ,                                     &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   smbr   , rovsdt )

!===============================================================================
! 6. Clipping
!===============================================================================

if (iturb.eq.32) then
  iclip = 1
else
  iclip = 2
endif

call clprij &
!==========
 ( ncelet , ncel   , nvar   ,                                     &
   iclip  ,                                                       &
   propce , rtpa   , rtp    )


! Free memory
deallocate(viscf, viscb)
deallocate(smbr, rovsdt)
if (allocated(gradro)) deallocate(gradro)
if (allocated(produc)) deallocate(produc)
deallocate(grdvel)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                      &
'   ** RESOLUTION DU Rij-EPSILON LRR'             ,/,&
'      -----------------------------'             ,/)
 1001 format(/,                                      &
'   ** RESOLUTION DU Rij-EPSILON SSG'             ,/,&
'      -----------------------------'             ,/)
 1002 format(/,                                      &
'   ** RESOLUTION DU Rij-EPSILON EBRSM'           ,/,&
'      -------------------------------'           ,/)

#else

 1000 format(/,                                      &
'   ** SOLVING Rij-EPSILON LRR'                   ,/,&
'      -----------------------'                   ,/)
 1001 format(/,                                      &
'   ** SOLVING Rij-EPSILON SSG'                   ,/,&
'      -----------------------'                   ,/)
 1002 format(/,                                      &
'   ** SOLVING Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------'                 ,/)

#endif

!----
! End
!----

return

end subroutine
