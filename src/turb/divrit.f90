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

!> \file divrit.f90
!>
!> \brief This subroutine perform  add the divergence of turbulent flux
!> to the transport equation of a scalar.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     coefa         boundary condition array for the variable
!> \param[in]     coefb         boundary condition array for the variable
!> \param[in]     xcpp          Cp
!> \param[out]    smbrs         Right hand side to update
!_______________________________________________________________________________

subroutine divrit &
 ( nvar   , nscal  ,                                              &
   iscal  , itspdv ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  ,                                              &
   xcpp   ,                                                       &
   smbrs )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum!
use pointe
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iscal  , itspdv
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision xcpp(ncelet)
double precision smbrs(ncelet)

! Local variables

integer          ifac, init, inc
integer          iccocg,iflmb0,imaspe
integer          ipcrom, ipbrom
integer          nswrgp, imligp, iwarnp
integer          itypfl
integer          ivar , iclvar, iel, ii, jj, isou
integer          itt
integer          f_id

double precision epsrgp, climgp, extrap
double precision xk, xe, xtt
double precision grav(3),xrij(3,3), temp(3)

logical          ilved

character*80     fname

double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:,:), pointer :: coefav
double precision, dimension(:,:,:), pointer :: coefbv
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: gradt
double precision, allocatable, dimension(:,:) :: coefat
double precision, allocatable, dimension(:,:,:) :: coefbt
double precision, allocatable, dimension(:) :: thflxf, thflxb
double precision, allocatable, dimension(:) :: divut
double precision, allocatable, dimension(:,:) :: w1

double precision, dimension(:,:), pointer :: cofarut
double precision, dimension(:,:,:), pointer :: cofbrut
double precision, dimension(:,:), pointer :: xut
double precision, dimension(:,:), pointer :: xuta

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Initializations to avoid compiler warnings
xtt = 0.d0

! First component is for x,y,z  and the 2nd for u,v,w
allocate(gradv(ncelet,3,3))
allocate(gradt(ncelet,3), thflxf(nfac), thflxb(nfabor))
ipcrom = ipproc(irom)
ipbrom = ipprob(irom)

! Compute scalar gradient
ivar = isca(iscal)
iccocg = 1
inc = 1

! Name of the scalar ivar
call field_get_name(ivarfl(ivar), fname)

! Index of the corresponding turbulent flux
call field_get_id(trim(fname)//'_turbulent_flux', f_id)

call field_get_val_v(f_id, xut)


nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

! Boundary condition pointers for gradients and advection
call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,         &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
   rtpa(1,ivar)    , coefap , coefbp ,                           &
   gradt  )

! Compute velocity gradient
iccocg = 1
inc    = 1
nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(iu)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)
iclvar = iclrtp(iu,icoef)

! Boundary condition pointers for gradients and advection
call field_get_coefa_v(ivarfl(iu), coefav)
call field_get_coefb_v(ivarfl(iu), coefbv)


ilved = .false.

call grdvec &
!==========
( iu     , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  ilved ,                                                        &
  rtp(1,iu)       , coefav , coefbv ,                            &
  gradv  )

! Find the variance of the thermal scalar
itt = -1
if (((abs(gx)+abs(gy)+abs(gz)).gt.epzero).and.irovar.gt.0.and.         &
    ((ityturt(iscal).eq.2).or.(ityturt(iscal).eq.3))) then
  grav(1) = gx
  grav(2) = gy
  grav(3) = gz
  do ii = 1, nscal
    if (iscavr(ii).eq.iscalt) itt = ii
  enddo
  if (itt.le.0) then
    write(nfecra,9999)
    call csexit(1)
  endif
endif

!===============================================================================
! 2. Agebraic models AFM
!===============================================================================
if (ityturt(iscal).ne.3) then

  allocate(w1(3,ncelet))

  do ifac = 1, nfac
    thflxf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    thflxb(ifac) = 0.d0
  enddo

  do iel = 1, ncel
    !Rij
    xrij(1,1) = rtpa(iel,ir11)
    xrij(2,2) = rtpa(iel,ir22)
    xrij(3,3) = rtpa(iel,ir33)
    xrij(1,2) = rtpa(iel,ir12)
    xrij(1,3) = rtpa(iel,ir13)
    xrij(2,3) = rtpa(iel,ir23)
    xrij(2,1) = xrij(1,2)
    xrij(3,1) = xrij(1,3)
    xrij(3,2) = xrij(2,3)
    ! Epsilon
    xe = rtpa(iel,iep)
    ! Kinetic turbulent energy
    xk = 0.5d0*(xrij(1,1)+xrij(2,2)+xrij(3,3))

    !  Turbulent time-scale (constant in AFM)
    if (iturt(iscal).eq.20) then
      xtt = xk/xe
    else
      xtt = xk/xe
    endif

    ! Compute thermal flux u'T'

    !FIXME compute u'T' for GGDH.
    do ii = 1, 3

      temp(ii) = 0.d0

      ! AFM and EB-AFM models
      !  "-C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)"
      if (ityturt(iscal).eq.2.and.ibeta.gt.0) then
        if (itt.gt.0) then
          temp(ii) = temp(ii) - ctheta(iscal)*xtt*                            &
                       etaafm*propce(iel,ipproc(ibeta))*grav(ii)*rtpa(iel,isca(itt))
        endif

        do jj = 1, 3
          if (ii.ne.jj) then
            temp(ii) = temp(ii)                                               &
                     - ctheta(iscal)*xtt*xiafm*gradv(iel,jj,ii)*xut(jj,iel)
          endif
        enddo
      endif

      ! Partial implicitation of "-C_theta*k/eps*( xi* uT'.Grad u )"
      if (iturt(iscal).eq.20) then
        temp(ii) = temp(ii)/(1.d0+ctheta(iscal)*xtt*xiafm*gradv(iel,ii,ii))
      endif

    enddo

    ! Add the term in "grad T" which is implicited by the GGDH part in covofi.
    !  "-C_theta*k/eps* R.grad T"
    do ii = 1, 3
      xut(ii,iel) = temp(ii) - ctheta(iscal)*xtt*( xrij(ii,1)*gradt(iel,1)  &
                                                 + xrij(ii,2)*gradt(iel,2)  &
                                                 + xrij(ii,3)*gradt(iel,3))
      ! In the next step, we compute the divergence of:
      !  "-Cp*C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)"
      !  The part "-C_theta*k/eps* R.Grad T" is computed by the GGDH part
      w1(ii,iel) = xcpp(iel)*temp(ii)
    enddo
  enddo

  itypfl = 1
  iflmb0 = 1
  init   = 1
  inc    = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  ! Local gradient boundaray conditions: homogenous Neumann
  allocate(coefat(3,ndimfb))
  allocate(coefbt(3,3,ndimfb))
  do ifac = 1, nfabor
    do ii = 1, 3
    coefat(ii,ifac) = 0.d0
      do jj = 1, 3
        if (ii.eq.jj) then
          coefbt(ii,jj,ifac) = 1.d0
        else
          coefbt(ii,jj,ifac) = 0.d0
        endif
      enddo
    enddo
  enddo

  call inimav &
  !==========
  ( ivar   , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp  , imligp, &
    iwarnp , nfecra ,                                     &
    epsrgp , climgp , extrap ,                            &
    propce(1,ipcrom), propfb(1,ipbrom),                   &
    w1     ,                                              &
    coefat , coefbt ,                                     &
    thflxf , thflxb )

  deallocate(coefat)
  deallocate(coefbt)
  deallocate(w1)

!===============================================================================
! 3. Transport equation on turbulent thermal fluxes (DFM)
!===============================================================================
else

  call field_get_val_prev_v(f_id, xuta)

  call resrit &
  !==========
( nscal  ,                                               &
  iscal  , xcpp   , xut    , xuta   ,                    &
  dt     , rtp    , rtpa   , propce ,                    &
  gradv  , gradt  )

  itypfl = 1
  iflmb0 = 1
  init   = 1
  inc    = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  do iel = 1, ncelet
    xuta(1,iel) = xut(1,iel)
    xuta(2,iel) = xut(2,iel)
    xuta(3,iel) = xut(3,iel)
  enddo

  allocate(w1(3, ncelet))

  do iel = 1, ncelet
    w1(1,iel) = xcpp(iel)*xut(1,iel)
    w1(2,iel) = xcpp(iel)*xut(2,iel)
    w1(3,iel) = xcpp(iel)*xut(3,iel)
  enddo

  ! Boundary Conditions on T'u' for the divergence term of
  ! the thermal transport equation
  call field_get_coefad_v(f_id,cofarut)
  call field_get_coefbd_v(f_id,cofbrut)

  call inimav &
  !==========
  ( ivar   , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp  , imligp, &
    iwarnp , nfecra ,                                     &
    epsrgp , climgp , extrap ,                            &
    propce(1,ipcrom), propfb(1,ipbrom),                   &
    w1     ,                                              &
    cofarut, cofbrut,                                     &
    thflxf , thflxb )

  deallocate(w1)

endif

!===============================================================================
! 4. Add the divergence of the thermal flux to the thermal transport equation
!===============================================================================

if ((ityturt(iscal).eq.2.or.ityturt(iscal).eq.3)) then
  allocate(divut(ncelet))

  init = 1

  call divmas &
  !==========
   ( ncelet , ncel   , nfac  , nfabor , init   , nfecra ,          &
     ifacel , ifabor ,                                             &
     thflxf , thflxb , divut )

  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) - divut(iel)
  enddo

  ! Free memory
  deallocate(divut)

endif

! Free memory
deallocate(gradv)
deallocate(gradt)
deallocate(thflxf)
deallocate(thflxb)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 9999 format( &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS' ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute'                             ,/,&
'@'                                                            ,/,&
'@  Le modele de flux thermique turbulent choisi        '      ,/,&
'@  necessite le calcul de la variance du scalaire thermique'  ,/,&
'@'                                                            ,/,&
'@  Verifier les donnees entrees dans l''interface'            ,/,&
'@    et dans les sous-programmes utilisateur.'                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 9999 format( &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,&
'@    ========'                                                ,/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run                  '         ,/,&
'@'                                                            ,/,&
'@  Turbulent heat flux model taken imposed that   '           ,/,&
'@  Thermal scalar variance has to be calculate.   '           ,/,&
'@'                                                            ,/,&
'@  Verify the provided data in the interface'                 ,/,&
'@    and in user subroutines.'                                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

end subroutine
