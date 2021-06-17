!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2021 EDF S.A.
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
! --------
!> \file visdyn.f90
!> \brief Calculation of turbulent viscosity for
!>        a dynamic Smagorinsky LES model
!>
!> \f[ smago = \dfrac{L_{ij}M_{ij}}{M_{ij}M_{ij}} \f]
!>
!> \f[ \mu_T = \rho smago L^2  \sqrt{2 S_{ij}S_{ij}} \f]
!> \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
!>
!> We have at edge faces types at previous time step
!>   (except at first time step, when tables itypfb and itrifb
!>   have not been filled).
!>
!> Please refer to the
!> <a href="../../theory.pdf#dynsmago"><b>dynamic Smagorinsky model</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        number of ncepdp cells with losses
!> \param[in]     icetsm        number of cells with mass source
!> \param[in]     itypsm        type of mass source for the variable
!>                               (cf. cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head losses
!> \param[in]     smacel        value of variables associated to the
!>                               mass source
!>                               for ivar = ipr, smacel = mass flux
!> \param[out]    gradv         the computed velocity gradients
!______________________________________________________________________________!

subroutine visdyn &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel, gradv )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstnum
use optcal
use cstphy
use ppincl
use entsor
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision, intent(inout) :: gradv(3,3,ncelet)

! Local variables

integer          ii, jj, iel, inc, ivar
integer          iprev
integer          iclipc
integer          iccocg, iclip
integer          key_turb_diff, key_sgs_sca_coef
integer          t_dif_id, sca_dync_id
integer          nswrgp, imligp, iwarnp

double precision coef, radeux, deux, delta, deltaf
double precision s11, s22, s33, s11f, s22f, s33f
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision dudyf, dudzf, dvdxf, dvdzf, dwdxf, dwdyf
double precision xfil, xa, xb, xfil2, xsmgmx, xsmgmn
double precision xl11, xl22, xl33, xl12, xl13, xl23
double precision xm11, xm22, xm33, xm12, xm13, xm23
double precision smagma, smagmi, smagmy
double precision scal1, scal2, scal3
double precision epsrgp, climgp, extrap

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: s_n, sf_n
double precision, allocatable, dimension(:) :: w0, xrof, xro
double precision, allocatable, dimension(:,:) :: grads, scami, scamif
double precision, allocatable, dimension(:,:) :: xmij, w61, w62, gradsf
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom, coefas, coefbs, cvar_sca
double precision, dimension(:), pointer :: cpro_turb_diff
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: visct, cpro_smago, cpro_sca_dync

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1.  Iniatilization
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_val_s(ivisct, visct)
call field_get_val_s(icrom, crom)

call field_get_val_s(ismago, cpro_smago)

! For the calculation of the viscosity of the sub-mesh
xfil   = xlesfl
xfil2  = xlesfd
xa     = ales
xb     = bles
deux   = 2.d0
radeux = sqrt(deux)
xsmgmx = smagmx
xsmgmn = smagmn

! Allocate some work arrays

allocate(w0(ncelet), w1(ncelet))
allocate(xmij(6,ncelet))
allocate(xro(ncelet), xrof(ncelet))

! Take into account variable density case: Favre filtering
! Constant density case: Reynolds filtering
if(irovar.eq.1) then
  do iel = 1, ncel
    xro(iel) = crom(iel)
  enddo
else
  xro(:) = 1.d0
endif

! In case of constant density, xrof always 1.d0
call les_filter(1, xro, xrof)
!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(s_n(ncelet), sf_n(ncelet))
allocate(w61(6,ncelet), w62(6,ncelet))

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, 0, inc, gradv)

do iel = 1, ncel

  ! gradv(iel, xyz, uvw)
  s11   = gradv(1, 1, iel)
  s22   = gradv(2, 2, iel)
  s33   = gradv(3, 3, iel)
  dudy  = gradv(2, 1, iel)
  dudz  = gradv(3, 1, iel)
  dvdx  = gradv(1, 2, iel)
  dvdz  = gradv(3, 2, iel)
  dwdx  = gradv(1, 3, iel)
  dwdy  = gradv(2, 3, iel)



! In the case of constant density, s11+s22+s33 is zero
  xmij(1,iel) = s11 - irovar*1.0d0/3.0d0*(s11+s22+s33)
  xmij(2,iel) = s22 - irovar*1.0d0/3.0d0*(s11+s22+s33)
  xmij(3,iel) = s33 - irovar*1.0d0/3.0d0*(s11+s22+s33)
  xmij(4,iel) = 0.5d0*(dudy+dvdx)
  xmij(5,iel) = 0.5d0*(dudz+dwdx)
  xmij(6,iel) = 0.5d0*(dvdz+dwdy)

  s_n(iel) = radeux*sqrt(                                      &
             xmij(1,iel)**2 + xmij(2,iel)**2 + xmij(3,iel)**2  &
           + 2.d0*(xmij(4,iel)**2 + xmij(5,iel)**2             &
           + xmij(6,iel)**2) )
  w62(:,iel) = xro(iel)*xmij(:,iel)
enddo

! w62 temperarily contains rho*S
call les_filter(6, w62, w61)

! w61 <rho*S>/<rho>, sf_n is ||<rho*S>/<rho>||
do iel = 1, ncel
  w61(:,iel) = w61(:,iel)/xrof(iel)
  sf_n(iel) = radeux*sqrt(                                      &
              w61(1,iel)**2 + w61(2,iel)**2 + w61(3,iel)**2     &
            + 2.d0*(w61(4,iel)**2 + w61(5,iel)**2               &
            + w61(6,iel)**2) )
enddo

!     Here XMIJ contains Sij
!         S_n contains ||S||
!            sqrt(2)*sqrt(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
!         Sf_n               contains ||SF||
!            sqrt(2)*sqrt(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2))

!===============================================================================
! 3.  Calculation of Mij
!===============================================================================

do iel = 1, ncel
  w0(iel) = xfil *(xa*volume(iel))**xb
enddo

! Reuse xmij as temporary array

do iel = 1, ncel
  delta = w0(iel)
  do ii = 1, 6
    xmij(ii,iel) = -deux*xro(iel)*delta**2*s_n(iel)*xmij(ii,iel)
  enddo
enddo

! w62 now contains <-2*rho*delta**2*||S||*S>
call les_filter(6, xmij, w62)

! Now compute final xmij value: M_ij = alpha_ij - beta_ij

do iel = 1, ncel
  delta = w0(iel)
  deltaf = xfil2*delta
  do ii = 1, 6
    xmij(ii,iel) = -deux*xrof(iel)*deltaf**2*sf_n(iel)*w61(ii,iel) - w62(ii,iel)
  enddo
enddo

deallocate(w61, w62)

!===============================================================================
! 4.  Calculation of the dynamic Smagorinsky constant
!===============================================================================

! Allocate work arrays
allocate(w2(ncelet), w3(ncelet), w4(ncelet))
allocate(w5(ncelet), w6(ncelet), w7(ncelet))
allocate(w8(ncelet), w9(ncelet))

! Filtering the velocity and its square

! U**2
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(1,iel)*vel(1,iel)
enddo
call les_filter(1, w0, w1)

! V**2
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(2,iel)*vel(2,iel)
enddo
call les_filter(1, w0, w2)

! W**2
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(3,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w3)

! UV
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(1,iel)*vel(2,iel)
enddo
call les_filter(1, w0, w4)

! UW
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(1,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w5)

! VW
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(2,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w6)

! U
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(1,iel)
enddo
call les_filter(1, w0, w7)
do iel = 1, ncel
  w7(iel) = w7(iel)/xrof(iel)
enddo

! V
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(2,iel)
enddo
call les_filter(1, w0, w8)
do iel = 1, ncel
  w8(iel) = w8(iel)/xrof(iel)
enddo

! W
do iel = 1,ncel
  w0(iel) = xro(iel)*vel(3,iel)
enddo
call les_filter(1, w0, w9)
do iel = 1, ncel
  w9(iel) = w9(iel)/xrof(iel)
enddo

do iel = 1, ncel

  ! Calculation of Lij
  xl11 = w1(iel) - xrof(iel) * w7(iel) * w7(iel)
  xl22 = w2(iel) - xrof(iel) * w8(iel) * w8(iel)
  xl33 = w3(iel) - xrof(iel) * w9(iel) * w9(iel)
  xl12 = w4(iel) - xrof(iel) * w7(iel) * w8(iel)
  xl13 = w5(iel) - xrof(iel) * w7(iel) * w9(iel)
  xl23 = w6(iel) - xrof(iel) * w8(iel) * w9(iel)

  xm11 = xmij(1,iel)
  xm22 = xmij(2,iel)
  xm33 = xmij(3,iel)
  xm12 = xmij(4,iel)
  xm13 = xmij(5,iel)
  xm23 = xmij(6,iel)
  ! Calculation of Mij :: Lij
  w1(iel) = xm11 * xl11 + 2.d0* xm12 * xl12 + 2.d0* xm13 * xl13  &
                        +       xm22 * xl22 + 2.d0* xm23 * xl23  &
                        +                           xm33 * xl33
  ! Calculation of Mij :: Mij
  w2(iel) = xm11 * xm11 + 2.d0* xm12 * xm12 + 2.d0* xm13 * xm13  &
                        +       xm22 * xm22 + 2.d0* xm23 * xm23  &
                        +                           xm33 * xm33

enddo

deallocate(xmij)

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w1)
  call synsca(w2)
endif

! By default we make a local average of numerator and of
! denominator, then only we make the quotient.
! The user can do otherwise in ussmag.

call les_filter(1, w1, w3)

call les_filter(1, w2, w4)

do iel = 1, ncel
  if(abs(w4(iel)).le.epzero) then
    cpro_smago(iel) = xsmgmx
  else
    cpro_smago(iel) = w3(iel)/w4(iel)
  endif
enddo

call ussmag                                                       &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   w1     , w2     )

iclipc = 0
do iel = 1, ncel
  if(cpro_smago(iel).ge.xsmgmx) then
    cpro_smago(iel) = xsmgmx
    iclipc = iclipc + 1
  elseif(cpro_smago(iel).le.xsmgmn) then
    cpro_smago(iel) = xsmgmn
    iclipc = iclipc + 1
  endif
enddo

!===============================================================================
! 5.  Calculation of (dynamic) viscosity
!===============================================================================

do iel = 1, ncel
  coef = cpro_smago(iel)
  delta  = xfil * (xa*volume(iel))**xb
  visct(iel) = crom(iel) * coef * delta**2 * s_n(iel)
enddo

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)

!     Some printings
if (vcopt%iwarni.ge.1) then

  smagma = -1.0d12
  smagmi =  1.0d12
  smagmy =  0.d0
  do iel = 1, ncel
    smagma = max(smagma,cpro_smago(iel))
    smagmi = min(smagmi,cpro_smago(iel))
    smagmy = smagmy + cpro_smago(iel)*volume(iel)
  enddo
  if (irangp.ge.0) then
    call parmax(smagma)
    call parmin(smagmi)
    call parsom(smagmy)
    call parcpt(iclipc)
  endif
  smagmy = smagmy / voltot
  write(nfecra,1000) iclipc
  write(nfecra,2001)
  write(nfecra,2002) smagmy, smagmi, smagma
  write(nfecra,2003)

endif

!===============================================================================
! 6.  Scalar turbulent model
!===============================================================================
! In case of gaz combustion, the SGS scalar flux constant and the turbulent
! diffusivity are only evaluated with the mixture fraction, then applied
! automatically to the other scalar equations

call field_get_key_id("turbulent_diffusivity_id", key_turb_diff)
call field_get_key_id("sgs_scalar_flux_coef_id", key_sgs_sca_coef)

do jj = 1, nscal
  ivar = isca(jj)

  ! For variance of a scalar, the turbulent diffsivity must not be computed
  if (iscavr(jj).gt.0) cycle

  ! For any scalar other than the mixture fraction in diffusion flames,
  ! Dt is not computed either. TODO Soot may be an exception
  if (ippmod(icod3p).ge.0) then
    if (jj.ne.ifm)  cycle
  endif

  call field_get_key_int(ivarfl(ivar), key_turb_diff, t_dif_id)
  call field_get_key_int(ivarfl(ivar), key_sgs_sca_coef, sca_dync_id)

  if (t_dif_id.ge.0.and.sca_dync_id.ge.0) then
    call field_get_val_s(t_dif_id, cpro_turb_diff)
    call field_get_val_s(sca_dync_id, cpro_sca_dync)

    !================================================================
    ! 6.1.  Compute the Mi for scalar
    !================================================================

    allocate(grads(3, ncelet))
    allocate(gradsf(3, ncelet))

    inc = 1
    iprev = 0
    call field_get_coefa_s(ivarfl(ivar), coefas)
    call field_get_coefb_s(ivarfl(ivar), coefbs)

    call field_get_val_s(ivarfl(ivar),cvar_sca)
    call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc, 0, grads)

    ! compute grad (<rho.Y>/<rho>)
    allocate(scami(3,ncelet))
    allocate(scamif(3,ncelet))

    do iel = 1, ncel
      w0(iel) = cvar_sca(iel)*xro(iel)
    enddo
    call les_filter(1, w0, w4)
    do iel = 1, ncel
      w4(iel) = w4(iel)/xrof(iel)
    enddo

    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
    inc = 1
    nswrgp = vcopt%nswrgr
    epsrgp = vcopt%epsrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    climgp = vcopt%climgr
    extrap = vcopt%extrag
    iccocg = 1

    call gradient_s                                                   &
      (-1     , imrgra , inc    , iccocg , nswrgp , imligp ,          &
      iwarnp  , epsrgp , climgp ,                                     &
      w4      , coefas , coefbs ,                                     &
      gradsf   )

    do iel = 1, ncel
      delta  = xfil * (xa*volume(iel))**xb
      do ii = 1, 3
        scami(ii,iel) = -xro(iel)*delta**2*s_n(iel)*grads(ii,iel)
      enddo
    enddo

    call les_filter(3, scami, scamif)
    do iel = 1, ncel
      deltaf = xfil2*xfil * (xa*volume(iel))**xb
      do ii = 1, 3
        scami(ii,iel) = -deltaf**2*xrof(iel)*sf_n(iel)*gradsf(ii,iel)  &
                      - scaMif(ii,iel)
      enddo
    enddo

    !================================================================
    ! 6.2.  Compute the Li for scalar
    !================================================================
    ! rho*U*Y
    do iel = 1, ncel
      w0(iel) = xro(iel)*vel(1,iel)*cvar_sca(iel)
    enddo
    call les_filter(1, w0, w1)

    ! rho*V*Y
    do iel = 1, ncel
      w0(iel) = xro(iel)*vel(2,iel)*cvar_sca(iel)
    enddo
    call les_filter(1, w0, w2)

    ! rho*W*Y
    do iel = 1, ncel
      w0(iel) = xro(iel)*vel(3,iel)*cvar_sca(iel)
    enddo
    call les_filter(1, w0, w3)

    do iel = 1, ncel
      scal1 = w1(iel) - xrof(iel)*w7(iel)*w4(iel)
      scal2 = w2(iel) - xrof(iel)*w8(iel)*w4(iel)
      scal3 = w3(iel) - xrof(iel)*w9(iel)*w4(iel)

      w1(iel) = scal1*scami(1,iel) + scal2*scami(2,iel) + scal3*scami(3,iel)
      w2(iel) = scami(1,iel)**2 + scami(2,iel)**2 + scami(3,iel)**2
    enddo

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(w1)
      call synsca(w2)
    endif

    call les_filter(1, w1, w3)
    call les_filter(1, w2, w4)

    !================================================================
    ! 6.3.  Compute the SGS flux coefficient and SGS diffusivity
    !       Cs >= 0, Dt >=0
    !================================================================

    do iel = 1, ncel
      if(abs(w4(iel)).le.epzero) then
        cpro_sca_dync(iel) = 0.d0
      else
        cpro_sca_dync(iel) = max(w3(iel)/w4(iel), 0.d0)
      endif

      delta  = xfil * (xa*volume(iel))**xb
      cpro_turb_diff(iel) = crom(iel) * cpro_sca_dync(iel) * delta**2 * s_n(iel)
    enddo

    deallocate(scami, scamif)
    deallocate(grads, gradsf)
  endif

enddo


! Free memory
deallocate(s_n, sf_n)
deallocate(w9, w8)
deallocate(w7, w6, w5)
deallocate(w4, w3, w2)
deallocate(w1, w0)
deallocate(xro, xrof)

!----
! Formats
!----

 1000 format(                                                           &
' Nb of clipping of the Smagorinsky constant by max values',I10,/)
 2001 format(                                                           &
' --- Informations on the squared Smagorinsky constant'        ,/,&
' --------------------------------'                            ,/,&
' Mean value  Min value  Max value'                            ,/,&
' --------------------------------'                              )
 2002 format(                                                           &
 e12.4    ,      e12.4,      e12.4                               )
 2003 format(                                                           &
' --------------------------------'                            ,/)

!----
! End
!----

return
end subroutine
