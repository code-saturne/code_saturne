!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!______________________________________________________________________________!

subroutine visdyn &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstnum
use optcal
use cstphy
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
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          ii, iel, inc
integer          iprev
integer          iclipc

double precision coef, radeux, deux, delta, deltaf
double precision s11, s22, s33, s11f, s22f, s33f
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision dudyf, dudzf, dvdxf, dvdzf, dwdxf, dwdyf
double precision xfil, xa, xb, xfil2, xsmgmx, xsmgmn
double precision xl11, xl22, xl33, xl12, xl13, xl23
double precision xm11, xm22, xm33, xm12, xm13, xm23
double precision smagma, smagmi, smagmy

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10, w0
double precision, allocatable, dimension(:,:) :: xmij, w61, w62
double precision, dimension(:,:,:), allocatable :: gradv, gradvf
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: visct, cpro_smago

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

!===============================================================================
! 2.  Calculation of velocity gradient and of
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet), gradvf(3,3,ncelet))

inc = 1
iprev = 0

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc, gradv)

! Filter the velocity gradient on the extended neighborhood

call les_filter(9, gradv, gradvf)

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

  s11f  = gradvf(1, 1, iel)
  s22f  = gradvf(2, 2, iel)
  s33f  = gradvf(3, 3, iel)
  dudyf = gradvf(2, 1, iel)
  dudzf = gradvf(3, 1, iel)
  dvdxf = gradvf(1, 2, iel)
  dvdzf = gradvf(3, 2, iel)
  dwdxf = gradvf(1, 3, iel)
  dwdyf = gradvf(2, 3, iel)

  xmij(1,iel) = s11
  xmij(2,iel) = s22
  xmij(3,iel) = s33
  xmij(4,iel) = 0.5d0*(dudy+dvdx)
  xmij(5,iel) = 0.5d0*(dudz+dwdx)
  xmij(6,iel) = 0.5d0*(dvdz+dwdy)

  visct(iel) = radeux*sqrt(                                       &
                       s11**2 + s22**2 + s33**2                   &
                     + 0.5d0*( (dudy+dvdx)**2                     &
                             + (dudz+dwdx)**2                     &
                             + (dvdz+dwdy)**2 )  )

  w1(iel) = radeux*sqrt(                                          &
                       s11f**2 + s22f**2 + s33f**2                &
                     + 0.5d0*( (dudyf+dvdxf)**2                   &
                             + (dudzf+dwdxf)**2                   &
                             + (dvdzf+dwdyf)**2 )  )
enddo

! Free memory
deallocate(gradv, gradvf)

!     Here XMIJ contains Sij
!         VISCT contains ||S||
!            sqrt(2)*sqrt(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
!         W1                 contains ||SF||
!            sqrt(2)*sqrt(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2))

!===============================================================================
! 3.  Calculation of Mij
!===============================================================================

do iel = 1, ncel
  w0(iel) = xfil *(xa*volume(iel))**xb
enddo

allocate(w61(6,ncelet), w62(6,ncelet))

call les_filter(6, xmij, w61)

! Reuse xmij as temporary array

do iel = 1, ncel
  delta = w0(iel)
  do ii = 1, 6
    xmij(ii,iel) = -deux*delta**2*visct(iel)*xmij(ii,iel)
  enddo
enddo

call les_filter(6, xmij, w62)

! Now compute final xmij value: M_ij = alpha_ij - beta_ij

do iel = 1, ncel
  delta = w0(iel)
  deltaf = xfil2*delta
  do ii = 1, 6
    xmij(ii,iel) = -deux*deltaf**2*w1(iel)*w61(ii,iel) - w62(ii,iel)
  enddo
enddo

deallocate(w61, w62)

!===============================================================================
! 4.  Calculation of the dynamic Smagorinsky constant
!===============================================================================

! Allocate work arrays
allocate(w2(ncelet), w3(ncelet), w4(ncelet))
allocate(w5(ncelet), w6(ncelet), w7(ncelet))
allocate(w8(ncelet), w9(ncelet), w10(ncelet))

! Filtering the velocity and its square

! U**2
do iel = 1,ncel
  w0(iel) = vel(1,iel)*vel(1,iel)
enddo
call les_filter(1, w0, w1)

! V**2
do iel = 1,ncel
  w0(iel) = vel(2,iel)*vel(2,iel)
enddo
call les_filter(1, w0, w2)

! W**2
do iel = 1,ncel
  w0(iel) = vel(3,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w3)

! UV
do iel = 1,ncel
  w0(iel) = vel(1,iel)*vel(2,iel)
enddo
call les_filter(1, w0, w4)

! UW
do iel = 1,ncel
  w0(iel) = vel(1,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w5)

! VW
do iel = 1,ncel
  w0(iel) = vel(2,iel)*vel(3,iel)
enddo
call les_filter(1, w0, w6)

! U
do iel = 1,ncel
  w0(iel) = vel(1,iel)
enddo
call les_filter(1, w0, w7)

! V
do iel = 1,ncel
  w0(iel) = vel(2,iel)
enddo
call les_filter(1, w0, w8)

! W
do iel = 1,ncel
  w0(iel) = vel(3,iel)
enddo
call les_filter(1, w0, w9)

do iel = 1, ncel

  ! Calculation of Lij
  xl11 = w1(iel) - w7(iel) * w7(iel)
  xl22 = w2(iel) - w8(iel) * w8(iel)
  xl33 = w3(iel) - w9(iel) * w9(iel)
  xl12 = w4(iel) - w7(iel) * w8(iel)
  xl13 = w5(iel) - w7(iel) * w9(iel)
  xl23 = w6(iel) - w8(iel) * w9(iel)

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
! 3.  Calculation of (dynamic) velocity
!===============================================================================

! Clipping in (mu + mu_t)>0 in phyvar

do iel = 1, ncel
  coef = cpro_smago(iel)
  delta  = xfil * (xa*volume(iel))**xb
  visct(iel) = crom(iel) * coef * delta**2 * visct(iel)
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

! Free memory
deallocate(w10, w9, w8)
deallocate(w7, w6, w5)
deallocate(w4, w3, w2)
deallocate(w1, w0)

!----
! Formats
!----

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
' Nb Clipping Constante  Smagorinsky par valeurs maximales',I10,/)
 2001 format(                                                           &
' --- Informations sur la constante de Smagorinsky^2             ',/,&
' ----------------------------------                          ',/,&
' Val. moy.      Val. min   Val. max                               ',/,&
' ----------------------------------                          '  )
 2002 format(                                                           &
 e12.4    ,      e12.4,      e12.4                               )
 2003 format(                                                           &
' ----------------------------------                          ',/)

#else

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

#endif

!----
! End
!----

return
end subroutine
