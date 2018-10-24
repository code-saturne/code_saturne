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

!> \file visecv.f90
!> \brief Computes the secondary viscosity contribution \f$\kappa
!> -\dfrac{2}{3} \mu\f$ in order to compute:
!> \f[
!> \grad\left( (\kappa -\dfrac{2}{3} \mu) \trace( \gradt(\vect{u})) \right)
!> \f]
!> with:
!>   - \f$ \mu = \mu_{laminar} + \mu_{turbulent} \f$
!>   - \f$ \kappa \f$ is the volume viscosity (generally zero)
!>
!> \remark
!> In LES, the tensor
!> \f$\overline{\left(\vect{u}-\overline{\vect{u}}\right)\otimes\left(\vect{u}
!>-\overline{\vect{u}}\right)}\f$
!> is modeled by \f$\mu_t \overline{\tens{S}}\f$
!> and not by
!> \f$\mu_t\overline{\tens{S}}-\dfrac{2}{3}\mu_t
!> \trace\left(\overline{\tens{S}}\right)\tens{1}+\dfrac{2}{3}k\tens{1}\f$
!> so that no term
!> \f$\mu_t \dive \left(\overline{\vect{u}}\right)\f$ is needed.
!>
!> Please refer to the
!> <a href="../../theory.pdf#visecv"><b>visecv</b></a> section
!> of the theory guide for more informations.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in,out] secvif        lambda*surface at interior faces
!> \param[in,out] secvib        lambda*surface at boundary faces
!______________________________________________________________________________

subroutine visecv &
 ( secvif , secvib )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use entsor
use field
use numvar
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

double precision secvif(nfac), secvib(nfabor)

! Local variables

integer          iel, ifac, ii, jj
integer          key_t_ext_id
integer          iviext

double precision d2s3m, secvsi, secvsj, pnd

double precision, allocatable, dimension(:) :: secvis
double precision, dimension(:), pointer :: porosi
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: cpro_viscv
double precision, dimension(:), pointer :: cproa_viscl, cproa_visct

!===============================================================================

!===============================================================================
! 1.  INITIALIZATION
!===============================================================================

! Allocate temporary arrays
allocate(secvis(ncelet))

call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

if (ippmod(icompf).ge.0) then
  if (iviscv.ge.0) then
    call field_get_val_s(iviscv, cpro_viscv)
  endif
endif

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

!===============================================================================
! 2. Computation of the second viscosity: lambda = K -2/3 mu
!===============================================================================

!  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...

d2s3m = -2.d0/3.d0

! Laminar viscosity
call field_get_key_int(iviscl, key_t_ext_id, iviext)
if(isno2t.gt.0 .and. iviext.gt.0) then
  call field_get_val_prev_s(iviscl, cproa_viscl)
  do iel = 1, ncel
    secvis(iel) = d2s3m*cproa_viscl(iel)
  enddo
else
  do iel = 1, ncel
    secvis(iel) = d2s3m*viscl(iel)
  enddo
endif

! Volume viscosity if present
if (ippmod(icompf).ge.0) then
  if (iviscv.ge.0) then
    do iel = 1, ncel
      secvis(iel) = secvis(iel) + cpro_viscv(iel)
    enddo
  else
    do iel = 1, ncel
      secvis(iel) = secvis(iel) + viscv0
    enddo
  endif
endif

! Turbulent viscosity (if not in Rij or LES)
call field_get_key_int(ivisct, key_t_ext_id, iviext)
if (itytur.ne.3 .and. itytur.ne.4) then
  if(isno2t.gt.0 .and. iviext.gt.0) then
    call field_get_val_prev_s(ivisct, cproa_visct)
    do iel = 1, ncel
      secvis(iel) = secvis(iel) + d2s3m*cproa_visct(iel)
    enddo
  else
    do iel = 1, ncel
      secvis(iel) = secvis(iel) + d2s3m*visct(iel)
    enddo
  endif
endif

! With porosity
if (iporos.eq.1.or.iporos.eq.2) then
  call field_get_val_s(ipori, porosi)
  do iel = 1, ncel
    secvis(iel) = secvis(iel)*porosi(iel)
  enddo
endif


! ---> Parallelism and periodicity treatmenT

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(secvis)
endif

! --- Interior faces
! TODO we should (re)test the weigthen walue for imvisf=0

if (imvisf.eq.0) then
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    secvsi = secvis(ii)
    secvsj = secvis(jj)
    pnd = pond(ifac)

    secvif(ifac) = 0.5d0*(secvsi+secvsj)
  enddo
else
  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    secvsi = secvis(ii)
    secvsj = secvis(jj)
    pnd = pond(ifac)

    secvif(ifac) = secvsi*secvsj/(pnd*secvsi+(1.d0-pnd)*secvsj)
  enddo
endif

! --- Boundary faces
! TODO shall we extrapolate this value?

do ifac = 1, nfabor
  ii = ifabor(ifac)
  secvib(ifac) = secvis(ii)
enddo

! --- TODO stresses at the wall?

! Free memory
deallocate(secvis)

return
end subroutine
