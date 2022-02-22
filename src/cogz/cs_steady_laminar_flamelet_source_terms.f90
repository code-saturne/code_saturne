!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_steady_laminar_flamelet_source_terms.f90
!>
!> \brief Specific physic subroutine: STE/VTE and progress variable equations.
!>
!> This subroutine defines the source terms for the soot mass fraction
!> and the precursor number for soot model of Moss et al for one time step.
!
!  The equations read: \f$ rovsdt \delta a = smbrs \f$
!
!  \f$ rovsdt \f$ et \f$ smbrs \f$ could already contain source term
!  and don't have to be erased but incremented.
!
!  For stability sake, only positive terms should be add in \f$ rovsdt \f$.
!  There is no constrain for \f$ smbrs \f$.
!
!  For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
!           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
!           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
!
!  Here are set \f$ rovsdt \f$ and \f$ smbrs \f$ containning \f$ \rho \Omega \f$
!   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
!     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.°C.s^{-1} \f$,
!     for enthalpy: \f$ J.s^{-1} \f$)
!   - \f$ rovsdt \f$ en \f$ kg.s^{-1} \f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar index
!> \param[in,out] smbrs         explicit right hand side
!> \param[in,out] rovsdt        implicit terms
!_______________________________________________________________________________

subroutine cs_steady_laminar_flamelet_source_terms                &
 ( iscal  ,                                                       &
   smbrs  , rovsdt )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use pointe, only: itypfb
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character(len=80) :: chaine
integer          ivar, iel, ifac, t_dif_id, key_turb_diff, ifcvsl
integer          imrgrp, nswrgp, imligp, iwarnp
integer          iprev, inc, iccocg

double precision  epsrgp, climgp
double precision  cexp, cimp, delta_les

double precision, dimension(:), pointer :: crom, fp2m
double precision, dimension(:), pointer :: cvara_scal, cvara_fm
double precision, dimension(:), pointer :: cpro_omegac, coefap, coefbp
double precision, dimension(:), pointer :: cpro_viscls, cpro_turb_diff

double precision, dimension(:,:), allocatable :: grad
double precision, dimension(:), allocatable :: coefa_p, coefb_p


type(var_cal_opt) :: vcopt_fm, vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

ivar = isca(iscal)
call field_get_label(ivarfl(ivar), chaine)
call field_get_val_s(icrom, crom)

call field_get_val_prev_s(ivarfl(ivar), cvara_scal)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

!===============================================================================
! 2. Writings
!===============================================================================

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

if(mode_fp2m.eq.0) then
  call field_get_val_prev_s(ivarfl(isca(ifp2m)), fp2m)
elseif(mode_fp2m.eq.1) then
  call field_get_val_s(irecvr, fp2m)
endif

if (iturb.eq.41) then
  ! Retrieve turbulent diffusivity value for the mixture fraction
  call field_get_key_id("turbulent_diffusivity_id", key_turb_diff)
  call field_get_key_int(ivarfl(isca(ifm)), key_turb_diff, t_dif_id)
  if (t_dif_id.ge.0) then
    call field_get_val_s(t_dif_id, cpro_turb_diff)
  endif

endif

!=======================================================================
! --- Cuenot et al.:
! STE: Prod := 0
!      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
! VTE: Prod := 2*rho*(D + Dtur)*|grad(Z)|**2
!      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
!
! --- Pierce:
! Progress variable equation:
!      Prod := flamelet_lib(fm, fp2m, ki, progvar)
!=======================================================================

cexp = 0.d0; cimp = 0.d0

! For the moment, this model for source computation is only available in LES
if (iturb.eq.41) then
  if (mode_fp2m .eq. 0) then
    if (ivar.eq.isca(ifp2m)) then
      ! Allocate a temporary array for the gradient reconstruction
      allocate(grad(3,ncelet))
      allocate(coefa_p(nfabor), coefb_p(nfabor))

      iprev = 1
      inc = 1
      iccocg = 1

      ! Homogeneous Neumann on convective inlet on the production term for the
      ! variance
      call field_get_val_prev_s(ivarfl(isca(ifm)), cvara_fm)
      call field_get_coefa_s (ivarfl(isca(ifm)), coefap)
      call field_get_coefb_s (ivarfl(isca(ifm)), coefbp)

      ! pas de diffusion en entree
      do ifac = 1, nfabor
        coefa_p(ifac) = coefap(ifac)
        coefb_p(ifac) = coefbp(ifac)
        if (itypfb(ifac).eq.i_convective_inlet) then
          coefa_p(ifac) = 0.d0
          coefb_p(ifac) = 1.d0
        endif
      enddo

      call field_get_key_struct_var_cal_opt(ivarfl(isca(ifm)), vcopt_fm)

      imrgrp = vcopt_fm%imrgra
      nswrgp = vcopt_fm%nswrgr
      imligp = vcopt_fm%imligr
      iwarnp = vcopt_fm%iwarni
      epsrgp = vcopt_fm%epsrgr
      climgp = vcopt_fm%climgr

      call gradient_s                                                         &
        ( ivarfl(isca(ifm))  , imrgrp , inc    , iccocg , nswrgp , imligp ,   &
        iwarnp          , epsrgp , climgp ,                                   &
        cvara_fm        , coefa_p, coefb_p,                                   &
        grad )

      deallocate (coefa_p, coefb_p)

      do iel = 1, ncel
        delta_les = xlesfl *(ales*volume(iel))**bles
        cexp = + 2.d0*(cpro_turb_diff(iel) + cpro_viscls(iel))*cell_f_vol(iel) &
          *(grad(1,iel)**2.d0 + grad(2,iel)**2.d0 + grad(3,iel)**2.d0)         &
          - ((cpro_turb_diff(iel) + cpro_viscls(iel))/(coef_k                  &
          * delta_les**2.d0) * fp2m(iel))*cell_f_vol(iel)

        cimp = 0.d0
        smbrs(iel)  = smbrs(iel)  + cexp + cimp*cvara_scal(iel)
        rovsdt(iel) = rovsdt(iel) + max(-cimp,0.d0)
      enddo
      deallocate(grad)

    endif


  elseif (mode_fp2m .eq. 1) then
    if (ivar .eq. isca(ifsqm)) then
      iprev = 1
      call combustion_reconstruct_variance(iprev)

      do iel = 1, ncel
        delta_les = xlesfl *(ales*volume(iel))**bles
        cexp =  - ((cpro_turb_diff(iel) + cpro_viscls(iel))/(coef_k       &
                * delta_les**2.d0)*fp2m(iel))*cell_f_vol(iel)

        cimp = 0.d0
        smbrs(iel)  = smbrs(iel)  + cexp + cimp*cvara_scal(iel)
        rovsdt(iel) = rovsdt(iel) + max(-cimp,0.d0)
      enddo
    endif
  endif
endif

if (ippmod(islfm).ge.2) then
  if (ivar.eq.isca(ipvm)) then

    call field_get_val_s(iomgc, cpro_omegac)
    do iel = 1, ncel

      cexp = cpro_omegac(iel)

      cimp = 0.d0

      smbrs(iel)  = smbrs(iel)  + cexp + cimp*cvara_scal(iel)
      rovsdt(iel) = rovsdt(iel) + max(-cimp,0.d0)

    enddo

  endif
endif


!--------
! Formats
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! End
!----

return

end subroutine
