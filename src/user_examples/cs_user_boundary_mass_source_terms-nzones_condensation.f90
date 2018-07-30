!-------------------------------------------------------------------------------

!VERS

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
! ---------

!> \file cs_user_boundary_mass_source_terms-nzones_condensation.f90
!>
!> \brief Source terms for each zone associated at the boundary faces and the
!> neighboring cells with surface condensation.
!>
!> \par Examples of settings for boundary condensation mass source terms
!>      Examples are available
!>      \ref condens_h_boundary "here".
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iappel        indicates which at which stage the routine is
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     itypcd        type of condensation source term for each ivar
!> \param[in]     izftcd        faces zone with condensation source terms imposed
!>                              (at previous and current time steps)
!> \param[out]    spcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, spcond is the flow rate
!>                              \f$ \Gamma_{cond}^n \f$)
!_______________________________________________________________________________

subroutine cs_user_boundary_mass_source_terms &
 ( nvar   , nscal  ,                                              &
   nfbpcd , iappel ,                                              &
   ifbpcd , itypcd , izftcd ,                                     &
   spcond , tpar)

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
use ppincl
use mesh
use field
use cs_c_bindings
use cs_f_interfaces
use cs_nz_condensation, only:nzones,izzftcd,izcophc,izcophg,iztag1d, ztpar
use cs_nz_tagmr

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iappel
integer          nfbpcd

integer          ifbpcd(nfbpcd), itypcd(nfbpcd,nvar)
integer          izftcd(ncel)

double precision spcond(nfbpcd,nvar)
double precision tpar

! Local variables

!< [loc_var_dec]
integer          ieltcd, ii, iz
integer          ifac, iel, iesp, iscal
integer          ivarh
integer          ilelt, nlelt
integer          izone
integer          f_id

double precision hvap, tk

type(gas_mix_species_prop) s_h2o_g

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: cvar_h
!< [loc_var_dec]

!===============================================================================

!< [init]
! Allocate a temporary array for cells selection
allocate(lstelt(nfabor))

call field_get_id_try("y_h2o_g", f_id)
if (f_id.ne.-1) &
  call field_get_key_struct_gas_mix_species_prop(f_id, s_h2o_g)
!< [init]

!< [zones_definition]
if (iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! 1. One or two calls
! -------------------
!  - iappel = 1: nfbpcd: calculation of the number of faces with
!                             condensation source term
!  - iappel = 2: ifbpcd: index number of faces with condensation source terms
!
! Remarks
! =======
!  - Do not use spcond in this section (it is set on the third call, iappel=3)
!  - Do not use ifbpcd in this section on the first call (iappel=1)
!  - This section (iappel=1 or 2) is only accessed at the beginning of a
!     calculation. Should the localization of the condensation source terms evolve
!     in time, the user must identify at the beginning all cells that can
!     potentially become condensation  source term.
!===============================================================================


  !--To Select the cells with condensation source term
  !---------------------------------------------------

  izone = 0
  ieltcd = 0

  ! Cells with a boundary face of color 60
  call getfbr('60 and box[-0.03,0.0,0.5,0.03,0.6,1.22]',nlelt,lstelt)

  izone = izone + 1

  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)
    ieltcd = ieltcd + 1
    if (iappel.eq.2) then
      ifbpcd(ieltcd) = ifac
      izzftcd(ieltcd) = izone
    endif
  enddo

  if (irangp.le.0.and.iappel.eq.2.and.ieltcd.gt.0) then
    write(*,*) "izzftcd(",izone,")= ", izzftcd(ieltcd)
  endif

  !========================================================
  ! Faces selection associated to the Concrete wall surface
  !========================================================
  call getfbr('60 and box[-0.03,0.0,1.22,0.03,0.6,2.55]', nlelt, lstelt)

  izone = izone + 1

  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)
    ieltcd = ieltcd + 1
    if (iappel.eq.2) then
      ifbpcd(ieltcd) = ifac
      izzftcd(ieltcd) = izone
    endif
  enddo

  if (irangp.le.0.and.iappel.eq.2.and.ieltcd.gt.0) then
    write(*,*) "izzftcd(",izone,")= ", izzftcd(ieltcd)
  endif

endif

! For iappel = 1,
! Specification of nfbpcd.
if (iappel.eq.1) then
  nfbpcd = ieltcd
  nzones = izone
endif
!< [zones_definition]

!< [model_settings]
!===============================================================================
! Parameters padding of the wall thermal model and condensation model
! ------------------------------------------------------------------
! The condensation model can be used alone with a constant temperature
! specified by the user at the cold wall (at iappel=3 tpar=tpar0 in this case)
! or together with a 0-D thermal model. In the latter case, the two models are
! coupled.
!===============================================================================

if (iappel.eq.2) then

  do ii = 1, nfbpcd
    iz = izzftcd(ii)
    if (iz.eq.1.and.icondb.eq.0) then
      ! Turbulent law and empiric correlations used to
      ! define the exchange coefficients of the sink
      ! source term and heat transfer to the cooling
      ! wall associated to the condensation phenomenon
      ! given by the COPAIN model.
      !-----------------------------------------------

      ! Choice the way to compute the exchange coefficient (hcond)
      ! associated to the condensation sink source term.

      ! With the parameter icophc defined below:
      ! ----------------------------------------
      !  1 : this one provided by the turbulent flow
      !        at the cooling wall (hcond = hcdt)
      !  2 : this one is given by the copain
      !        correlation (hcond = hcdcop)
      !  3 : this one is obtained by the estimation of
      !        the maximal value between those two previous
      !        coefficient as below : hcond= max(hcdt, hcdcop)
      izcophc(iz) = 3

      ! Choice the way to compute the thermal exchange coefficient
      ! associated to the heat transfer at the cooling wall,
      ! due to the energy loss by condensation phenomenon.

      ! With the parameter icophg defined below:
      ! ----------------------------------------
      !  2 : this one is given by the copain
      !      correlation (hpcond = hw_cop)
      !  3 : this one is obtained by the estimation of
      !      the maximal value between the current
      !      and previous value of hcdcop given by
      !      the copain correlation as below:
      !         hpcond= max(hw_cop^n, hw_cop^n+1)
      izcophg(iz) = 3

      ! Choice the way to impose the wall temperature (tpar)
      ! at the solid/fluid interface:
      !
      ! with the parameter itag1d defined below:
      ! ----------------------------------------
      !  0 : A constant wall temperature imposed is given by the user
      !     ( tpar = tpar0 used as the wall temperature by the condensation model )
      !  1 : A variable wall temperature is imposed with a 1-D thermal model
      !     ( tpar = tmur(ii,1) computed by tagmro.f90 and used as the
      !       wall temperature by the condensation model )
      iztag1d(iz) = 1

      ! Wall temperature computed by a 1-D thermal model
      ! with a implicit scheme and variable over time.
      ! ------------------------------------------------
      ! Remark : the wall temperature is in unit [°C].
      if (iztag1d(iz).eq.1) then
        !---------------------------------------------------
        ! Numerical parameters used by the 1-D thermal model
        !---------------------------------------------------
        ! (theta) parameter of the implicit scheme
        ztheta(iz) = 1.d0
        ! (dxmin) First cell size of the 1D mesh
        ! -> with (dxmin.eq.0) the Dx is constant
        ! -> with (dxmin.ne.0) the Dx is variable
        zdxmin(iz) = 0.d0
        ! (nmur) space steps number of the 1D mesh
        znmur(iz) = 10
        ! (epais) thickness of the 1D wall
        zepais(iz) = 0.024d0

        !-------------------------------------------
        !Initial condition of the 1-D thermal model
        !-------------------------------------------
        ztpar0(iz) = 26.57d0
      endif
    endif

  enddo

  if (irangp.ge.0) then
    call parimx(nzones, izcophc)
    call parimx(nzones, izcophg)
    call parimx(nzones, iztag1d)
    call parrmx(nzones, ztheta)
    call parrmx(nzones, zdxmin)
    call parimx(nzones, znmur )
    call parrmx(nzones, zepais)
    call parrmx(nzones, ztpar0)
  endif

!< [model_settings]

!< [source_terms_values]
elseif (iappel.eq.3) then

!===============================================================================
! 2. For nfbpcd > 0 , third call
!    iappel = 3 : itypcd: type of condensation source term
!                  spcond: condensation source term
! Remark
! ======
! If itypcd(ieltcd,ivar) is set to 1, spcond(ieltcd,ivar) must be set.
!===============================================================================

  !-- pointer to the specific heat
  if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

  !-- pointer to the enthalpy value
  ivarh = isca(iscalt)
  call field_get_val_s(ivarfl(ivarh), cvar_h)

  do ii = 1, nfbpcd

    ifac = ifbpcd(ii)
    iel  = ifabor(ifac)

    iz = izzftcd(ii)

    if (icondb.eq.0) then
      if (iztag1d(iz).eq.1) then
        !-------------------------------------------
        !Boundary conditions of the 1-D thermal model
        !-------------------------------------------
        zhext(iz) =  1.d+8 ; ztext(iz) = 26.57d0
        ! --------------------------------------------
        ! Physical properties of the concrete material
        ! --------------------------------------------
        ! (rob) density (kg.m-3)
        zrob(iz)   = 8000.d0
        ! (condb) thermal conductivity (W.m-1.C-1)
        zcondb(iz) = 12.8d0
        ! (cpb)   Specific heat (J.kg-1.C-1)
        zcpb(iz) = 500.0d0
      else
        ! Wall temperature imposed as constant
        ! with a value specified by the user
        ztpar(iz) = 26.57d0
      endif
    elseif (iz.eq.2.and.icondb.eq.0) then
      if (iztag1d(iz).eq.1) then
        !-------------------------------------------
        !Boundary conditions of the 1-D thermal model
        !-------------------------------------------
        zhext(iz) =  1.d+8 ; ztext(iz) = 26.57d0
        ! --------------------------------------------
        ! Physical properties of the concrete material
        ! --------------------------------------------
        ! (rob) density (kg.m-3)
        zrob(iz)   = 8000.d0
        ! (condb) thermal conductivity (W.m-1.C-1)
        zcondb(iz) = 12.8d0
        ! (cpb)   Specific heat (J.kg-1.C-1)
        zcpb(iz) = 500.0d0
      else
        ! Wall temperature imposed as constant
        ! with a value specified by the user
        ztpar(iz) = 26.57d0
      endif
    endif

    ! To fill the spcond(nfbpcd,ivar) array
    ! if we want to specify a variable value
    !---------------------------------------

    ! Compute the enthalpy value of vapor gas
    if (ntcabs.le.1) then
      tk = t0
    else
      tk = cvar_h(iel)/cpro_cp(iel)
    endif
    hvap = s_h2o_g%cp*tk

    ! any condensation source term
    ! associated to each velocity component
    ! momentum equation in this case.
    !----------------------------------------
    itypcd(ii,iu) = 0
    spcond(ii,iu) = 0.d0
    itypcd(ii,iv) = 0
    spcond(ii,iv) = 0.d0
    itypcd(ii,iw) = 0
    spcond(ii,iw) = 0.d0

    ! any condensation source term
    ! associated to each turbulent variables
    ! for (k -eps) standrad turbulence model
    !----------------------------------------
    if (itytur.eq.2) then
      itypcd(ii,ik ) = 0
      spcond(ii,ik ) = 0.d0
      itypcd(ii,iep) = 0
      spcond(ii,iep) = 0.d0
    endif
    if (nscal.gt.0) then
      do iscal = 1, nscal
        if (iscal.eq.iscalt) then

          ! enthalpy value used for
          ! the explicit condensation term
          itypcd(ii,isca(iscalt)) = 1
          spcond(ii,isca(iscalt)) = hvap
        else

          ! scalar values used for
          ! the explicit condensation term
          itypcd(ii,isca(iscal)) = 1
          spcond(ii,isca(iscal)) = 0.d0
        endif
      enddo
    endif

  enddo

  if (irangp.ge.0) then
    call parrmx(nzones, zhext)
    call parrmx(nzones, ztext)
    call parrmx(nzones, zrob )
    call parrmx(nzones, zcondb)
    call parrmx(nzones, zcpb )
    call parrmx(nzones, ztpar)
  endif

endif
!< [source_terms_values]

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_boundary_mass_source_terms
