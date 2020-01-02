!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cs_user_boundary_mass_source_terms.f90
!>
!> \brief Source terms associated at the boundary faces and the neighboring
!>        cells with surface condensation.
!>
!> \par Examples of settings for boundary condensation mass source terms
!>      Examples are available
!>      \ref condens_h_boundary "here".
!>
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
use cs_tagmr

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
integer          ieltcd
integer          ifac, iel, iscal
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
!     potentially becomea condensation  source term.
!===============================================================================


  !--To Select the cells with condensation source term
  !---------------------------------------------------

  izone = 0
  ieltcd = 0

  ! Cells with a boundary face of color 60
  call getfbr('60',nlelt,lstelt)

  izone = izone + 1

  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)
    izftcd(iel) = izone
    ieltcd = ieltcd + 1
    if (iappel.eq.2) ifbpcd(ieltcd) = ifac
  enddo

endif

! For iappel = 1,
! Specification of nfbpcd.
if (iappel.eq.1) then
  nfbpcd = ieltcd
endif
!< [zones_definition]

!< [model_settings]
!===============================================================================
! Parameters padding of the 1-D thermal model and condensation model
! ------------------------------------------------------------------
! Both models can be enabled and coupled together or
! the condensation model can be used without enabling the 1-D thermal model.
! In this latter case a constant wall temperature must be specified by the
! user at the cold wall (at the third call (iappel=3) tpar=tpar0).
!===============================================================================
if (iappel.eq.2) then

  if (icondb.eq.0) then

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

    icophc = 3

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

    icophg = 3

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

    itag1d = 1

    ! Wall temperature computed by a 1-D thermal model
    ! with a implicit scheme and variable over time.
    ! ------------------------------------------------
    ! Remark : the wall temperature is in unit [°C].
    if(itag1d.eq.1) then

      !---------------------------------------------------
      ! Numerical parameters used by the 1-D thermal model
      !---------------------------------------------------

      ! (theta) parameter of the implicit scheme
      theta = 1.d0
      ! (dxmin) First cell size of the 1D mesh
      ! -> with (dxmin.eq.0) the Dx is constant
      ! -> with (dxmin.ne.0) the Dx is variable
      dxmin = 0.d0
      ! (nmur) space steps number of the 1D mesh
      nmur = 10
      ! (epais) thickness of the 1D wall
      epais = 0.024d0

      !-------------------------------------------
      !Initial condition of the 1-D thermal model
      !-------------------------------------------
      tpar0 = 26.57d0

   endif

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

  if (icondb.eq.0) then
    if(itag1d.eq.1) then
      !-------------------------------------------
      !Boundary conditions of the 1-D thermal model
      !-------------------------------------------
      hext =  1.d+8 ; text = 26.57d0
      ! --------------------------------------------
      ! Physical properties of the concrete material
      ! --------------------------------------------
      ! (rob) density (kg.m-3)
      rob   = 8000.d0
      ! (condb) thermal conductivity (W.m-1.C-1)
      condb = 12.8d0
      ! (cpb)   Specific heat (J.kg-1.C-1)
      cpb = 500.0d0

    else
      ! Wall temperature imposed as constant
      ! with a value specified by the user
      tpar = 26.57d0
    endif
  endif

  ! To fill the spcond(nfbpcd,ivar) array
  ! if we want to specify a variable value
  !---------------------------------------
  do ieltcd = 1, nfbpcd

    ifac = ifbpcd(ieltcd)
    iel  = ifabor(ifac)

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
    itypcd(ieltcd,iu) = 0
    spcond(ieltcd,iu) = 0.d0
    itypcd(ieltcd,iv) = 0
    spcond(ieltcd,iv) = 0.d0
    itypcd(ieltcd,iw) = 0
    spcond(ieltcd,iw) = 0.d0

    ! any condensation source term
    ! associated to each turbulent variables
    ! for (k -eps) standrad turbulence model
    !----------------------------------------
    if (itytur.eq.2) then
      itypcd(ieltcd,ik ) = 0
      spcond(ieltcd,ik ) = 0.d0
      itypcd(ieltcd,iep) = 0
      spcond(ieltcd,iep) = 0.d0
    endif
    if (nscal.gt.0) then
      do iscal = 1, nscal
        if (iscal.eq.iscalt) then

          ! enthalpy value used for
          ! the explicit condensation term
          itypcd(ieltcd,isca(iscalt)) = 1
          spcond(ieltcd,isca(iscalt)) = hvap
        else

          ! scalar values used for
          ! the explicit condensation term
          itypcd(ieltcd,isca(iscal)) = 1
          spcond(ieltcd,isca(iscal)) = 0.d0
        endif
      enddo
    endif


  enddo

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
