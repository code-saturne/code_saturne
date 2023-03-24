!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file pthrbm.f90
!
!> \brief Update the density \f$ \rho^{n+1}\f$  with the
!> \f$ \rho^{n-\frac{1}{2}} \f$ density with the state law and a thermodynamic
!> pressure \f$ p_{ther}^{n+1} \f$ estimated from the integral over the total
!> fluid domain of the mass conservation equation.
!>
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     dt            time step (per cell)
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, spcond is the flow rate
!>                              \f$ \Gamma_{s, cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!_______________________________________________________________________________

subroutine pthrbm &
 ( nvar   , ncesmp , nfbpcd, ncmast,                                    &
   dt     , smacel , spcond, svcond)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar , ncesmp , nfbpcd , ncmast

double precision dt(ncelet)
double precision smacel(ncesmp,nvar)
double precision spcond(nfbpcd,nvar), svcond(ncelet,nvar)

! Local variables

integer          iel , ifac

double precision new_pther
double precision ro0moy

type(var_cal_opt) :: vcopt

double precision, dimension(:), pointer :: brom, crom

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

!   pointers for the different variables

call field_get_val_s(icrom, crom)

call field_get_val_s(ibrom, brom)

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

!===============================================================================
! 2. Flow rate computation for the inlet and oulet conditions
!===============================================================================

!----------------------------------
! Update the thermodynamic pressure
! for the previous time step
!----------------------------------

pthera = pther

call compute_td_pressure_perfect_gas(nvar, ncesmp, nfbpcd, ncmast,    &
                                     dt, smacel, spcond, svcond,      &
                                     new_pther)

call cs_user_physical_properties_td_pressure(new_pther)

pther = new_pther

!===============================================================================
! 6. Thermodynamic pressure and density computation
!===============================================================================

! Update the density rho[p_ther^(n+1)] = rho^(n+1)
do iel = 1, ncel
  crom(iel) = pther/pthera*crom(iel)
enddo

! Synchronize density array
if (irangp.ge.0 .or. iperio.eq.1) then
  call synsca(crom)
endif

! Update the density at the boundary face
! with cell value for severe accident low-Mach algorithm
if (idilat.eq.3) then
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    brom(ifac) = crom(iel)
  enddo
! else with boundary values
else
  do ifac = 1, nfabor
    brom(ifac) = pther/pthera*brom(ifac)
  enddo
endif

!===============================================================================
! 7. Change the reference variable rho0
!===============================================================================

!initialization
ro0moy = 0.d0

do iel=1,ncel
  ro0moy = ro0moy + crom(iel)*volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(ro0moy)
endif

ro0 = ro0moy/voltot

!===============================================================================
! 8. Printing for severe accident low-Mach algorithm
!===============================================================================

if (vcopt%iwarni.ge.1) then
  if (cs_log_default_is_active() .eqv. .true.) then
    write (nfecra, 2002) ttcabs, pther, pthera,                            &
                         (pther-pthera)/dt(1), ro0
  endif

endif

!--------
! Formats
!--------

2002 format(/,                                                     &
   3x,'** Thermodynamic pressure computation (pthrbm): '       , /,&
   3x,'   -------------------------------------------- '       , /,&
   /,&
   6x,'-------------------------------------------------------',  &
   '---------------------------------------------------'    ,/,&
   16x,'time', 7x,'pther(n+1)',5x,'pther(n)',3x,'Dpther/Dt',       &
   7x,'ro0',                                                    /,&
   6x,'-------------------------------------------------------',  &
   '---------------------------------------------------',      &
   5(e12.5, 1x))

!----
! End
!----

return

end subroutine

!===============================================================================

!> \brief Compute the thermodynamic pressure for a perfect gas.
!>
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     dt            time step (per cell)
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, spcond is the flow rate
!>                              \f$ \Gamma_{s, cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!_______________________________________________________________________________

subroutine compute_td_pressure_perfect_gas &
  (nvar   , ncesmp , nfbpcd, ncmast, &
   dt     , smacel , spcond, svcond, &
   new_pther)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use field
use entsor
use parall
use optcal, only: isuite, ntcabs, idilat
use mesh, only: nfabor, volume, surfbn, ifabor, ncelet, ncel
use pointe, only: itypfb, icetsm, ltmast
use ppincl, only: icondv
use numvar, only: ivarfl, ipr, icrom, kbmasf
use paramx, only: ientre, i_convective_inlet, isolib, ifrent
use cstphy, only: pther, pthermax, sleak, kleak, roref, voltot, ro0, p0
use cs_tagms, only:s_metal
use cs_nz_condensation, only: ifbpcd
use cs_c_bindings

use, intrinsic :: iso_c_binding

integer, intent(in) :: nvar, ncesmp, nfbpcd, ncmast
double precision, intent(in) :: dt(ncelet)
double precision, intent(in) :: smacel(ncesmp,nvar)
double precision, intent(in) :: spcond(nfbpcd,nvar), svcond(ncelet,nvar)
double precision, intent(out) :: new_pther

integer iflmab
integer ifac, iel, ieltsm, ipcd, icmet
logical lromo
double precision rho
double precision dp
double precision debin, debout, debtot
double precision roamoy, romoy
type(var_cal_opt) :: vcopt

double precision, allocatable, dimension(:) :: surfbm

double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:), pointer :: crom, cromo

call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
call field_get_val_s(iflmab, bmasfl)
call field_get_val_s(icrom, crom)

call field_have_previous(icrom, lromo)
if (lromo) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

!---------------
! Initialization
!---------------

debin  = 0.d0
debout = 0.d0
debtot = 0.d0

! Computation of mass flux imposed on the boundary faces
do ifac = 1, nfabor
  if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then
    ! the inlet mass flux integrated in space
    debin = debin - bmasfl(ifac)
  else if (itypfb(ifac).eq.isolib.or.itypfb(ifac).eq.ifrent) then
    ! the outlet mass flux integrated in space:
    debout = debout - bmasfl(ifac)
  endif
enddo

debtot = debin + debout

! Computation of the inlet mass flux imposed on the cells volume
if (ncesmp.gt.0) then
  do ieltsm = 1, ncesmp
    iel = icetsm(ieltsm)
    debtot = debtot + smacel(ieltsm,ipr)*volume(iel)
  enddo
endif

!===============================================================================
! 3. Flow rate computation associated to the condensation phenomena
!===============================================================================

! Sink source term associated to
! the surface condensation modelling
if (nfbpcd.gt.0) then
  do ipcd = 1, nfbpcd
    ifac= ifbpcd(ipcd) + 1 ! C numbering
    iel = ifabor(ifac)
    debtot = debtot + surfbn(ifac) * spcond(ipcd,ipr)
  enddo
endif

! Sink source term associated to
! the metal structures condensation modelling
if (icondv.eq.0) then
  allocate(surfbm(ncelet))
  surfbm(:) = 0.d0

  do icmet = 1, ncmast
    iel= ltmast(icmet)
    surfbm(iel) = s_metal*volume(iel)/voltot
    debtot = debtot + surfbm(iel)*svcond(iel,ipr)
  enddo

  deallocate(surfbm)
endif

!===============================================================================
! 4. Parallelism processing
!===============================================================================

if (irangp.ge.0) then
  call parsom (debtot)
endif

!===============================================================================
! 5. Global leak
!===============================================================================

dp = pther-p0

if (dp.gt.0.d0) then
  rho = ro0
else
  rho = roref
endif

debtot = debtot - sign(1.d0,dp) * sleak * sqrt(2.d0*rho/kleak*abs(dp))

! for the first time step : rho^(n-1) = rho^(n)
if (isuite.eq.0 .and. ntcabs.eq.1) then
  do iel = 1, ncel
    cromo(iel) = crom(iel)
  enddo
endif

! initialization
roamoy = 0.d0
romoy  = 0.d0

do iel = 1, ncel
  roamoy = roamoy + cromo(iel)*volume(iel)
  romoy  = romoy  + crom(iel)*volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(roamoy)
  call parsom(romoy)
endif

! Compute the thermodynamic pressure p_ther^(n+1)
new_pther = pther*(roamoy/romoy + dt(1)*debtot/romoy)

! pthermodynamic pressure clipping in case of user venting
if (pthermax.gt.0) then
  new_pther = min(new_pther, pthermax)
endif

if (idilat.eq.3.and. vcopt%iwarni.ge.1) then
  if (cs_log_default_is_active() .eqv. .true.) then
    write (nfecra, 2003) ttcabs, roamoy/romoy, (dt(1)*debtot)/romoy
  endif
endif

!--------
! Formats
!--------

2003 format                                                            &
  (/,                                                                  &
   '   ** Perfect gas computation of average td_pressure:', /,         &
   '      -----------------------------------------------', /,         &
   /,                                                                  &
   '      -------------------------------------------------------', /, &
   '      time       rho(n-1)/rho(n)     dt.debtot/ro(n)', /,          &
   '      -------------------------------------------------------', /, &
   3(e12.5, 1x))

!----
! End
!----

return

end subroutine
