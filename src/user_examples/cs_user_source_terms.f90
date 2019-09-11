!===============================================================================
! User source terms definition.
!
! 1) Momentum equation (coupled solver)
! 2) Species transport
! 3) Turbulence (k-epsilon, k-omega, Rij-epsilon, v2-f, Spalart-Allmaras)
!===============================================================================

!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_source_terms.f90
!>
!> \brief Additional right-hand side source terms
!>
!> See \subpage cs_user_source_terms and \subpage cs_user_source_terms-scalar_in_a_channel
!> for examples.
!>
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \section ustsnv_use  Usage
!>
!> The additional source term is decomposed into an explicit part (\c crvexp) and
!> an implicit part (\c crvimp) that must be provided here.
!> The resulting equation solved by the code for a velocity is:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
!>   = \tens{crvimp} \cdot \vect{u} + \vect{crvexp}
!> \f]
!>
!> Note that \c crvexp and \c crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m/s2
!>   - crvimp is expressed in kg/s
!>
!> The \c crvexp and \c crvimp arrays are already initialized to 0
!> before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> \remark The additional force on \f$ x_i \f$ direction is given by
!>  \c crvexp(i, iel) + vel(j, iel)* crvimp(j, i).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a
!> \ref getcel command. For more info on the syntax of the \ref getcel command,
!> refer to the user manual or to the comments on the similar command
!> \ref getfbr in the routine \ref cs_user_boundary_conditions.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     ivar          index number of the current variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see
!>                               \ref cs_user_mass_source_terms)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none
!< [arg_1]
! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)
!< [arg_1]
!< [loc_var_dec_1]
! Local variables

character*80     chaine
integer          iel
double precision ckp, qdm, beta

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cvar_temperature
double precision, dimension(:), pointer ::  cpro_rom

type(var_cal_opt) :: vcopt

!< [loc_var_dec_1]

!===============================================================================
! 1. Initialization
!===============================================================================

!< [allocate_1]
! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

! --- Get variable calculation options
call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), chaine)
  write(nfecra,1000) chaine(1:8)
endif

call field_get_val_s(icrom, cpro_rom)
!< [allocate_1]

!===============================================================================
! 2. Example of arbitrary source term for component u:

!                             S = A * u + B

!            appearing in the equation under the form:

!                       rho*du/dt = S (+ standard Navier-Stokes terms)


!In the following example:
!  A = -rho*CKP
!  B =  XMMT
!
!with:
!  CKP = 1.d0 [1/s       ] (return term on velocity)
!  MMT = 100.d0 [kg/m2/s2] (momentum production by volume and time unit)
!
!which yields:
!     crvimp(1, 1, iel) = cell_f_vol(iel)* A = - cell_f_vol(iel)*(rho*CKP )
!     crvexp(1, iel) = cell_f_vol(iel)* B = cell_f_vol(iel)*(XMMT)

! ----------------------------------------------

!< [remaining_1]
ckp  = 10.d0
qdm  = 100.d0

do iel = 1, ncel
  crvimp(1, 1, iel) = - cell_f_vol(iel)*cpro_rom(iel)*ckp
enddo

do iel = 1, ncel
  crvexp(1, iel) = cell_f_vol(iel)*qdm
enddo

!< [remaining_1]

!< [boussinesq_st]

! Expension coefficient
beta = 1.d0

! Get temperature field
call field_get_val_s_by_name("temperature", cvar_temperature)

do iel = 1, ncel
  crvexp(3, iel) = cell_f_vol(iel) * ro0 * beta * (cvar_temperature(iel) - t0)
enddo

!< [boussinesq_st]

!--------
! Formats
!--------

!< [format_1]
 1000 format(' User source terms for variable ',A8,/)
!< [format_1]
!----
! End
!----
!< [deallocate_1]
! Deallocate the temporary array
deallocate(lstelt)
!< [deallocate_1]
return
end subroutine ustsnv


!===============================================================================


!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for scalar equations (user
!     scalars and specific physics scalars).

!
! Usage
! -----
! The routine is called for each scalar, user or specific physisc. It is
! therefore necessary to test the value of the scalar number iscal to separate
! the treatments of the different scalars (if (iscal.eq.p) then ....).
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code for a scalar f is:
!
!  rho*volume*df/dt + .... = crvimp*f + crvexp
!
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.[scal]/s, where [scal] is the unit of the scalar
!   - crvimp is expressed in kg/s
!
!
! The crvexp and crvimp arrays are already initialized to 0 before entering the
! the routine. It is not needed to do it in the routine (waste of CPU time).
!
! For stability reasons, Code_Saturne will not add -crvimp directly to the
! diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
! treated implicitely only if it strengthens the diagonal of the matrix.
! However, when using the second-order in time scheme, this limitation cannot
! be done anymore and -crvimp is added directly. The user should therefore test
! the negativity of crvimp by himself.
!
! When using the second-order in time scheme, one should supply:
!   - crvexp at time n
!   - crvimp at time n+1/2
!
!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command \ref getfbr in the
! routine cs_user_boundary_conditions.

! WARNING: If scalar is the temperature, the resulting equation
!          solved by the code is:
!
!  rho*Cp*cell_f_vol*dT/dt + .... = crvimp*T + crvexp
!
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in W
!   - crvimp is expressed in W/K
!

!
! STEEP SOURCE TERMS
!===================
! In case of a complex, non-linear source term, say F(f), for scalar f, the
! easiest method is to implement the source term explicitely.
!
!   df/dt = .... + F(f(n))
!   where f(n) is the value of f at time tn, the beginning of the time step.
!
! This yields :
!   crvexp = cell_f_vol*F(f(n))
!   crvimp = 0
!
! However, if the source term is potentially steep, this fully explicit
! method will probably generate instabilities. It is therefore wiser to
! partially implicit the term by writing:
!
!   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!
! This yields:
!   crvexp = cell_f_vol*( F(f(n)) - dF/df*f(n) )
!   crvimp = cell_f_vol*dF/df

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________________.
!  mode            name                 role                                           !
! _____________________________________________________________________________________!
!> \param[in]      nvar               total number of variables
!> \param[in]      nscal              total number of scalars
!> \param[in]      ncepdp             number of cells with head loss terms
!> \param[in]      ncesmp             number of cells with mass source terms
!> \param[in]      iscal              index number of the current scalar
!> \param[in]      icepdc             index number of cells with head loss terms
!> \param[in]      icetsm             index number of cells with mass source terms
!> \param[in]      itypsm             type of mass source term for each variable
!>                                     (see cs_user_mass_source_terms)
!> \param[in]      dt                 time step (per cell)
!> \param[in]      ckupdc             head loss coefficient
!> \param[in]      smacel             value associated to each variable in the mass
!>                                    source terms or mass rate
!>                                    (see cs_user_mass_source_terms)
!> \param[out]     crvexp             explicit part of the source term
!> \param[out]     crvimp             implicit part of the source term
!______________________________________________________________________________________!

subroutine ustssc &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

!< [arg_2]
! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
!< [arg_2]

!< [loc_var_dec_2]
! Local variables

character(len=80) :: chaine
integer          ivar, iiscvr,  iel
integer          ilelt, nlelt

double precision tauf, prodf, voltf, pwatt

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom

type(var_cal_opt) :: vcopt

!< [loc_var_dec_2]

!===============================================================================
! 1. Initialization
!===============================================================================

!< [allocate_2]
! Allocate a temporary array for cells selection
allocate(lstelt(ncel))
!< [allocate_2]

! --- Index number of the variable associated to scalar iscal
!< [index_2]
ivar = isca(iscal)
!< [index_2]

! --- Name of the the variable associated to scalar iscal
!< [name_2]
call field_get_label(ivarfl(ivar), chaine)
!< [name_2]

! --- Indicateur of variance scalars
!         If iscavr(iscal) = 0:
!           the scalar iscal is not a variance
!         If iscavr(iscal) > 0 and iscavr(iscal) < nscal + 1 :
!           the scalar iscal is the variance of the scalar iscavr(iscal)
!< [test_2]
iiscvr = iscavr(iscal)
!< [test_2]

! --- Density
!< [density_2]
call field_get_val_s(icrom, cpro_rom)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!< [density_2]

!===============================================================================
! 2. Example of arbitrary source term for the scalar f, 2nd scalar in the
!    calculation

!                             S = A * f + B

!            appearing in the equation under the form

!                       rho*df/dt = S (+ regular terms in the equation)


!In the following example:
!     A = - rho / tauf
!     B =   rho * prodf
!        with
!     tauf   = 10.d0  [ s  ] (dissipation time for f)
!     prodf  = 100.d0 [ [f]/s ] (production of f by unit of time)

!which yields
!     crvimp(iel) = cell_f_vol(iel)* A = - cell_f_vol(iel)*rho/tauf
!     crvexp(iel) = cell_f_vol(iel)* B =   cell_f_vol(iel)*rho*prodf

!===============================================================================

! ----------------------------------------------

!Source term applied to second scalar
!< [src_term_applied]
if (iscal.eq.2) then

   tauf  = 10.d0
   prodf = 100.d0

   do iel = 1, ncel
      crvimp(iel) = - cell_f_vol(iel)*cpro_rom(iel)/tauf
   enddo

   do iel = 1, ncel
      crvexp(iel) =   cell_f_vol(iel)*cpro_rom(iel)*prodf
   enddo

endif
!< [src_term_applied]

!===============================================================================
! 3. Example of arbitrary volumic heat term in the equation for enthalpy h

! In the considered example, a uniform volumic source of heating is imposed
! in the cells with coordinate X in [0;1.2] and Y in [3.1;4]

! The global heating power if Pwatt (in W) and the total volume of the selected
! cells is voltf (in m3)

! This yields
!     crvimp(iel) = 0
!     crvexp(iel) = cell_f_vol(iel)* Pwatt/voltf

!===============================================================================

! ----------------------------------------------

! WARNING :
! It is assumed here that the thermal scalar is an enthalpy.
! If the scalar is a temperature, PWatt does not need to be divided
! by Cp because Cp is put outside the diffusion term and multiplied
! in the temperature equation as follows:
!
!  rho*Cp*cell_f_vol*dT/dt + .... =  cell_f_vol(iel)* Pwatt/voltf

pwatt = 100.d0

! calculation of voltf

!< [ex_3_compute_voltf]
voltf  = 0.d0
call getcel('x > 0.0 and x < 1.2 and y > 3.1 and '//               &
            'y < 4.0', nlelt, lstelt)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  voltf = voltf + cell_f_vol(iel)
enddo

if (irangp.ge.0) then
  call parsom(voltf)
endif
!< [ex_3_compute_voltf]

!< [ex_3_apply]
do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! No implicit source term
  crvimp(iel) = 0.d0
! Explicit source term
  crvexp(iel) = cell_f_vol(iel)*pwatt/voltf
enddo
!< [ex_3_apply]

!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
!< [deallocate_2]
deallocate(lstelt)
!< [deallocate_2]

return
end subroutine ustssc


!===============================================================================


!===============================================================================
! Purpose:
! -------

!> \brief Additional right-hand side source terms for turbulence models
!>
!> \section cs_user_turbulence_source_terms_use  Usage
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equations solved by the code are:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\varia} + ....
!>   = \tens{crvimp} \varia + \vect{crvexp}
!> \f]
!> where \f$ \varia \f$ is the turbulence field of index \c f_id
!>
!> Note that crvexp, crvimp are defined after the Finite Volume
!> integration over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.m2/s2
!>   - crvimp is expressed in kg/s
!>
!> The crvexp, crvimp arrays are already initialized to 0 before
!> entering the routine. It is not needed to do it in the routine (waste of CPU time).
!>
!> For stability reasons, Code_Saturne will not add -crvimp directly to the
!> diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
!> treated implicitely only if it strengthens the diagonal of the matrix.
!> However, when using the second-order in time scheme, this limitation cannot
!> be done anymore and -crvimp is added directly. The user should therefore test
!> the negativity of crvimp by himself.
!>
!> When using the second-order in time scheme, one should supply:
!>   - crvexp at time n
!>   - crvimp at time n+1/2
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the \ref getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the routine
!> \ref cs_user_boundary_conditions.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     f_id          field index of the current turbulent variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref cs_user_mass_source_terms)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see
!>                               \ref cs_user_mass_source_terms)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   f_id   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

!< [arg_3]
! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          f_id

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
!< [arg_3]

!< [loc_var_dec_3]
! Local variables

integer          iel
double precision ff, tau

type(var_cal_opt) :: vcopt

character(len=80) :: fname

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom
double precision, dimension(:), pointer ::  cvar_var
!< [loc_var_dec_3]

!===============================================================================
! 1. Initialization
!===============================================================================

!< [allocate_3]
! Allocate a temporary array for cells selection
allocate(lstelt(ncel))
!< [allocate_3]

! --- Get the density array in cpro_rom
!< [dens_array_3]
call field_get_val_s(icrom, cpro_rom)
!< [dens_array_3]

! --- Get the array of the current turbulent variable and its name
!< [current_turb_3]
call field_get_val_s(f_id, cvar_var)
call field_get_name(f_id, fname)

! --- Get variable calculation options
call field_get_key_struct_var_cal_opt(f_id, vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000)
endif
!< [current_turb_3]

!===============================================================================
! 2. Example of arbitrary additional source term for turbulence models
!    (Source term on the TKE 'k' here)

!      Source term for cvar_var:
!         rho cell_f_vol d(cvar_var)/dt       = ...
!                        ... - rho*cell_f_vol*ff - rho*cell_f_vol*cvar_var/tau

!      With ff=3.d0 and tau = 4.d0

!===============================================================================

! NB the turbulence variable names are:
! - 'k' and 'epsilon' for the k-epsilon models
! - 'r11', 'r22', 'r33', 'r12', 'r13', 'r23' and 'epsilon'
!    for the Rij-epsilon LRR and SSG
! - 'r11', 'r22', 'r33', 'r12', 'r13', 'r23', 'epsilon' and 'alpha' for the EBRSM
! - 'k', 'epsilon', 'phi' and 'f_bar' for the phi-model
! - 'k', 'epsilon', 'phi' and 'alpha' for the Bl-v2-k model
! - 'k' and 'omega' for the k-omega turbulence model
! - 'nu_tilda' for the Spalart Allmaras model

!< [rem_code_3]

if (trim(fname).eq.'k') then

  ff  = 3.d0
  tau = 4.d0

  ! --- Explicit source terms
  do iel = 1, ncel
    crvexp(iel) = -cpro_rom(iel)*cell_f_vol(iel)*ff
  enddo

  ! --- Implicit source terms
  !        crvimp is already initialized to 0, no need to set it here
  do iel = 1, ncel
    crvimp(iel) = -cpro_rom(iel)*cell_f_vol(iel)/tau
  enddo

endif

!< [rem_code_3]

!--------
! Formats
!--------

!< [format_3]
 1000 format(' User source terms for turbulence model',/)
!< [format_3]

!----
! End
!----
!< [deallocate_3]
! Deallocate the temporary array
deallocate(lstelt)
!< [deallocate_3]

return
end subroutine cs_user_turbulence_source_terms
