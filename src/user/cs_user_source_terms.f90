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
! Copyright (C) 1998-2014 EDF S.A.
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
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \section use  Usage
!>
!> The additional source term is decomposed into an explicit part (\c crvexp) and
!> an implicit part (\c crvimp) that must be provided here.
!> The resulting equation solved by the code for a velocity is:
!> \f[
!>  \rho \norm{\vol{\celli}} \DP{\vect{u}} + ....
!>   = \tens{crvimp} \vect{u} + \vect{crvexp}
!> \f]
!>
!> Note that \c crvexp and \c crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volf" term. More precisely:
!>   - crvexp is expressed in kg.m/s2
!>   - crvimp is expressed in kg/s
!>
!> The \c crvexp and \c crvimp arrays are already initialized to 0
!> before entering the
!> the routine. It is not needed to do it in the routine (waste of CPU time).
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
!>                               source terms or mass rate (see \ref cs_user_mass_source_terms)
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables

character*80     chaine
integer          iel
double precision ckp, qdm

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

if (iwarni(ivar).ge.1) then
  call field_get_label(ivarfl(ivar), chaine)
  write(nfecra,1000) chaine(1:8)
endif

call field_get_val_s(icrom, cpro_rom)

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
!     crvimp(1, 1, iel) = volf(iel)* A = - volf(iel)*(rho*CKP )
!     crvexp(1, iel) = volf(iel)* B = volf(iel)*(XMMT)

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

ckp  = 10.d0
qdm  = 100.d0

do iel = 1, ncel
  crvimp(1, 1, iel) = - volf(iel)*cpro_rom(iel)*ckp
enddo

do iel = 1, ncel
  crvexp(1, iel) = volf(iel)*qdm
enddo

!--------
! Formats
!--------

 1000 format(' User source termes for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsnv


!===============================================================================


!===============================================================================


subroutine ustssc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

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
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

! WARNING: If scalar is the temperature, the resulting equation
!          solved by the code is:
!
!  rho*Cp*volf*dT/dt + .... = crvimp*T + crvexp
!
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volf" term. More precisely:
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
!   crvexp = volf*F(f(n))
!   crvimp = 0
!
! However, if the source term is potentially steep, this fully explicit
! method will probably generate instabilities. It is therefore wiser to
! partially implicit the term by writing:
!
!   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!
! This yields:
!   crvexp = volf*( F(f(n)) - dF/df*f(n) )
!   crvimp = volf*dF/df

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! iscal            ! i  ! <-- ! index number of the current scalar             !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see cs_user_mass_source_terms)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see cs_user_mass_source_terms)        !
! crvexp           ! ra ! --> ! explicit part of the source term               !
! crvimp           ! ra ! --> ! implicit part of the source term               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

character*80     chaine
integer          ivar, iiscvr,  iel
integer          ilelt, nlelt

double precision tauf, prodf, voltf, pwatt

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))


! --- Index number of the variable associated to scalar iscal
ivar = isca(iscal)

! --- Name of the the variable associated to scalar iscal
call field_get_label(ivarfl(ivar), chaine)

! --- Indicateur of variance scalars
!         If iscavr(iscal) = 0:
!           the scalar iscal is not a variance
!         If iscavr(iscal) > 0 and iscavr(iscal) < nscal + 1 :
!           the scalar iscal is the variance of the scalar iscavr(iscal)
iiscvr = iscavr(iscal)

! --- Density
call field_get_val_s(icrom, cpro_rom)

if (iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif


!===============================================================================
! 2. Example of arbitrary source term for the scalar f, 2nd scalar in the
!    calculation

!                             S = A * f + B

!            appearing in the equation under the form

!                       rho*df/dt = S (+ regular terms in the equation)


!In the following example:
!     A = - rho / tauf
!     B =   rho * prodf
!        AVEC
!     tauf   = 10.d0  [ s  ] (dissipation time for f)
!     prodf  = 100.d0 [ [f]/s ] (production of f by unit of time)

!which yields
!     crvimp(iel) = volf(iel)* A = - volf(iel)*rho/tauf
!     crvexp(iel) = volf(iel)* B =   volf(iel)*rho*prodf

!===============================================================================


! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

!Source term applied to second scalar
if (iscal.eq.2) then

   tauf  = 10.d0
   prodf = 100.d0

   do iel = 1, ncel
      crvimp(iel) = - volf(iel)*cpro_rom(iel)/tauf
   enddo

   do iel = 1, ncel
      crvexp(iel) =   volf(iel)*cpro_rom(iel)*prodf
   enddo

endif

!===============================================================================
! 3. Example of arbitrary volumic heat term in the equation for enthalpy h

! In the considered example, a uniform volumic source of heating is imposed
! in the cells with coordinate X in [0;1.2] and Y in [3.1;4]

! The global heating power if Pwatt (in W) and the total volume of the concerned
! cells is voltf (in m3)

! This yields
!     crvimp(iel) = 0
!     crvexp(iel) = volf(iel)* Pwatt/voltf

!===============================================================================


! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

! WARNING :
! It is assumed here that the thermal scalar is an enthalpy.
! If the scalar is a temperature, PWatt does not need to be devided
! by Cp because Cp is put outside the diffusion term and multiply
! the temperature equation as follows:
!
!  rho*Cp*volf*dT/dt + .... =  volf(iel)* Pwatt/voltf
!

pwatt = 100.d0

! calculation of voltf

voltf  = 0.d0
call getcel('x > 0.0 and x < 1.2 and y > 3.1 and '//               &
            'y < 4.0',nlelt,lstelt)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  voltf = voltf + volf(iel)
enddo

if (irangp.ge.0) then
  call parsom(voltf)
endif

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! No implicit source term
  crvimp(iel) = 0.d0
! Explicit source term
  crvexp(iel) = volf(iel)*pwatt/voltf
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustssc


!===============================================================================


!===============================================================================
! Purpose:
! -------

!> \brief Additional right-hand side source terms for turbulence models
!>
!> \section use  Usage
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
!>                               source terms or mass rate (see \ref cs_user_mass_source_terms)
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

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          f_id

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

integer          iel
double precision ff, tau

type(var_cal_opt) vcopt

character*80     fname

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  cpro_rom
double precision, dimension(:), pointer ::  cvar_var
!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (.true.) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

! --- Get the density array in cpro_rom
call field_get_val_s(icrom, cpro_rom)


! --- Get the array of the current turbulent variable and its name
call field_get_val_s(f_id, cvar_var)
call field_get_name(f_id, fname)

! --- Get variable calculation options
call field_get_key_struct_var_cal_opt(f_id, vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for turbulence models
!    (Source term on the TKE 'k' here)

!      Source term for cvar_var:
!         rho volf d(cvar_var)/dt       = ...
!                        ... - rho*volf*ff - rho*volf*cvar_var/tau

!      With ff=3.d0 and tau = 4.d0

!===============================================================================

! NB the turbulence varaible names are:
! - 'k' and 'epsilon' for the k-epsilon models
! - 'r11', 'r22', 'r33', 'r12', 'r13', 'r23' and 'epsilon'
!    for the Rij-epsilon LRR and SSG
! - 'r11', 'r22', 'r33', 'r12', 'r13', 'r23', 'epsilon' and 'alpha' for the EBRSM
! - 'k', 'epsilon', 'phi' and 'f_bar' for the phi-model
! - 'k', 'epsilon', 'phi' and 'alpha' for the Bl-v2-k model
! - 'k' and 'omega' for the k-omega turbulence model
! - 'nu_tilda' for the Spalart Allmaras model

if (.false.) then
  if (trim(fname).eq.'k') then

    ff  = 3.d0
    tau = 4.d0

    ! --- Explicit source terms
    do iel = 1, ncel
      crvexp(iel) = -cpro_rom(iel)*volf(iel)*ff
    enddo

    ! --- Implicit source terms
    !        crvimp is already initialized to 0, no need to set it here
    do iel = 1, ncel
      crvimp(iel) = -cpro_rom(iel)*volf(iel)/tau
    enddo

  endif
endif

!--------
! Formats
!--------

 1000 format(' User source terms for turbulence model',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_turbulence_source_terms
