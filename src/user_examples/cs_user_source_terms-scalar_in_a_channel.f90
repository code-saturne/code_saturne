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

!> \file cs_user_source_terms-scalar_in_a_channel.f90
!>
!> \brief User source terms for a scalar in a channel example.
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
!> \param[in]     ncssmp        number of cells with mass source terms
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
!================

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

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables

character*80     chaine
integer          iel
double precision ckp, qdm

double precision, dimension(:), pointer :: cpro_rom

type(var_cal_opt) :: vcopt

!===============================================================================


!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
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
!     crvimp(1, 1, iel) = volume(iel)* A = - volume(iel)*(rho*CKP )
!     crvexp(1, iel) = volume(iel)* B = volume(iel)*(XMMT)

! ----------------------------------------------

ckp  = 0.d0
qdm  = 1.d0

do iel = 1, ncel
  crvimp(1, 1, iel) = - volume(iel)*cpro_rom(iel)*ckp
enddo

do iel = 1, ncel
  crvexp(1, iel) = volume(iel)*qdm
enddo

!--------
! Formats
!--------

 1000 format(' User source termes for variable ',A8,/)

!----
! End
!----

return
end subroutine ustsnv


!===============================================================================



!===============================================================================
! Purpose:
! -------

!    User subroutine.
!> \file cs_user_source_terms-scalar_in_a_channel.f90
!> \brief    Additional right-hand side source terms for scalar equations (user
!>     scalars and specific physics scalars).
!>
!>
!> Usage
!> -----
!> The routine is called for each scalar, user or specific physisc. It is
!> therefore necessary to test the value of the scalar number iscal to separate
!> the treatments of the different scalars (if (iscal.eq.p) then ....).
!>
!> The additional source term is decomposed into an explicit part (crvexp) and
!> an implicit part (crvimp) that must be provided here.
!> The resulting equation solved by the code for a scalar f is:
!>
!>   \f[ \rho*volume*\frac{df}{dt} + .... = crvimp*f + crvexp \f]
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in kg.[scal]/s, where [scal] is the unit of the scalar
!>   - crvimp is expressed in kg/s
!>
!>
!> The crvexp and crvimp arrays are already initialized to 0 before entering the
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
!>
!> The selection of cells where to apply the source terms is based on a getcel
!> command. For more info on the syntax of the getcel command, refer to the
!> user manual or to the comments on the similar command \ref getfbr in the
!> routine cs_user_boundary_conditions.

!> WARNING: If scalar is the temperature, the resulting equation
!>          solved by the code is:
!>
!>  rho*Cp*volume*dT/dt + .... = crvimp*T + crvexp
!>
!>
!> Note that crvexp and crvimp are defined after the Finite Volume integration
!> over the cells, so they include the "volume" term. More precisely:
!>   - crvexp is expressed in W
!>   - crvimp is expressed in W/K
!>

!>
!> STEEP SOURCE TERMS
!>===================
!> In case of a complex, non-linear source term, say F(f), for scalar f, the
!> easiest method is to implement the source term explicitely.
!>
!>   df/dt = .... + F(f(n))
!>   where f(n) is the value of f at time tn, the beginning of the time step.
!>
!> This yields :
!>   crvexp = volume*F(f(n))
!>   crvimp = 0
!>
!> However, if the source term is potentially steep, this fully explicit
!> method will probably generate instabilities. It is therefore wiser to
!> partially implicit the term by writing:
!>
!>   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!>
!> This yields:
!>   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!>   crvimp = volume*dF/df
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncssmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]                    (see cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                    source terms or mass rate (see cs_user_mass_source_terms)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!

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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
!< [loc_var_scal]
! Local variables

integer          iel

double precision, dimension(:,:), pointer :: cvar_vel

double precision ubulk, flux_tot
integer          ifac
double precision ubulkm
data             ubulkm /0.d0/
save             ubulkm
!< [loc_var_scal]
!===============================================================================

if (iscal.ne.iscalt) return

!===============================================================================
! 1. Initialization
!===============================================================================
!< [init_scal]
ubulk= 0.d0
!< [init_scal]
! Map field arrays
!< [map_field_scal]
call field_get_val_v(ivarfl(iu), cvar_vel)

!< [map_field_scal
! Bulk mean velocity (x component)
do iel = 1, ncel
  ubulk = ubulk + cvar_vel(1, iel) * cell_f_vol(iel)
enddo
!< [map_field_scal]
! Parallelism
!< [parallelism_scal]
if (irangp.ge.0) then
  call parsom(ubulk)
endif

ubulk = max(abs(ubulk) / voltot, 1d-12)

! Flux x Total surface / (rho Cp)
flux_tot = 1.d0
do iel = 1, ncel
  crvimp(iel) = 0.d0
  crvexp(iel) = cell_f_vol(iel)*cvar_vel(1, iel) * flux_tot / ubulk
enddo
!< [parallelism_scal]
!--------
! Formats
!--------
!< [format_scal]
 1000 format(' User source terms for variable ',A8,/)
!< [format_scal]
!----
! End
!----
!< [end_scal]
return
end subroutine ustssc
!< [end_scal]
