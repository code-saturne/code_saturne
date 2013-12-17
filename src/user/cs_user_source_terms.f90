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
! Copyright (C) 1998-2013 EDF S.A.
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
!> over the cells, so they include the "volume" term. More precisely:
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
!>                               (see \ref ustsma)
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtpa          calculated variables at cell centers
!>                               (preceding time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see \ref ustsma)
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
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

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)
double precision, dimension(:), pointer ::  crom
! Local variables

character*80     chaine
integer          iel
double precision ckp, qdm

integer, allocatable, dimension(:) :: lstelt

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

call field_get_val_s(icrom, crom)

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

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

ckp  = 10.d0
qdm  = 100.d0

do iel = 1, ncel
  crvimp(1, 1, iel) = - volume(iel)*crom(iel)*ckp
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
   dt     , rtpa   , rtp    , propce ,                            &
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
!  rho*Cp*volume*dT/dt + .... = crvimp*T + crvexp
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
!   crvexp = volume*F(f(n))
!   crvimp = 0
!
! However, if the source term is potentially steep, this fully explicit
! method will probably generate instabilities. It is therefore wiser to
! partially implicit the term by writing:
!
!   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!
! This yields:
!   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!   crvimp = volume*dF/df

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
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time step)                         !
! rtp              ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (current time step)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
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

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

character*80     chaine
integer          ivar, iiscvr,  iel
integer          ilelt, nlelt

double precision tauf, prodf, volf, pwatt

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom

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

! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

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
!     crvimp(iel) = volume(iel)* A = - volume(iel)*rho/tauf
!     crvexp(iel) = volume(iel)* B =   volume(iel)*rho*prodf

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
      crvimp(iel) = - volume(iel)*crom(iel)/tauf
   enddo

   do iel = 1, ncel
      crvexp(iel) =   volume(iel)*crom(iel)*prodf
   enddo

endif

!===============================================================================
! 3. Example of arbitrary volumic heat term in the equation for enthalpy h

! In the considered example, a uniform volumic source of heating is imposed
! in the cells with coordinate X in [0;1.2] and Y in [3.1;4]

! The global heating power if Pwatt (in W) and the total volume of the concerned
! cells is volf (in m3)

! This yields
!     crvimp(iel) = 0
!     crvexp(iel) = volume(iel)* Pwatt/volf

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
!  rho*Cp*volume*dT/dt + .... =  volume(iel)* Pwatt/volf
!

pwatt = 100.d0

! calculation of volf

volf  = 0.d0
call getcel('x > 0.0 and x < 1.2 and y > 3.1 and '//               &
            'y < 4.0',nlelt,lstelt)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  volf = volf + volume(iel)
enddo

if (irangp.ge.0) then
  call parsom(volf)
endif

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! No implicit source term
  crvimp(iel) = 0.d0
! Explicit source term
  crvexp(iel) = volume(iel)*pwatt/volf
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


subroutine ustske &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , tinstk , divu   ,                            &
   crkexp , creexp , crkimp , creimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for k and epsilon equations
!    when using:
!     - k-epsilon model (ITURB=20)
!     - k-epsilon Linear Production model (ITURB=21)
!     - v2-f phi-model (ITURB=50)
!
! Usage
! -----
!
! The additional source term is decomposed into an explicit part (crkexp,creexp) and
! an implicit part (crkimp,creimp) that must be provided here.
! The resulting equations solved by the code are:
!
!  rho*volume*d(k)/dt   + .... = crkimp*k   + crkexp

!  rho*volume*d(eps)/dt + .... = creimp*eps + creexp
!
! Note that crkexp, crkimp, creexp and creimp are defined after the Finite Volume
! integration over the cells, so they include the "volume" term. More precisely:
!   - crkexp is expressed in kg.m2/s3
!   - creexp is expressed in kg.m2/s4
!   - crkimp is expressed in kg/s
!   - creimp is expressed in kg/s
!
! The crkexp, crkimp, creexp and creimp arrays are already initialized to 0 before
! entering the routine. It is not needed to do it in the routine (waste of CPU time).
!
! For stability reasons, Code_Saturne will not add -crkimp directly to the
! diagonal of the matrix, but Max(-crkimp,0). This way, the crkimp term is
! treated implicitely only if it strengthens the diagonal of the matrix.
! However, when using the second-order in time scheme, this limitation cannot
! be done anymore and -crkimp is added directly. The user should therefore test
! the negativity of crkimp by himself.
! The same mechanism applies to cveimp.
!
! When using the second-order in time scheme, one should supply:
!   - crkexp and creexp at time n
!   - crkimp and creimp at time n+1/2
!
! When entering the routine, two additional work arrays are already set for
! potential user need:
!   tinstk =  2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!          +  (2 S12)**2 + (2 S13)**2 + (2 S23)**2
!
!          where Sij = (dUi/dxj+dUj/dxi)/2
!
!   divu = du/dx + dv/dy + dw/dz



!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! tinstk           ! ra ! <-- ! tubulent production term (see comment above)   !
! divu             ! ra ! <-- ! velocity divergence (see comment above)        !
! crkexp           ! ra ! --> ! explicit part of the source term for k         !
! creexp           ! ra ! --> ! explicit part of the source term for epsilon   !
! crkimp           ! ra ! --> ! implicit part of the source term for k         !
! creimp           ! ra ! --> ! implicit part of the source term for epsilon   !
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

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision tinstk(ncelet), divu(ncelet)
double precision crkexp(ncelet), crkimp(ncelet)
double precision creexp(ncelet), creimp(ncelet)

! Local variables

integer          iel
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom
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


! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

if (iwarni(ik).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for k and epsilon

!      Source term for k :
!         rho volume d(k)/dt       = ...
!                        ... - rho*volume*ff*epsilon - rho*volume*k/tau

!      Source term for epsilon :
!         rho volume d(epsilon)/dt = ...
!                        ... + rho*volume*xx

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

! --- Explicit source terms

ff  = 3.d0
tau = 4.d0
xx  = 2.d0

do iel = 1, ncel
  crkexp(iel) = -crom(iel)*volume(iel)*ff*rtpa(iel,iep)
  creexp(iel) =  crom(iel)*volume(iel)*xx
enddo

! --- Implicit source terms
!        creimp is already initialized to 0, no need to set it here

do iel = 1, ncel
  crkimp(iel) = -crom(iel)*volume(iel)/tau
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for k and epsilon',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustske


!===============================================================================


subroutine ustskw &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , s2kw   , divukw ,                            &
   gkgw   , ggrho  , xf1    ,                                     &
   crkexp , crwexp , crkimp , crwimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for k and omega equations
!    when using k-omega SST (ITURB=60)
!
! Usage
! -----
!
! The additional source term is decomposed into an explicit part (crkexp,crwexp) and
! an implicit part (crkimp,crwimp) that must be provided here.
! The resulting equations solved by the code are:
!
!  rho*volume*d(k)/dt     + .... = crkimp*k     + crkexp

!  rho*volume*d(omega)/dt + .... = crwimp*omega + crwexp
!
! Note that crkexp, crkimp, crwexp and crwimp are defined after the Finite Volume
! integration over the cells, so they include the "volume" term. More precisely:
!   - crkexp is expressed in kg.m2/s3
!   - crwexp is expressed in kg/s2
!   - crkimp is expressed in kg/s
!   - crwimp is expressed in kg/s
!
! The crkexp, crkimp, crwexp and crwimp arrays are already initialized to 0 before
! entering the routine. It is not needed to do it in the routine (waste of CPU time).
!
! For stability reasons, Code_Saturne will not add -crkimp directly to the
! diagonal of the matrix, but Max(-crkimp,0). This way, the crkimp term is
! treated implicitely only if it strengthens the diagonal of the matrix.
! However, when using the second-order in time scheme, this limitation cannot
! be done anymore and -crkimp is added directly. The user should therefore test
! the negativity of crkimp by himself.
! The same mechanism applies to cvwimp.
!
! When using the second-order in time scheme, one should supply:
!   - crkexp and crwexp at time n
!   - crkimp and crwimp at time n+1/2
!
! When entering the routine, some additional work arrays are already set for
! potential user need:
!   s2kw   =  2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!          +  (2 S12)**2 + (2 S13)**2 + (2 S23)**2
!
!            where Sij = (dUi/dxj+dUj/dxi)/2
!
!   divukw = du/dx + dv/dy + dw/dz

!   gkgw = grad(k).grad(omega)

!   ggrho = -g.grad(rho)/prdtur/rho if igrake>0
!         = 0                       if igrake=0
!            where prdtur is the turbulent Prandtl number

!   xf1   = k-eps/k-omega blending coefficient

!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! s2kw             ! ra ! <-- ! turbulent production term (see comment above)  !
! divukw           ! ra ! <-- ! velocity divergence (see comment above)        !
! gkgw             ! ra ! <-- ! grad(k).grad(omega)                            !
! ggrho            ! ra ! <-- ! gravity source term (see comment above)        !
! xf1              ! ra ! <-- ! k-eps/k-w blending function (see comment above)!
! crkexp           ! ra ! --> ! explicit part of the source term for k         !
! crwexp           ! ra ! --> ! explicit part of the source term for omega     !
! crkimp           ! ra ! --> ! implicit part of the source term for k         !
! crwimp           ! ra ! --> ! implicit part of the source term for omega     !
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

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision s2kw(ncelet)  , divukw(ncelet)
double precision gkgw(ncelet)  , ggrho(ncelet), xf1(ncelet)
double precision crkexp(ncelet), crkimp(ncelet)
double precision crwexp(ncelet), crwimp(ncelet)

! Local variables

integer          iel
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom
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


! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

if (iwarni(ik).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for k and omega

!      Source term for k :
!         rho volume d(k)/dt       = ...
!                        ... - rho*volume*ff*omega - rho*volume*k/tau

!      Source term for omega :
!         rho volume d(omega)/dt = ...
!                        ... + rho*volume*xx

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

! --- Explicit source terms

ff  = 3.d0
tau = 4.d0
xx  = 2.d0

do iel = 1, ncel
  crkexp(iel) = -crom(iel)*volume(iel)*ff*rtpa(iel,iomg)
  crwexp(iel) =  crom(iel)*volume(iel)*xx
enddo

! --- Implicit source terms
!        crwimp is already initialized to 0, no need to set it here

do iel = 1, ncel
  crkimp(iel) = -crom(iel)*volume(iel)/tau
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for k and omega',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustskw


!===============================================================================


subroutine ustsri &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smcelp , gamma  , gradv  , produc ,                   &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for the equations of Rij and
!    epsilon in Rij-LRR model (ITURB=30) or Rij-SSG model (ITURB=31)

!
! Usage
! -----
! The routine is called for each variable Rij and epsilon. It is therefore
! necessary to test the value of the variable ivar to separate the treatments
! of the variables ivar=ir11, ir22, ir33, ir12,
! ir13, ir23, or iep.
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code for a variable var is:
!
!  rho*volume*d(var)/dt + .... = crvimp*(var) + crvexp
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.m2/s3 for Rij
!   - crvexp is expressed in kg.m2/s4 for epsilon
!   - crvimp is expressed in kg/s for Rij and epsilon
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
! When entering the routine, two additional work arrays are already set for
! potential user need:
!   produc is the turbulent production term
!          produc(6,ncelet) = (P11, P22, P33, P12, P13, P23)
!          with Pij=-Rik.dUj/dxk - Rjk.dUi/dxk
!
!   gradv is the
!          gradv(j,i,ncelet) = dUi/dxj


! WARNING
! =======
! produc is only allocated and defined for the Rij-LRR model (ITURB=30)
! gradv is only allocated and defined for the Rij-SSG model (ITURB=31)
! DO NOT USE produc WITH THE SSG MODEL
! DO NOT USE gradv WITH THE LRR MODEL


!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! ivar             ! i  ! <-- ! index number of the current variable           !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smcelp(ncelet)   ! ra ! <-- ! value of variable ivar associated to mass      !
!                  ! ra !     !  source term (see ustsma)                      !
! gamma(ncelet)    ! ra ! <-- ! volumic rate of mass source term               !
! gradv            ! ra ! <-- ! velocity gradient (only for iturb=31)          !
! produc(6,ncelet) ! ra ! <-- ! turbulent production term (only for iturb=30)  !
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
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itpsmp(ncesmp)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision gradv(3, 3, ncelet), produc(6,ncelet)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

integer          iel
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom
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


! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

if (iwarni(ir11).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for R11 and epsilon

!      Source term for R11 :
!         rho volume d(R11)/dt       = ...
!                        ... - rho*volume*ff*epsilon - rho*volume*R11/tau

!      Source term for epsilon :
!         rho volume d(epsilon)/dt = ...
!                        ... + rho*volume*xx

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------


! ---  For R11
!      -------

if (ivar.eq.ir11) then

  ff  = 3.d0
  tau = 4.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) = -crom(iel)*volume(iel)                 &
                                               *ff*rtpa(iel,iep)
  enddo

!    -- Implicit source term

  do iel = 1, ncel
    crvimp(iel) = -crom(iel)*volume(iel)/tau
  enddo


! ---  For epsilon
!      -----------

elseif (ivar.eq.iep) then

  xx  = 2.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) =  crom(iel)*volume(iel)*xx
  enddo

!    -- Implicit source term
!        crvimp is already initialized to 0, no need to set it here


endif

!--------
! Formats
!--------

 1000 format(' User source terms for Rij and epsilon',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsri


!===============================================================================


subroutine ustsv2 &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , produc , gphigk ,                            &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for the equations of phi and f_bar
!    with the v2-f phi-model turbulence model (ITURB=50)
!
! Usage
! -----
! The routine is called for both phi and f_bar. It is therefore necessary
! to test the value of the variable ivar to separate the treatments of the
! ivar=iphi or ivar=ifb.
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code are as follows:
!
! For f_bar (ivar=ifb):
!  volume*div(grad(f_bar))= ( volume*f_bar + ..... + crvimp*f_bar + crvexp )/L^2
!
! For phi (ivar=iphi)
!  rho*volume*d(phi)/dt + .... = crvimp*phi + crvexp

!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in m3/s for f_bar
!   - crvexp is expressed in kg/s for phi
!   - crvimp is expressed in m3 for f_bar
!   - crvimp is expressed in kg/s for phi
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

! When entering the routine, two additional work arrays are already set for
! potential user need:
!   produc = 2*mu_t*Sij*Sij -2/3*rho*k*div(u) -2/3*mu_t*div(u)**2
!          + gravity term (when activated)
!          (produc is the production term in the equation for k)
!
!   gphigk = grad(phi).grad(k)

!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! ivar             ! i  ! <-- ! index number of the current variable           !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! produc(ncelet)   ! ra ! <-- ! production term for k                          !
! gphigk(ncelet)   ! ra ! <-- ! grad(phi).grad(k)                              !
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
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision produc(ncelet), gphigk(ncelet)

! Local variables

integer          iel
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom
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


! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

if (iwarni(ifb).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for k and epsilon

!      Source term for f_bar :
!         volume div(grad f_barre) = (... + volume*xx )/L^2

!      Source term for phi :
!         rho volume d(phi)/dt       = ...
!                        ... + rho*volume*ff*f_bar - rho*volume*phi/tau

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

! ---  For f_bar
!      ---------

if (ivar.eq.ifb) then

  xx  = 2.d0

!    -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) =  volume(iel)*xx
  enddo

!    -- Implicit source term
!
!       crvimp is already initialized to 0, no need to set it here


! ---  For phi
!      -------

elseif (ivar.eq.iphi) then

  ff  = 3.d0
  tau = 4.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) = crom(iel)*volume(iel)*ff*rtpa(iel,ifb)
  enddo

!    -- Implicit source term

  do iel = 1, ncel
    crvimp(iel) = -crom(iel)*volume(iel)/tau
  enddo


endif

!--------
! Formats
!--------

 1000 format(' User source terms for phi and f_bar',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsv2


!===============================================================================


subroutine ustssa &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , tinssa , divu   ,                            &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for the viscosity-like variable
!    when using the Spalart-Allmaras model (ITURB=70)
!
! Usage
! -----
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equations solved by the code are:
!
!  rho*volume*d(nusa)/dt   + .... = crvimp*nusa   + crvexp
!
! Note that crvexp, crvimp are defined after the Finite Volume
! integration over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.m2/s2
!   - crvimp is expressed in kg/s
!
! The crvexp, crvimp arrays are already initialized to 0 before
! entering the routine. It is not needed to do it in the routine (waste of CPU time).
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
! When entering the routine, two additional work arrays are already set for
! potential user need:
!   tinssa =  2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!          +  (2 S12)**2 + (2 S13)**2 + (2 S23)**2
!
!          where Sij = (dUi/dxj+dUj/dxi)/2
!
!   divu = du/dx + dv/dy + dw/dz



!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! cs_user_boundary_conditions.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! tinssa           ! ra ! <-- ! tubulent production term (see comment above)   !
! divu             ! ra ! <-- ! velocity divergence (see comment above)        !
! crvexp           ! ra ! --> ! explicit part of the source term for k         !
! crvimp           ! ra ! --> ! implicit part of the source term for k         !
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

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision tinssa(ncelet), divu(ncelet)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

integer          iel
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer ::  crom
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


! --- Index number of the density in the propce array
call field_get_val_s(icrom, crom)

if (iwarni(inusa).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for nusa

!      Source term for k :
!         rho volume d(k)/dt       = ...
!                        ... - rho*volume*ff - rho*volume*nusa/tau

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

if (.true.) return

! ----------------------------------------------

! --- Explicit source terms

ff  = 3.d0
tau = 4.d0
xx  = 2.d0

do iel = 1, ncel
  crvexp(iel) = -crom(iel)*volume(iel)*ff
enddo

! --- Implicit source terms
!        crvimp is already initialized to 0, no need to set it here

do iel = 1, ncel
  crvimp(iel) = -crom(iel)*volume(iel)/tau
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for Spalart-Allmaras',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustssa
