!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!> \file cs_user_physical_properties.f90
!> \brief Definition of physical variable laws.
!>
!> usphyv
!> \brief Definition of physical variable laws.
!>
!> \section Warning
!>
!> It is \b forbidden to modify turbulent viscosity \c visct here
!> (a specific subroutine is dedicated to that: \ref usvist)
!>
!> - icp = 1 must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable specific heat
!>    cpro_cp (otherwise: memory overwrite).
!>
!> - the kivisl field integer key (scalar_diffusivity_id)
!>    must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable viscosity
!>    \c viscls.
!>
!>
!> \remarks
!>  - This routine is called at the beginning of each time step
!>    Thus, <b> AT THE FIRST TIME STEP </b> (non-restart case), the only
!>    values initialized before this call are those defined
!>      - in the GUI or  \ref usipsu (cs_user_parameters.f90)
!>             - density    (initialized at \c ro0)
!>             - viscosity  (initialized at \c viscl0)
!>      - in the GUI or \ref cs_user_initialization
!>             - calculation variables (initialized at 0 by defaut
!>             or to the value given in the GUI or in \ref cs_user_initialization)
!>
!>  - We may define here variation laws for cell properties, for:
!>     - density:                                    rom    kg/m3
!>     - density at boundary faces:                  romb   kg/m3)
!>     - molecular viscosity:                        cpro_viscl  kg/(m s)
!>     - specific heat:                              cpro_cp     J/(kg degrees)
!>     - diffusivities associated with scalars:      cpro_vscalt kg/(m s)
!>
!> \b Warning: if the scalar is the temperature, cpro_vscalt corresponds
!> to its conductivity (Lambda) in W/(m K)
!>
!>
!> The types of boundary faces at the previous time step are available
!>   (except at the first time step, where arrays \c itypfb and \c itrifb have
!>   not been initialized yet)
!>
!> It is recommended to keep only the minimum necessary in this file
!>   (i.e. remove all unused example code)
!>
!>
!> \section usphyv_cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh
use lagran
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

integer          ivart, iel, ifac
integer          iscal, ifcvsl
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc
double precision xvart

double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:), pointer :: bfpro_rom, cpro_rom
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscalt, cpro_cp, cpro_beta
double precision, dimension(:), pointer :: cvar_scalt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0. Initializations to keep
!===============================================================================

!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!  Each example is bounded by a test using .false. as a precaution.
!  Replace .false. by .true to activate the example.

!  It is recommended to keep only the minimum necessary in this file
!  (i.e. remove all unused example code)


!  example 1: variable density as a function of temperature
!  example 2: variable viscosity as a function of tempeprature
!  example 3: variable specific heat as a function of tempeprature
!  example 4: variable Lambda/CP as a function of temperature
!             for temperature or enthalpy
!  example 5: variable sclalars diffusivity as a function of temperature
!===============================================================================


!===============================================================================
!  Example 1: variable density as a function of temperature
!  =========
!    Below, we define the same density law
!    Values of this property must be defined at cell centers
!      (and optionally, at boundary faces).
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       (and of its boundary conditions)
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Pointers to density values

  call field_get_val_s(icrom, cpro_rom)
  call field_get_val_s(ibrom, bfpro_rom)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  vara  = -4.0668d-3
  varb  = -5.0754d-2
  varc  =  1000.9d0

  ! Density at cell centers
  !------------------------
  ! law                    rho  = t  * ( a *  t +  b) +   c
  ! so      cpro_rom(iel) = xvart * (vara*xvart+varb) + varc

  ! Volumic thermal expansion coefficient
  !--------------------------------------
  ! law                     cpro_beta  = -1/rho * (d rho / d T)
  ! so cpro_beta(iel) = (-1.d0/cpro_rom(iel))*(2.d0*vara*xvart+varb)

  call field_get_val_s(iprpfl(ibeta), cpro_beta)

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_rom(iel) = xvart * (vara*xvart+varb) + varc
    cpro_beta(iel)= (-1.d0/cpro_rom(iel))*(2.d0*vara*xvart+varb)
  enddo


  ! Density at boundary faces
  !---------------------------

  ! By default, the value of rho at the boundary is the value taken
  !   at the center of adjacent cells. This is the recommended approach.
  ! To be in this case, nothing needs to be done:
  !   do not prescribe a value for bfpro_rom(ifac) and
  !   do not modify mbrom

  ! For users who do not wish to follow this recommendation, we
  !   note that the boundary temperature may be fictitious, simply
  !   defined so as to conserve a flux (this is especially the case
  !   at walls). The value of rho which is computed at the boundary
  !   when introducing this fictitious temperature in a physical law
  !   may thus be completely false (negative for example).

  ! If we wish to specify a law anyways:
  !                        rho  = t  * ( a *  t +  b) +   c
  ! so      bfpro_rom(ifac) = xvart * (vara*xvart+varb) + varc

  ! 't' being the temperature at boundary face centers, we may use the
  ! following lines of code (voluntarily deactived, as the must be used
  ! with caution):

  ! Note that when we prescribe the density at the boundary, it must be done
  ! at ALL boundary faces.
  !    ===

  if (.false.) then

    ! Boundary condition coefficients
    call field_get_coefa_s(ivarfl(ivart), coefap)
    call field_get_coefb_s(ivarfl(ivart), coefbp)

    ! Caution: mbrom = 1 is necessary for the law to be taken
    !                           into account.
    mbrom = 1

    do ifac = 1, nfabor

      ! ifabor(ifac) is the cell adjacent to the boundary face
      iel = ifabor(ifac)
      xvart = coefap(ifac) + cvar_scalt(iel)*coefbp(ifac)
      bfpro_rom(ifac) = xvart * (vara*xvart+varb) + varc
    enddo

  endif ! --- Test on .false.

endif ! --- Test on .false.


!===============================================================================
!  Example 2: variable viscosity as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit(1)
  endif

  ! --- Molecular dynamic viscosity

  call field_get_val_s(iprpfl(iviscl), cpro_viscl)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

  ! Molecular dynamic viscosity in kg/(m.s) at cell centers
  !--------------------------------------------------------
  ! law                    mu   = t * (t * (am * t + bm) + cm) + dm
  ! so      cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm
  enddo

endif ! --- Test on .false.


!===============================================================================
!  Example 3: specific heat as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Specific heat

  if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

  ! --- Stop if Cp is not variable

  if (icp.le.0) then
    write(nfecra,1000) icp
    call csexit (1)
  endif

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varac = 0.00001d0
  varbc = 1000.0d0

  ! Specific heat in J/(kg.degrees) at cell centers
  !------------------------------------------------
  ! law                    cpro_cp  = ac * t + bm
  ! so          cpro_cp(iel) = varac*xvart + varbc

  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_cp(iel) = varac*xvart + varbc
  enddo

endif ! --- Test on .false.


!======================================================================================
!  Example 4: Lambda/Cp a function of temperature for enthalpy or
!             Lambda    a function of temperature for temperature because Cp is put
!                       outside the divergence term
!  =========
!  Below, we define the same lambda/Cp ratio law (or lambda law if temperature is used)
!  Values of this property must be defined at cell centers
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Lambda/Cp of the thermal (or Lambda if temperature is used)

  call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)

  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_vscalt)
  else
    cpro_vscalt => NULL()
  endif

  ! --- Stop if Lambda/CP (or Lambda if temperature is used) is not variable

  if (ifcvsl.lt.0) then
    write(nfecra,1010) iscalt
    call csexit (1)
  endif

  ! if thermal variable is not temperature
  if (iscacp(iscal).le.0) then

    ! --- Specific heat

    if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)

    ! --- Coefficients of laws chosen by the user
    !       Values given here are fictitious

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! Lambda/Cp in kg/(m.s) at cell centers
    !--------------------------------------
    ! law    Lambda/Cp = {t * (t * (al * t +  bl) + cl) + dl} / Cp
    ! so     cpro_vscalt(iel) &
    !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)/cp0

    ! We assume Cp has been defined previously.

    if (icp.le.0) then

      ! --- If Cp is uniform, we use cp0
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) =                                           &
             (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)         &
             /cp0
      enddo

    else

      ! --- If Cp is not uniform, we use cpro_vscalt above
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) =                                           &
             (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)         &
             /cpro_cp(iel)
      enddo

    endif

  ! if iscalt is temperature, the Cp division is not needed
  else

    ! --- Coefficients of laws chosen by the user
    !       Values given here are fictitious

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! Lambda in W/(m.K) at cell centers
    !--------------------------------------
    ! law    Lambda = {t * (t * (al * t +  bl) + cl) + dl}
    ! so     cpro_vscalt(iel) &
    !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)

    do iel = 1, ncel
      xvart = cvar_scalt(iel)
      cpro_vscalt(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
    enddo

  endif

endif ! --- Test on .false.


!===============================================================================
!  Example 5: Diffusivity as a function of temperature for user scalars
!  =========
!    Excluding:
!      - temperature, enthalpy (handled above)
!      - fluctuation variances (property equal to that of the associated scalar)
!
!    Below, we define the same diffusivity law for all scalars (except the
!      ones excluded above).
!    Values of this property must be defined at cell centers
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  do iscal = 1, nscaus ! Loop on scalars

    ! --- If the variable is a fluctuation, its diffusivity is the same
    !       as that of the scalar to which it is attached:
    !       there is nothing to do here, we move on to the next variable
    !       without setting cpro_vscalt(iel).

    ! We only handle here variables which are not fluctuations
    if (iscavr(iscal).le.0) then

      ! --- Number of the thermal variable
      !     To use user scalar 2 instead, write 'ivart = isca(2)'

      if (iscalt.gt.0) then
        ivart = isca(iscalt)
        call field_get_val_s(ivarfl(ivart), cvar_scalt)
      else
        write(nfecra,9010) iscalt
        call csexit (1)
      endif

      ! --- Scalar's Lambda

      call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0) then
        call field_get_val_s(ifcvsl, cpro_vscalt)
      else
        cpro_vscalt => NULL()
      endif

      ! --- Stop if Lambda is not variable

      if (ifcvsl.lt.0) then
        write(nfecra,1010) iscal
        call csexit (1)
      endif

      ! --- Coefficients of laws chosen by the user
      !     Values given here are fictitious

      varal = -3.3283d-7
      varbl =  3.6021d-5
      varcl =  1.2527d-4
      vardl =  0.58923d0

      ! Lambda in kg/(m.s) at cell centers
      !--------------------------------------
      ! law    Lambda = {t * (t * (al * t +  bl) + cl) + dl}
      ! so     cpro_vscalt(iel) &
      !             = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)

      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
      enddo

    endif ! --- Tests on 'iscavr'

  enddo ! --- Loop on scalars
endif ! --- Test on .false.

!===============================================================================
!  Example 6: Diffusivity as a function of temperature for user scalars
!  =========
!    Excluding:
!      - temperature, enthalpy (handled above)
!      - fluctuation variances (property equal to that of the associated scalar)
!
!    Below, we define the same diffusivity law for all scalars (except the
!      ones excluded above).
!    Values of this property must be defined at cell centers
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  call field_get_val_s(iprpfl(iviscl), cpro_viscl)
  call field_get_val_s(icrom, cpro_rom)

  do iscal = 1, nscaus ! Loop on scalars

    ! --- If the variable is a fluctuation, its diffusivity is the same
    !       as that of the scalar to which it is attached:
    !       there is nothing to do here, we move on to the next variable
    !       without setting cpro_vscalt(iel).

    ! We only handle here variables which are not fluctuations
    if (iscavr(iscal).le.0) then

      ! --- Number of the thermal variable
      !     To use user scalar 2 instead, write 'ivart = isca(2)'

      if (iscalt.gt.0) then
        ivart = isca(iscalt)
        call field_get_val_s(ivarfl(ivart), cvar_scalt)
      else
        write(nfecra,9010) iscalt
        call csexit (1)
      endif

      ! --- Scalar's Lambda

      call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0) then
        call field_get_val_s(ifcvsl, cpro_vscalt)
      else
        cpro_vscalt => NULL()
      endif

      ! --- Stop if Lambda is not variable

      if (ifcvsl.lt.0) then
        write(nfecra,1010) iscal
        call csexit (1)
      endif

      ! Stokes-Einstein equation

      varal =  1.38d-23 ! Boltzmann constant
      varbl =  0.3d-9   ! Diameter of soluble ions

      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_vscalt(iel) = (cpro_rom(iel) * varal * xvart)       &
                           /(3.d0 * pi * cpro_viscl(iel) * varbl)
      enddo

    endif ! --- Tests on 'iscavr'

  enddo ! --- Loop on scalars
endif ! --- Test on .false.

!===============================================================================
!  Example 7: Solubility as a function of temperature for user scalars
!  =========
!
!    Excluding:
!      - temperature, enthalpy (handled above)
!      - fluctuation variances (property equal to that of the associated scalar)
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then
  if (ipreci == 1) then
    do iel = 1, ncel
      solub(iel) = -0.0088d0 *  cvar_scalt(iel) + 3.8839d0 ! ppb (ug/L)
      solub(iel) = solub(iel)*1.d-9 * cpro_rom(iel) ! kg/m3
   enddo
  endif
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@      usipsu indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ', i10                                  ,/,&
'@      la diffusivite est uniforme alors que                 ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      the diffusivity is uniform while',/,                      &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@',/,                                                            &
'@    The variable on which physical properties depend does',/,   &
'@      seem to be a calculation variable.',/,                    &
'@    Indeed, we are trying to use the temperature while',/,      &
'@      iscalt = ',i10                                  ,/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Check the programming in usphyv (and the test when',/,      &
'@      defining ivart).',/,                                      &
'@    Check the definition of calculation variables in',/,        &
'@      usipsu. If a scalar should represent the,',/,             &
'@      temperature, check that iscalt has been defined',/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

#endif

!----
! End
!----

return
end subroutine usphyv


!===============================================================================
!===============================================================================
! Purpose:
! -------

!> usvist
!> \brief Modify turbulent viscosity
!>
!> This subroutine is called at beginning of each time step
!> after the computation of the turbulent viscosity
!> (physical quantities have already been computed in \ref usphyv).
!>
!> Turbulent viscosity \f$ \mu_T \f$ (kg/(m s)) can be modified.
!>
!> A modification of the turbulent viscosity can lead to very
!> significant differences betwwen solutions and even give wrong
!> results.
!>
!> This subroutine is therefore reserved to expert users.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        head loss cell numbering
!> \param[in]     icetsm        numbering of cells with mass source term
!> \param[in]     itypsm        kind of mass source for each variable
!>                               (cf. \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head loss terms
!> \param[in]     smacel        values of variables related to mass source
!>                              term. If ivar=ipr, smacel=mass flux
!_______________________________________________________________________________

subroutine usvist &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh
use field
use field_operator

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

integer          iel, inc, iprev
double precision dudx, dudy, dudz, sqdu, visct, rom

double precision, allocatable, dimension(:,:,:) :: gradv
double precision, dimension(:), pointer :: cpro_rom
double precision, dimension(:), pointer :: cpro_visct

!===============================================================================

!===============================================================================
! 1.  Example :
!                visct = max(visct, rom * sqrt(dudx**2 + dudy**2 + dudz**2)
!                (intentionally fancyful relation)
!                Remark: incomming viscosity is consistent with the selected
!                turbulence modelling

!===============================================================================

!  The test below allows deactivating instructions (which are defined
!     only as a starting example)
!  Replace .true. with .false. or remove this test to activate the example.

if (.true.) return

!=============================================================================
! 1.2 Initialization
!=============================================================================

! Allocate work arrays
! First component is for x,y,z  and the 2nd for u,v,w
allocate(gradv(3,3,ncelet))

call field_get_val_s(iprpfl(ivisct), cpro_visct)
call field_get_val_s(icrom, cpro_rom)

!===============================================================================
! 1.3 Compute velocity gradient
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,     &
                           gradv)

!===============================================================================
! 1.4 Computation of the dynamic viscosity
!===============================================================================

do iel = 1, ncel

  ! --- Current dynamic viscosity and fluid density
  visct = cpro_visct(iel)
  rom   = cpro_rom(iel)
  ! --- Various computations
  dudx = gradv(1,1,iel)
  dudy = gradv(1,2,iel)
  dudz = gradv(1,3,iel)
  sqdu = sqrt(dudx**2+dudy**2+dudz**2)

  ! --- Computation of the new dynamic viscosity
  visct = max (visct,rom*sqdu)

  ! --- Store the new computed dynamic viscosity
  cpro_visct(iel) = visct

enddo

! Free memory
deallocate(gradv)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine usvist


!===============================================================================


subroutine ussmag &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   smagor , mijlij , mijmij )

!===============================================================================
! FONCTION :
! --------

! MODIFICATION UTILISATEUR DE LA CONSTANTE DE SMAGORINSKY
! DANS LE CAS DE L'UTILISATION D'UN MODELE DYNAMIQUE

!              SMAGOR = Mij.Lij / Mij.Mij

! EN FAIT, DES MOYENNES LOCALES DU NUMERATEUR ET DU DENOMINATEUR
! SONT REALISEES AVANT L'APPEL A USSMAG, SOIT

!              SMAGOR = < Mij.Lij > / < Mij.Mij >

! DANS CET ROUTINE, Mij.Lij ET Mij.Mij SONT PASSES EN ARGUMENT
! AVANT LA MOYENNE LOCALE.
! DANS L'EXEMPLE CI-DESSOUS ON REALISE UNE MOYENNE LOCALE DU
! RAPPORT.
!              SMAGOR = < Mij.Lij / Mij.Mij >

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. cs_user_mass_source_terms)     !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smagor(ncelet)   ! tr ! <-- ! constante de smagorinsky dans le cas           !
!                  !    !     ! d'un modlele dynamique                         !
! mijlij(ncelet    ! tr ! <-- ! mij.lij avant moyenne locale                   !
! mijmij(ncelet    ! tr ! <-- ! mij.mij avant moyenne locale                   !
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
use cstnum
use optcal
use cstphy
use entsor
use parall
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smagor(ncelet), mijlij(ncelet), mijmij(ncelet)

! Local variables

integer          iel

double precision, allocatable, dimension(:) :: w1, w2, w3

!===============================================================================

!  The test below allows deactivating instructions (which are defined
!     only as a starting example)
!  Replace .true. with .false. or remove this test to activate the example.

if (.true.) return

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

!===============================================================================
! 2.  MOYENNE SPATIALE SUR LE VOISINAGE ETENDU

!     Dans le cas ou l'utilisateur souhaite utilise le voisinage
!       etendu, il est fortement conseille de passer le mode de
!       calcul des gradients en IMRGRA = 2, afin de conserver
!       la totalite du voisinage etendu. En effet, le calcul de
!       moyenne locale est generalement degradee en voisinage reduit
!       (IMRGRA = 3).

!===============================================================================

!     On calcule le rapport
do iel = 1, ncel
  if (abs(mijmij(iel)).le.epzero) then
    w1(iel) = smagmx**2
  else
    w1(iel) = mijlij(iel)/mijmij(iel)
  endif
enddo

!     On passe dans le filtre local
call cfiltr ( w1     , smagor , w2     , w3     )
!==========

! Free memory
deallocate(w1, w2, w3)

!----
! End
!----

return
end subroutine ussmag


!===============================================================================

!===============================================================================
! Purpose:
! -------

!> usvima
!> \brief User subroutine dedicated the use of ALE
!>  (Arbitrary Lagrangian Eulerian Method): fills mesh viscosity arrays.
!>
!> This subroutine is called at the beginning of each time step.
!>
!> Here one can modify mesh viscosity value to prevent cells and nodes
!> from huge displacements in awkward areas, such as boundary layer for example.
!>
!> IF variable IORTVM = 0, mesh viscosity modeling is isotropic therefore VISCMX
!> array only needs to be filled.
!> IF variable IORTVM = 1, mesh viscosity modeling is orthotropic therefore
!> all arrays VISCMX, VISCMY and VISCMZ need to be filled.
!>
!> Note that VISCMX, VISCMY and VISCMZ arrays are initialized at the first time step
!> to the value of 1.
!>
!> \section usvima_cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[out]    viscmx        mesh viscosity in X direction
!> \param[out]    viscmy        mesh viscosity in Y direction
!> \param[out]    viscmz        mesh viscosity in Z direction
!_______________________________________________________________________________

subroutine usvima &
 ( nvar   , nscal  ,                                              &
   dt     ,                                                       &
   viscmx , viscmy , viscmz )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)
double precision viscmx(ncelet), viscmy(ncelet), viscmz(ncelet)

! Local variables

integer          iel
double precision rad, xr2, xcen, ycen, zcen

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
!  1. Example :
!       One gives a huge mesh viscosity value to the cells included in a circle
!       with a radius 'rad' and a center '(xcen, ycen, zcen)'

!     In general it appears quite much easier to fill mesh viscosity arrays at
!     the beginning of the calculations basing on the initial geometry.

!  The test on .false. allows deactivating instructions (which are defined
!     only as a starting example)

if (.false.) then

  if (ntcabs.eq.0) then
    rad = (1.d-3)**2
    xcen  = 1.d0
    ycen  = 0.d0
    zcen  = 0.d0

    do iel = 1, ncel
      xr2 = (xyzcen(1,iel)-xcen)**2 + (xyzcen(2,iel)-ycen)**2       &
           + (xyzcen(3,iel)-zcen)**2
      if (xr2.lt.rad) viscmx(iel) = 1.d10
    enddo

    ! 2. In case of orthotropic mesh viscosity modeling one can choose
    !    to submit nodes to a lower stress in Z direction

    if (iortvm.eq.1) then
      do iel = 1, ncel
        viscmy(iel) = viscmx(iel)
        viscmz(iel) = 1.d0
      enddo
    endif

  endif

endif ! --- Test on .false.

!----
! Formats
!----

!----
! End
!----

return
end subroutine usvima

!===============================================================================
! Purpose:
! -------

!> usatph
!> \brief User subroutine dedicated to modifie physical properties of the
!>        atmospheric module
!>
!> This subroutine is called at beginning of each time step at the end of
!> atphyv.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine usatph

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use field
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!----
! Formats
!----

!----
! End
!----

return
end subroutine usatph
