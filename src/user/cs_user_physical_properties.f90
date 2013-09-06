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
!>    in \ref usipph if we wish to define a variable specific heat
!>    cp (otherwise: memory overwrite).
!>
!> - ivisls = 1 must <b> have been specified </b>
!>    in \ref usipsc if we wish to define a variable viscosity
!>    \c viscls (otherwise: memory overwrite).
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
!>     - molecular viscosity:                        viscl  kg/(m s)
!>     - specific heat:                              cp     J/(kg degrees)
!>     - diffusivities associated with scalars:      visls kg/(m s)
!>
!> \b Warning: if the scalar is the temperature, visls corresponds
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
!> \section cell_id Cells identification
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
!> \param[in]     ibrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!> \param[in]                    (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   ibrom  ,                                                       &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          ibrom

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)

! Local variables

integer          ivart, iclvar, iel, ifac
integer          ipcrom, ipbrom, ipcvis, ipccp
integer          ipcvsl, ith, iscal, ii
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc
double precision xrtp

double precision, dimension(:), pointer :: coefap, coefbp

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
!  ===================================================================

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
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Position of boundary conditions for variable 'ivart'

  iclvar = iclrtp(ivart,icoef)

  ! --- Rank of density
  !     in 'propce', physical properties at element centers:       'ipcrom'
  !     in 'propfb', physical properties at boundary face centers: 'ipbrom'

  ipcrom = ipproc(irom)
  ipbrom = ipprob(irom)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  vara  = -4.0668d-3
  varb  = -5.0754d-2
  varc  =  1000.9d0

  ! Density at cell centers
  !------------------------
  ! law                    rho  = t  * ( a *  t +  b) +   c
  ! so      propce(iel, ipcrom) = xrtp * (vara*xrtp+varb) + varc

  ! Volumic thermal expansion coefficient
  !--------------------------------------
  ! law                     beta  = -1/rho * (d rho / d T)
  ! so propce(iel, ipproc(ibeta)) = (-1.d0/propce(iel,ipcrom))*(2.d0*vara*xrtp+varb)

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcrom) = xrtp * (vara*xrtp+varb) + varc
    propce(iel,ipproc(ibeta))= (-1.d0/propce(iel,ipcrom))*(2.d0*vara*xrtp+varb)
  enddo


  ! Density at boundary faces
  !---------------------------

  ! By default, the value of rho at the boundary is the value taken
  !   at the center of adjacent cells. This is the recommended approach.
  ! To be in this case, nothing needs to be done:
  !   do not prescribe a value for propfb(ifac, ipbrom) and
  !   do not modify ibrom

  ! For users who do not wish to follow this recommendation, we
  !   note that the boundary temperature may be fictitious, simply
  !   defined so as to conserve a flux (this is especially the case
  !   at walls). The value of rho which is computed at the boundary
  !   when introducing this fictitious temperature in a physical law
  !   may thus be completely false (negative for example).

  ! If we wish to specify a law anyways:
  !                        rho  = t  * ( a *  t +  b) +   c
  ! so      propfb(iel, ipbrom) = xrtp * (vara*xrtp+varb) + varc

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

    ! Caution: ibrom = 1 is necessary for the law to be taken
    !                           into account.
    ibrom = 1

    do ifac = 1, nfabor

      ! ifabor(ifac) is the cell adjacent to the boundary face
      iel = ifabor(ifac)
      xrtp = coefap(ifac) + rtp(iel, ivart)*coefbp(ifac)
      propfb(ifac, ipbrom) = xrtp * (vara*xrtp+varb) + varc
    enddo

  endif ! --- Test on .false.

endif ! --- Test on .false.


!===============================================================================
!  Example 2: variable viscosity as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
  else
    write(nfecra,9010) iscalt
    call csexit(1)
  endif

  ! --- Rank of molecular dynamic viscosity
  !     in 'propce', physical properties at element centers: 'ipcvis'

  ipcvis = ipproc(iviscl)

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

  ! Molecular dynamic viscosity in kg/(m.s) at cell centers
  !--------------------------------------------------------
  ! law                    mu   = t * (t * (am * t + bm) + cm) + dm
  ! so      propce(iel, ipcvis) = xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvis) =                                        &
         xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
  enddo

endif ! --- Test on .false.


!===============================================================================
!  Example 3: specific heat as a function of temperature
!  =========
!    Below, we define the same viscosity law
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Rank of the specific heat
  !     in 'propce', physical properties at element centers: 'ipccp'

  if (icp.gt.0) then
    ipccp  = ipproc(icp   )
  else
    ipccp  = 0
  endif

  ! --- Stop if Cp is not variable

  if (ipccp.le.0) then
    write(nfecra,1000) icp
    call csexit (1)
  endif

  ! --- Coefficients of laws chosen by the user
  !       Values given here are fictitious

  varac = 0.00001d0
  varbc = 1000.0d0

  ! Specific heat in J/(kg.degrees) at cell centers
  !------------------------------------------------
  ! law                    cp  = ac * t + bm
  ! so      propce(iel, ipccp) = varac*xrtp + varbc

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipccp ) = varac*xrtp + varbc
  enddo

endif ! --- Test on .false.


!======================================================================================
!  Example 4: Lambda/Cp a function of temperature for enthalpy or
!             Lambda    a function of temperature for temperature because Cp is put
!                       outside the divergence term
!  =========
!  Below, we define the same lambda/Cp ratio law (or lambda law if temperature is used)
!  Values of this property must be defined at cell centers
!  ====================================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the thermal variable
  !       To use user scalar 2 instead, write 'ivart = isca(2)'

  if (iscalt.gt.0) then
    ivart = isca(iscalt)
  else
    write(nfecra,9010) iscalt
    call csexit (1)
  endif

  ! --- Rank of Lambda/Cp of the thermal (or Lambda if temperature is used)
  !     in 'propce', physical properties at element centers: 'ipcvsl'

  if (ivisls(iscalt).gt.0) then
    ipcvsl = ipproc(ivisls(iscalt))
  else
    ipcvsl = 0
  endif

  ! --- Stop if Lambda/CP (or Lambda if temperature is used) is not variable

  if (ipcvsl.le.0) then
    write(nfecra,1010)                                          &
         iscalt, iscalt, ivisls(iscalt)
    call csexit (1)
  endif

  ! if iscalt is not temperature
  if (abs(iscsth(iscalt)).ne.1) then

    ! --- Rank of the specific heat
    !     in 'propce', physical properties at element centers: 'ipccp'

    if (icp.gt.0) then
      ipccp  = ipproc(icp   )
    else
      ipccp  = 0
    endif

    ! --- Coefficients of laws chosen by the user
    !       Values given here are fictitious

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! Lambda/Cp in kg/(m.s) at cell centers
    !--------------------------------------
    ! law    Lambda/Cp = {t * (t * (al * t +  bl) + cl) + dl} / Cp
    ! so     propce(iel,ipcvsl) &
    !             = (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)/cp0

    ! We assume Cp has been defined previously.

    if (ipccp.le.0) then

      ! --- If Cp is uniform, we use cp0
      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)         &
             /cp0
      enddo

    else

      ! --- If Cp is not uniform, we use propce above
      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)         &
             /propce(iel,ipccp)
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
    ! so     propce(iel,ipcvsl) &
    !             = (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)

    do iel = 1, ncel
      xrtp = rtp(iel,ivart)
      propce(iel,ipcvsl) =                                      &
           (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
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
!  ===================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  do ii = 1, nscaus ! Loop on scalars

    ! --- Number of user scalar 'ii' in the lsit of scalars
    iscal = ii

    ! --- If it is a thermal variable, it has already been handled above
    ith = 0
    if (iscal.eq.iscalt) ith = 1

    ! --- If the variable is a fluctuation, its diffusivity is the same
    !       as that of the scalar to which it is attached:
    !       there is nothing to do here, we move on to the next variable
    !       without setting propce(iel,ipcvsl).

    ! We only handle here non-thermal variables which are not fluctuations
    if (ith.eq.0.and.iscavr(iscal).le.0) then

      ! Position of variables, coefficients
      ! -----------------------------------

      ! --- Number of the thermal variable
      !       To use user scalar 2 instead, write 'ivart = isca(2)'

      if (iscalt.gt.0) then
        ivart = isca(iscalt)
      else
        write(nfecra,9010) iscalt
        call csexit (1)
      endif

      ! --- Rank of scalar's Lambda
      !     in 'propce', physical properties at element centers: 'ipcvsl'

      if (ivisls(iscal).gt.0) then
        ipcvsl = ipproc(ivisls(iscal))
      else
        ipcvsl = 0
      endif

      ! --- Stop if Lambda is not variable

      if (ipcvsl.le.0) then
        write(nfecra,1010) iscal, iscal, ivisls(iscal)
        call csexit (1)
      endif

      ! --- Coefficients of laws chosen by the user
      !       Values given here are fictitious

      varal = -3.3283d-7
      varbl =  3.6021d-5
      varcl =  1.2527d-4
      vardl =  0.58923d0

      ! Lambda in kg/(m.s) at cell centers
      !--------------------------------------
      ! law    Lambda = {t * (t * (al * t +  bl) + cl) + dl}
      ! so     propce(iel,ipcvsl) &
      !             = (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)

      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
      enddo

    endif ! --- Tests on 'ith' and 'iscavr'

  enddo ! --- Loop on scalars
endif ! --- Test on .false.


!===============================================================================

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
'@      usipph indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipph ou usphyv.                              ',/,&
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
'@    Pour le scalaire ',I10                                   ,/,&
'@      usipsc indique que la diffusivite est uniforme        ',/,&
'@        IVISLS(',I10   ,') = ',I10   ,' alors que           ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsc ou usphyv.                              ',/,&
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
'@      usipph specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipph or usphyv.',/,                                &
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
'@      usipsu specifies that the diffusivity is uniform',/,      &
'@        ivislc(',i10   ,') = ',i10   ,' while',/,               &
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

!> uscfpv
!> \brief  Set (variable) physical properties for the compressible flow scheme.
!>
!> \section des Description
!>
!> This subroutine replaces the user subroutine \ref usphyv for the
!> compressible flow scheme.
!>
!> This subroutine is called at the beginning of each time step.
!>
!> At the very first time step (not at restart), the only variables that
!> have been initialized are those provided:
!>   - in the GUI and in the user subroutines \ref usipsu and \ref uscfx2; ex.:
!>     - the density             (set to ro0)
!>     - the molecular viscosity (set to viscl0)
!>     - the volumetric molecular viscosity (set to viscv0)
!>     - the molecular thermal conductivity (set to visls0(itempk))
!>   - in the user subroutine \ref cs_user_initialization; ex.:
!>     - the unknown variables (null by default)
!>
!> This subroutine allows the user to set the cell values for:
!>   - the molecular viscosity:                            viscl  kg/(m s)
!>   - the isobaric specific heat
!>   (\f$ C_p = \left. \dfrac{\dd h}{\dd T}\right|_P \f$): cp     J/(kg degree)
!>   - the molecular thermal conductivity:                 lambda W/(m degree)
!>   - the molecular diffusivity for user-defined scalars: viscls kg/(m s)
!>
!> \section Warnings
!>
!> The density <b> must not </b> be set here: for the compressible scheme,
!> it is one of the unknowns, and it can be initialized as such in the user
!> subroutine \ref cs_user_initialization (rtp array).
!>
!> The turbulent viscosity <b> must not </b> be modified here (to modify this
!> variable, use the user subroutine \ref usvist)
!>
!> To set a variable isobaric specific heat, the integer \c icp must
!> have been set to 1: the value for \c icp is set automatically in the
!> subroutine \ref cfther, depending on the thermodynamics laws selected
!> by the user.
!>
!> To set a variable diffusivity for a given user-defined scalar, the
!> variable \c ivisls(scalar_number) must have been set to 1 in the user
!> subroutine \ref usipsc or in the GUI (otherwise, a memory problem is
!> expected).
!>
!> Examples are provided in the present subroutine (but they do not have
!> any physical signification).
!>
!> \section cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!>
!> The type of the boundary faces at the previous time step is available
!> (except at the first time step, since the arrays \c itypfb and \c itrifb have
!> not yet been set);
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
!> \param[in]     rtpR, rtpa    calculated variables at cell centers
!>                               (at current and preceding time steps)
!> \param[in,out] propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!_______________________________________________________________________________

subroutine uscfpv &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use period
use ppppar
use ppthch
use ppincl
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)

! Local variables

integer          ivart, iel
integer          ipcvis, ipcvsv, ipccp
integer          ipcvsl, ith, iscal, ii, iccfth, imodif
double precision varam, varbm, varcm, vardm
double precision varal, varbl, varcl, vardl
double precision varac, varbc
double precision xrtp

double precision rvoid(1)
double precision, allocatable, dimension(:) :: w1, w2, w3

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     However, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (1.eq.1) then
  return
endif


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Mandatory initializations
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

!===============================================================================

! Warning: the examples provided below are physically meaningless.
! =======

! These examples must be adapted by the user. Hence, the default
! (library reference) version stops immediately after each example
! (the 'call csexit(1)' directive must be discarded to use the
! portion of code).

! It is adviced to discard all the examples that are not necessary, so
! as to minimize the risk of error.

! List of examples
! ================

! Ex. 1: molecular viscosity varying with temperature
! Ex. 2: molecular volumetric viscosity varying with temperature
! Ex. 3: isobaric specific heat varying with temperature
! Ex. 4: molecular thermal conductivity varying with temperature
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature

!===============================================================================


!===============================================================================
! Ex. 1: molecular viscosity varying with temperature
! =====
!    The values of the molecular viscosity are provided as a function of
!    the temperature. All variables are evaluated at the cell centres.
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! --- Rank of the temperature in the array 'rtp'
  !     To refer to the user-defined scalar number 2 instead, for example, use
  !     ivart = isca(2)

  ivart = isca(itempk)

  ! --- Rank 'ipcvis' of the molecular dynamic viscosity
  !     in the array 'propce' (physical properties at the cell centers)

  ipcvis = ipproc(iviscl)

  ! --- User-defined coefficients for the selected law.
  !     The values hereafter are provided as a mere example. They
  !     are physically meaningless.

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

  ! --- Molecular dynamic viscosity mu at the cell centres, kg/(m s)
  !     In this example, mu is provided as a function of the temperature T:
  !       mu(T)              =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
  !     that is:
  !       propce(iel,ipcvis) =   xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvis) = xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
  enddo

endif ! --- Test on .false.

! --- Discard the following test so that the code does not stop
if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 2: molecular volumetric viscosity varying with temperature
! =====
!    The values of the molecular volumetric viscosity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! --- Rank of the temperature in the array 'rtp'
  !     To refer to the user-defined scalar number 2 instead, for example, use
  !     ivart = isca(2)

  ivart = isca(itempk)

  ! --- Rank 'ipcvsv' of the molecular dynamic viscosity
  !     in the array 'propce' (physical properties at the cell centers)

  if (iviscv.gt.0) then
    ipcvsv = ipproc(iviscv)
  else
    ipcvsv = 0
  endif

  ! --- Stop if the viscosity has not been defined as variable

  if (ipcvsv.le.0) then
    write(nfecra,2000) iviscv
    call csexit (1)
  endif

  ! --- User-defined coefficients for the selected law.
  !     The values provided hereafter are provided as a mere example. They
  !     are physically meaningless.

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

  ! --- Molecular dynamic volumetric viscosity kappa at the cell centres, kg/(m s)
  !     In this example, kappa is provided as a function of the temperature T:
  !       kappa(T)           =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
  !     that is:
  !       propce(iel,ipcvsv) =   xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvsv) = xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
  enddo

endif ! --- Test on .false.

! --- Discard the following test so that the code do not stop
if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 3: isobaric specific heat varying with temperature
! =====
!    The values of the isobaric specific heat values are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! Warning:
  ! =======
  ! do not discard the call to the subroutine 'usthht' at the end of this
  ! example: its purpose is to calculate the isochoric specific heat.
  ! Indeed, this variable needs to be computed from the isobaric specific heat
  ! using the thermodynamics laws.

  ! --- Rank of the temperature in the array 'rtp'
  !     To refer to the user-defined scalar number 2 instead, for example, use
  !     ivart = isca(2)

  ivart = isca(itempk)

  ! --- Rank 'ipcpp' of the isobaric specific heat
  !     in the array 'propce' (physical properties at the cell
  !     centers)

  if (icp.gt.0) then
    ipccp  = ipproc(icp   )
  else
    ipccp  = 0
  endif

  ! --- Stop if the iobaric or iochoric specific heat (cp or cv) has not
  !     been defined as variable

  if (ipccp.le.0) then
    write(nfecra,1000) icp
    call csexit (1)
  endif
  if (icv.le.0) then
    write(nfecra,1001) icv
    call csexit (1)
  endif

  ! --- User-defined coefficients for the selected law.
  !     The values provided hereafter are provided as a mere example. They
  !     are physically meaningless.

  varac = 0.00001d0
  varbc = 1000.0d0

  ! --- Isobaric specific heat cp at the cell centres, J/(kg degree)
  !     In this example, cp is provided as a function of the temperature T:
  !       cp(T)              =      ac * T  + ab
  !     that is:
  !       propce(iel,ipccp ) =    varac*xrtp+varbc

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipccp ) = varac*xrtp + varbc
  enddo

  ! --- The isochoric specific heat is deduced from the isobaric specific
  !     heat using the subroutine 'cfther'.

  iccfth = 432
  imodif = 0

  call cfther                                                       &
  !==========
   ( nvar   ,                                                       &
     iccfth , imodif ,                                              &
     dt     , rtp    , rtpa  , propce , propfb ,                    &
     propce(1, ipproc(icv))  , w1     , w2     , w3     ,           &
     rvoid  , rvoid )

endif ! --- Test on .false.

! --- Discard the following test so that the code do not stop
if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 4: molecular thermal conductivity varying with temperature
! =====
!    The values of the molecular thermal conductivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! --- Rank of the temperature in the array 'rtp'
  !     To refer to the user-defined scalar number 2 instead, for example, use
  !     ivart = isca(2)

  ivart = isca(itempk)

  ! --- Rank 'ipcvsl' of the olecular thermal conductivity
  !     in the array 'propce' (physical properties at the cell
  !     centers)

  if (ivisls(itempk).gt.0) then
    ipcvsl = ipproc(ivisls(itempk))
  else
    ipcvsl = 0
  endif

  ! --- Stop if the molecular thermal conductivity has not
  !     been defined as variable

  if (ipcvsl.le.0) then
    write(nfecra,1010) itempk, itempk, ivisls(itempk)
    call csexit (1)
  endif

  ! --- User-defined coefficients for the selected law.
  !     The values provided hereafter are provided as a mere example. They
  !     are physically meaningless.

  varal = -3.3283d-7
  varbl =  3.6021d-5
  varcl =  1.2527d-4
  vardl =  0.58923d0

  ! --- Molecular thermal conductivity lambda at the cell centres, W/(m degree)
  !     In this example, lambda is provided as a function of the temperature T:
  !       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
  !     that is:
  !       propce(iel,ipcvsl) =   xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvsl) = (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
  enddo

endif ! --- Test on .false.

! --- Discard the following test so that the code do not stop
if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature
! =====
!    The molecular diffusivity can be set for all the user-defined scalars
!    ** except **:
!      - temperature and enthalpy (already dealt with above: for these
!        variables, the 'diffusivity' is the thermal conductivity)
!      - variances of the fluctuations of another scalar variable (the
!        diffusivity is assumed to be equal to that of the associated
!        scalar)
!    The values of the molecular diffusivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

!    The test on .false. allows deactivating instructions (which are defined
!       only as a starting example)

if (.false.) then

  ! --- Loop on the scalars
  do ii = 1, nscaus

    ! --- Rank of the ii-th scalar in the list of all scalars
    iscal = ii


    ! --- If the scalar is the temperature, it is marked by ith = 1
    !     so that it will be skipped.

    ith = 0
    if (iscal.eq.itempk) ith = 1

    ! --- If the variable represents the variance of the fluctuations of
    !     another scalar variable (iscavr <= 0), it is simply skipped.

    if (ith.eq.0.and.iscavr(iscal).le.0) then

      ! --- Here, iscal points to any scalar variable except the temperature,
      !     the enthalpy and the variance of the fluctuations of another
      !     scalar variable.

      ! --- Rank of the temperature in the array 'rtp'
      !     To refer to the user-defined scalar number 2 instead, for example, use
      !     ivart = isca(2)

      ivart = isca(itempk)

      ! --- Rank 'ipcvsl' of the molecular diffusivity of the current scalar iscal
      !     in the array 'propce' (physical properties at the cell centers)

      if (ivisls(iscal).gt.0) then
        ipcvsl = ipproc(ivisls(iscal))
      else
        ipcvsl = 0
      endif

      ! --- Stop if the molecular diffusivity has not been defined as variable

      if (ipcvsl.le.0) then
        write(nfecra,1010) iscal, iscal, ivisls(iscal)
        call csexit (1)
      endif

      ! --- User-defined coefficients for the selected law.
      !     The values provided hereafter are provided as a mere example. They
      !     are physically meaningless.

      varal = -3.3283d-7
      varbl =  3.6021d-5
      varcl =  1.2527d-4
      vardl =  0.58923d0

      ! --- Molecular diffusivity lambda at the cell centres, kg/(m s)
      !     In this example, lambda is provided as a function of the temperature T:
      !       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
      !     that is:
      !       propce(iel,ipcvsl) =   xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl

      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) = (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
      enddo


    endif
    ! --- End of the tests on ith and iscavr

  enddo
  ! --- End of the loop on the scalars

endif ! --- Test on .false.

! Free memory
deallocate(w1, w2, w3)

! --- Discard the following test so that the code do not stop
if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usipph'', the',/, &
'@         isobaric specific heat is declared as a property',/,   &
'@         uniform in space: icp = ',i10   ,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipph'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1001 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         isochoric specific heat is declared as a property',/,  &
'@         uniform in space: icv = ',i10   ,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@     For the scalar ',i10,/,                                    &
'@       in the GUI or in the user subroutine ''usipsc'', the',/, &
'@         molecular diffusivity is declared as a property',/,    &
'@         uniform in space: ivisls(',i10   ,') = ',i10   ,/,     &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsc'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the user subroutine ''uscfx2'', the molecular',/,     &
'@         volumetric viscosity is declared as a property',/,     &
'@         uniform in space: iviscv = ',i10,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the user subroutines',/,      &
'@    ''uscfx2'' or ''uscfpv''.',/,                               &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     Call to ''csexit'' from the user subroutine ''uscfpv''.',/,&
'@',/,                                                            &
'@     The subroutine ''csexit'' (run stop) was called from ',/,  &
'@       within the user subroutine ''uscfpv''. The user shall',/,&
'@       ensure that all the default examples provided in the',/, &
'@       reference version of the user subroutine have been',/,   &
'@       discarded. It shall also be checked that there is no',/, &
'@       remaining stopping test at the end of the examples ',/,  &
'@       that have been retained.',/,                             &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Check and modify the user subroutine ''uscfpv''.',/,          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine uscfpv


!===============================================================================


subroutine uselph &
!================

 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb )

!===============================================================================
! FONCTION :
! --------

!   REMPLISSAGE DES VARIABLES PHYSIQUES POUR LE MODULE ELECTRIQUE

!     ----> Effet Joule
!     ----> Arc Electrique
!     ----> Conduction Ionique

!      1) Masse Volumique
!      2) Viscosite moleculaire
!      3) Chaleur massique Cp
!      4) Lambda/Cp moleculaire
!      4) Diffusivite moleculaire



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)

! Pour le module electrique, toutes les proprietes physiques sont
!   supposees variables et contenues dans le tableau PROPCE
!   (meme si elles sont physiquement constantes)


! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usiniv/useliv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a l
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
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
use optcal
use cstphy
use entsor
use ppppar
use ppthch
use ppincl
use elincl
use field
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          ibrom
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)

! Local variables

integer          iel
integer          ipcrom, ipcvis, ipccp , ipcvsl, ipcsig
integer          mode

double precision tp
double precision xkr   , xbr
double precision rom0  , temp0 , dilar , aa    , bb    , cc
double precision srrom1, rhonp1

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 0 - INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire


ipass = ipass + 1

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================


!     En Joule, on s'arrete : il faut que l'utilisateur
!       donne les proprietes physiques
if ( ippmod(ieljou).ge.1 ) then

  write(nfecra,9010)
  call csexit (1)

!     En Arc on continue car on a un fichier de donnees
!       Un message indique que l'utilisateur n'a rien fourni
elseif (ippmod(ielarc).ge.1) then

  if (ipass.eq.1) then
    write(nfecra,9011)
  endif

  return

endif

 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES PROP. PHYSIQUES   ',/,&
'@    =========                                               ',/,&
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uselph DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       proprietes physiques. Il est indispensable.          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(/,                                                   &
' Module arc electrique: pas d''intervention utilisateur pour ',/,&
'                          le calcul des proprietes physiques.',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!     Message au premier passage pour indiquer que l'utilisateur a
!       rapatrie le sous-programme.
if (ipass.eq.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 1 - EFFET JOULE
!===============================================================================

if (ippmod(ieljou).ge.1) then


!     Attention, dans les modules electriques, la chaleur massique, la
!       conductivite thermique et la conductivite electriques sont
!       toujours dans le tableau PROPCE
!       qu'elles soient physiquement variables ou non.

!       On n'utilisera donc PAS les variables
!          =====================
!                                CP0, VISLS0(ISCALT)
!                                VISLS0(IPOTR) et VISLS0(IPOTI)

!       Informatiquement, ceci se traduit par le fait que
!                                ICP>0, IVISLS(ISCALT)>0,
!                                IVISLS(IPOTR)>0 et IVISLS(IPOTI)>0





!       Calcul de la temperature a partir de l'enthalpie
!       ------------------------------------------------

!       Ceci depend largement des choix utilisateur en
!         matiere de loi H-T (T en Kelvin)

!       On demande de fournir cette loi dans le sous programme usthht
!          (USERS/base/usthht.F)
!           usthht fournit en particulier un exemple d'interpolation
!            a partir d'une tabulation utilisateur
!           usthht en mode T->H sera utilise pour l'initialisation
!            de l'enthalpie dans useliv.

!       MODE = 1 : H=RTP(IEL,ISCA(IHM)) -> T=PROPCE(IEL,IPPROC(ITEMP))
  mode = 1

  do iel = 1, ncel
    call usthht (mode,                                            &
         rtp(iel,isca(ihm)),propce(iel,ipproc(itemp)))
  enddo


!       Masse volumique au centre des cellules
!       --------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la masse volumique ici
!       en renseignant PROPCE(IEL,IPCROM)
!       (meme si elle est uniforme ou constante).


!     Masse Vol : RO = ROM0 / (1+DILAR*(T-T0)
!         (Choudhary) semblable a Plard (HE-25/94/017)

!          avec sous-relaxation (sauf au premier pas de temps)

  temp0  = 300.d0
  rom0   = 2500.d0
  dilar  = 7.5d-5
  if (ntcabs.gt.1) then
    srrom1 = srrom
  else
    srrom1 = 0.d0
  endif

  ipcrom = ipproc(irom)
  do iel = 1, ncel
    rhonp1 = rom0 /                                               &
            (1.d0+ dilar * (propce(iel,ipproc(itemp))-temp0) )
    propce(iel,ipcrom) =                                          &
         srrom1*propce(iel,ipcrom)+(1.d0-srrom1)*rhonp1
  enddo


!       Viscosite moleculaire dynamique en kg/(m s)
!        ------------------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la viscosite ici
!       en renseignant PROPCE(IEL,IPCVIS)
!       (meme si elle est uniforme ou constante).


!     Viscosite : MU = EXP((AA/T-BB)-CC)
!          (Choudhary)
!      Plard (HE-25/94/017) ; limite a 1173K par C Delalondre

  ipcvis = ipproc(iviscl)
  aa     = 10425.d0
  bb     =   500.d0
  cc     =-6.0917d0

  do iel = 1, ncel
    if ( propce(iel,ipproc(itemp)) .gt. 1173.d0 ) then
      tp = propce(iel,ipproc(itemp))
    else
      tp= 1173.d0
    endif
    propce(iel,ipcvis) = exp( (aa/(tp-bb))+cc )
  enddo


!       Chaleur specifique J/(kg degres)
!       --------------------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la chaleur massique ici
!       en renseignant PROPCE(IEL,IPCPP)
!       (meme si elle est uniforme ou constante).


!        CP = 1381 (Choudhary)
!          coherent avec Plard (HE-25/94/017)

  ipccp  = ipproc(icp)
  do iel = 1, ncel
    propce(iel,ipccp) = 1381.d0
  enddo


!       Lambda/Cp en kg/(m s)
!       ---------------------

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCVSL)
!       (meme si elle est uniforme ou constante).


!         Lambda
!          On suppose Cp renseigne au prealable.

!          Plard (HE-25/94/017)

  ipcvsl = ipproc(ivisls(iscalt))

  do iel = 1, ncel
    xbr = 85.25d0                                                 &
         -5.93d-2*(propce(iel,ipproc(itemp))-tkelvi)              &
         +2.39d-5*(propce(iel,ipproc(itemp))-tkelvi)**2
    xkr = 16.d0*stephn*(1.4d0)**2*(propce(iel,ipproc(itemp)))**3  &
         /(3.d0*xbr)

    propce(iel,ipcvsl) = 1.73d0 + xkr
  enddo

! --- On utilise CP calcule  dans PROPCE ci dessus
  do iel = 1, ncel
    propce(iel,ipcvsl) = propce(iel,ipcvsl)/propce(iel,ipccp)
  enddo


!       Conductivite electrique en S/m
!       ==============================

!     ATTENTION :
!     =========
!       Dans le module electrique effet Joule, on fournira
!       OBLIGATOIREMENT la loi de variation de la conductivite ici
!       en renseignant PROPCE(IEL,IPCSIG)
!       (meme si elle est uniforme ou constante).


!         SIGMA  (Plard HE-25/94/017)

  ipcsig = ipproc(ivisls(ipotr))
  do iel = 1, ncel
    propce(iel,ipcsig) =                                          &
         exp(7.605d0-7200.d0/propce(iel,ipproc(itemp)))
  enddo

!     La conductivite electrique pour le potentiel imaginaire est
!       toujours implicitement prise egale a la conductivite
!       utilisee pour le potentiel reel.
!       IL NE FAUT PAS la renseigner.

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       Ce choix est fait en dur dans varpos.
!       Les pointeurs pour les deux existent quand meme.
!     Sinon, on pourrait faire ceci :
  if (1.eq.0) then
    if ( ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
      do iel = 1, ncel
        propce(iel,ipproc(ivisls(ipoti))) =                       &
             propce(iel,ipproc(ivisls(ipotr)))
      enddo
    endif
  endif
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!       Diffusivite variable a l'exclusion de l'enthalpie et du potentiel
!       -----------------------------------------------------------------
!     Pour le moment, il n'y a pas d'autres scalaires et
!                                                  on ne fait donc rien

endif

!===============================================================================
! 2 - ARC ELECTRIQUE
!===============================================================================

!     Les proprietes physiques sont a priori fournies par fichier
!       de donnees. IL n'y a donc rien a faire ici.

!      IF ( IPPMOD(IELARC).GE.1 ) THEN
!      ENDIF


!===============================================================================
! 3 - CONDUCTION IONIQUE
!===============================================================================

!     CETTE OPTION N'EST PAS ACTIVABLE

!--------
! Formats
!--------

 1000 format(/,                                                   &
' Module electrique: intervention utilisateur pour        ',/,    &
'                      le calcul des proprietes physiques.',/)

!----
! End
!----

return
end subroutine uselph


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
!>                               (cf. \ref ustsma)
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in,out] propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     ckupdc        work array for head loss terms
!> \param[in]     smacel        values of variables related to mass source
!>                              term. If ivar=ipr, smacel=mass flux
!_______________________________________________________________________________

subroutine usvist &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          iel, iccocg, inc
integer          ipcrom, ipcvst
double precision dudx, dudy, dudz, sqdu, visct, rom

double precision, dimension(:,:), pointer :: coefav
double precision, dimension(:,:,:), pointer :: coefbv
double precision, allocatable, dimension(:,:,:) :: gradv

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
allocate(gradv(ncelet,3,3))

! --- Physical quantity numbers in PROPCE (physical quantities defined
!     at each cell center)
ipcvst = ipproc(ivisct)
ipcrom = ipproc(irom  )

!===============================================================================
! 1.3 Compute velocity gradient
!===============================================================================

iccocg = 1
inc = 1

! Note: this example should produce an error if used with ivelco = 0,
!       so it should be updated

! Boundary condition pointers for gradients and advection
call field_get_coefa_v(ivarfl(iu), coefav)
call field_get_coefb_v(ivarfl(iu), coefbv)

call grdvec &
!==========
 ( iu  , imrgra , inc    , iccocg ,                      &
   nswrgr(iu) , imligr(iu) ,                             &
   iwarni(iu) , nfecra ,                                 &
   epsrgr(iu) , climgr(iu) , extrag(iu) ,                &
   rtpa(1,iu) , coefav , coefbv ,                        &
   gradv  )

!===============================================================================
! 1.4 Computation of the dynamic viscosity
!===============================================================================

do iel = 1, ncel

  ! --- Current dynamic viscosity and fluid density
  visct = propce(iel,ipcvst)
  rom   = propce(iel,ipcrom)

  ! --- Various computations
  dudx = gradv(iel,1,1)
  dudy = gradv(iel,2,1)
  dudz = gradv(iel,3,1)
  sqdu = sqrt(dudx**2+dudy**2+dudz**2)

  ! --- Computation of the new dynamic viscosity
  visct = max (visct,rom*sqdu)

  ! --- Store the new computed dynamic viscosity
  propce(iel,ipcvst) = visct

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
   dt     , rtp    , rtpa   , propce , propfb ,                   &
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
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
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

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)
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
!> \section cell_id Cells identification
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
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and preceding time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[out]    viscmx        mesh viscosity in X direction
!> \param[out]    viscmy        mesh viscosity in Y direction
!> \param[out]    viscmz        mesh viscosity in Z direction
!_______________________________________________________________________________

subroutine usvima &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
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

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(ndimfb,*)
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
