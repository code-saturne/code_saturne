!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nphmx  ,                                                       &
   ibrom  ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Definition of physical variable laws.

! Warning:
! -------

! It is forbidden to modify turbulent viscosity "visct" here
!       =========
!    (a specifi subroutine is dedicated to that: usvist)


! icp = 1 must have been specified
!                ========================
!    in usini1 if we wish to define a varible specific heat
!    cp for phase iphas (otherwise: memory overwrite).


! ivisls = 1 must have been specified
!                   ========================
!    in usini1 if we wish to define a variable viscosity
!    viscls for phase iphas (otherwise: memory overwrite).


! Notes:
! -----

! This routine is called at the beginning of each time step

!    Thus, AT THE FIRST TIME STEP (non-restart case), the only
!    values initialized before this call are those defined
!      - in usini1 :
!             . density    (initialized at ro0)
!             . viscosity  (initialized at viscl0)
!      - in usiniv :
!             . calculation variables (initialized at 0 by defaut
!             or to the value given in the GUI or in usiniv)

! We may define here variation laws for cell properties, for:
!     - density                                    rom    kg/m3
!         (possibly also at boundary faces         romb   kg/m3)
!     - molecular viscosity                        viscl  kg/(m s)
!     - specific heat                              cp     J/(kg degrees)
!     - "diffusivities" associated with sclalars   viscls kg/(m s)


! The types of boundary faces at the previous time step are available
!   (except at the first time step, where arrays itypfb and itrifb have
!   not been initialized yet)


! It is recommended to keep only the minimum necessary in this file
!   (i.e. remove all unused example code)


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nphmx            ! i  ! <-- ! nphsmx                                         !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! ra ! --- ! work array                                     !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nphmx

integer          ibrom
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ivart, iclvar, iel, iphas
integer          ipcrom, ipbrom, ipcvis, ipccp
integer          ipcvsl, ith, iscal, ii
integer          iutile
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc
double precision xrtp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0. Initializations to keep
!===============================================================================

! --- Memory initialization

idebia = idbia0
idebra = idbra0


!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!  Each example is bounded by a test on "iutile", as a precaution.
!  Set iutile to 1 to activate the example.

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
!    Below, we define the same density law for all phases
!    Values of this property must be defined at cell centers
!      (and optionally, at boundary faces).
!  ===================================================================

!    The test on 'iutile' allows deactivating instructions (which are defined
!       only as a starting example)

iutile = 0
if(iutile.eq.1) then

  do iphas = 1, nphas ! Loop on phases

    ! Position of variables, coefficients
    ! -----------------------------------

    ! --- Number of the thermal variable for the current phase 'iphas'
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

    ! --- Rank of density for current phase 'iphas'
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

    do iel = 1, ncel
      xrtp = rtp(iel,ivart)
      propce(iel,ipcrom) = xrtp * (vara*xrtp+varb) + varc
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

    ! Note that when we prscribe the density at the boundary, it must be done
    ! at ALL boundary faces.
    !    ===

    ! ibrom = 1
    ! do ifac = 1, nfabor
    !   iel = ifabor(ifac)
    !   xrtp = coefa(ifac, iclvar)+rtp(iel, ivart)*coefb(ifac, iclvar)
    !   propfb(ifac, ipbrom) = xrtp * (vara*xrtp+varb) + varc
    ! enddo

    ! ifabor(ifac) is the cell adjacent to the boundary face

    ! Caution: ibrom = 1 is necessary for the law to be taken
    !                           into account.

  enddo ! --- Loop on phases
endif ! --- Test on 'iutile'


!===============================================================================
!  Example 2: variable viscosity as a function of temperature
!  =========
!    Below, we define the same viscosity law for all phases
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on 'iutile' allows deactivating instructions (which are defined
!       only as a starting example)

iutile = 0
if(iutile.eq.1) then

  do iphas = 1, nphas ! Loop on phases

    ! Position of variables, coefficients
    ! -----------------------------------

    ! --- Number of the thermal variable for the current phase 'iphas'
    !       To use user scalar 2 instead, write 'ivart = isca(2)'

    if (iscalt.gt.0) then
      ivart = isca(iscalt)
    else
      write(nfecra,9010) iscalt
      call csexit(1)
    endif

    ! --- Rank of molecular dynamic viscosity for current phase 'iphas'
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

  enddo ! --- Loop on phases
endif ! --- Test on 'iutile'


!===============================================================================
!  Example 3: specific heat as a function of temperature
!  =========
!    Below, we define the same viscosity law for all phases
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on 'iutile' allows deactivating instructions (which are defined
!       only as a starting example)

iutile = 0
if(iutile.eq.1) then

  do iphas = 1, nphas ! Loop on phases

    ! Position of variables, coefficients
    ! -----------------------------------

    ! --- Number of the thermal variable for the current phase 'iphas'
    !       To use user scalar 2 instead, write 'ivart = isca(2)'

    if (iscalt.gt.0) then
      ivart = isca(iscalt)
    else
      write(nfecra,9010) iscalt
      call csexit (1)
    endif

    ! --- Rank of the specific heat for current phase 'iphas'
    !     in 'propce', physical properties at element centers: 'ipccp'

    if(icp.gt.0) then
      ipccp  = ipproc(icp   )
    else
      ipccp  = 0
    endif

    ! --- Stop if Cp is not variable

    if(ipccp.le.0) then
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

  enddo ! --- Loop on phases
endif ! --- Test on 'iutile'


!===============================================================================
!  Example 4: Lambda/Cp a function of temperature for temperature or enthalpy
!  =========
!    Below, we define the same lambda/Cp ratio law for all phases
!    Values of this property must be defined at cell centers
!  ===================================================================

!    The test on 'iutile' allows deactivating instructions (which are defined
!       only as a starting example)

iutile = 0
if(iutile.eq.1) then

  do iphas = 1, nphas ! Loop on phases

    ! Position of variables, coefficients
    ! -----------------------------------

    ! --- Number of the thermal variable for the current phase 'iphas'
    !       To use user scalar 2 instead, write 'ivart = isca(2)'

    if (iscalt.gt.0) then
      ivart = isca(iscalt)
    else
      write(nfecra,9010) iscalt
      call csexit (1)
    endif

    ! --- Rank of Lambda/Cp of the thermal variable for current phase 'iphas'
    !     in 'propce', physical properties at element centers: 'ipcvsl'

    if(ivisls(iscalt).gt.0) then
      ipcvsl = ipproc(ivisls(iscalt))
    else
      ipcvsl = 0
    endif

    ! --- Stop if Lambda/CP is not variable

    if(ipcvsl.le.0) then
      write(nfecra,1010)                                          &
           iscalt, iscalt, ivisls(iscalt)
      call csexit (1)
    endif

    ! --- Rank of the specific heat for current phase 'iphas'
    !     in 'propce', physical properties at element centers: 'ipccp'

    if(icp.gt.0) then
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

    if(ipccp.le.0) then

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

  enddo ! --- Loop on phases
endif ! --- Test on 'iutile'


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

!    The test on 'iutile' allows deactivating instructions (which are defined
!       only as a starting example)

iutile = 0
if(iutile.eq.1) then

  do ii = 1, nscaus ! Loop on scalars

    ! --- Number of user scalar 'ii' in the lsit of scalars
    iscal = ii

    ! --- If it is a thermal variable, it has already been handled above
    ith = 0
    do iphas = 1, nphas
      if (iscal.eq.iscalt) ith = 1
    enddo

    ! --- If the variable is a fluctuation, its diffusivity is the same
    !       as that of the scalar to which it is attached:
    !       there is nothing to do here, we move on to the next variable
    !       without settign propce(iel,ipcvsl).

    ! We only handle here non-thermal variables which are not fluctuations
    if (ith.eq.0.and.iscavr(iscal).le.0) then

      ! Position of variables, coefficients
      ! -----------------------------------

      ! --- Number of the thermal variable for the current phase 'iphas'
      !       To use user scalar 2 instead, write 'ivart = isca(2)'

      iphas = 1
      if (iscalt.gt.0) then
        ivart = isca(iscalt)
      else
        write(nfecra,9010) iscalt
        call csexit (1)
      endif

      ! --- Rank of scalar's Lambda
      !     in 'propce', physical properties at element centers: 'ipcvsl'

      if(ivisls(iscal).gt.0) then
        ipcvsl = ipproc(ivisls(iscal))
      else
        ipcvsl = 0
      endif

      ! --- Stop if Lambda is not variable

      if(ipcvsl.le.0) then
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

  enddo ! --- Loop on phases
endif ! --- Test on 'iutile'


!===============================================================================

!===============================================================================
! Formats
!----

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@      usini1 indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ',I10                                   ,/,&
'@      usini1 indique que la diffusivite est uniforme        ',/,&
'@        IVISLS(',I10   ,') = ',I10   ,' alors que           ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                           &
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
'@      ISCALT = ',I10                                  ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usini1. Si un scalaire doit jouer le role de la       ',/,&
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
'@      usini1 specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usini1 or usphyv.',/,                                &
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
'@      usini1 specifies that the diffusivity is uniform',/,      &
'@        ivislc(',i10   ,') = ',i10   ,' while',/,               &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usini1 or usphyv.',/,                                &
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
'@      usini1. If a scalar should represent the,',/,             &
'@      temperature, check that iscalt has been defined',/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

#endif

!----
! End
!----

return
end subroutine
