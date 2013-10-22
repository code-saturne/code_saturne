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

! Purpose:
! -------

! User subroutines for input of calculation parameters (Fortran modules).
!   These subroutines are called in all cases.

! If the Code_Saturne GUI is used, this file is not required (but may be
!   used to override parameters entered through the GUI, and to set
!   parameters not accessible through the GUI).

! Several routines are present in the file, each destined to defined
!   specific parameters.

! To modify the default value of parameters which do not appear in the
!   examples provided, code should be placed as follows:
!   - usipsu   for numerical and physical options
!   - usipes   for input-output related options

! As a convention, "specific physics" defers to the following modules only:
!   pulverized coal, gas combustion, electric arcs.

! In addition, specific routines are provided for the definition of some
!   "specific physics" options.
!   These routines are described at the end of this file and will be activated
!   when the corresponding option is selected in the usppmo routine.

!-------------------------------------------------------------------------------


!===============================================================================


subroutine usipes &
!================

 ( nmodpp )


!===============================================================================
! Purpose:
! --------

! User subroutine for the input of additional user parameters for
! input/output.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nmodpp           ! i  ! <-- ! number of active specific physics models       !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl

use coincl
use cs_coal_incl
use cs_fuel_incl
use cpincl
use elincl
use ppcpfu
use radiat

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, ipp, imom, idirac, icla, icha
integer idimve, iesp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipes'' must be completed',/,       &
'@       in file cs_user_parameters.f90',/,                       &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================


!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

!===============================================================================
! 1. Input-output (entsor)
!===============================================================================

! Frequency of log output

if (.false.) then
!< [init_01]
  ntlist = 1
!< [init_01]


endif

! Log (listing) verbosity

if (.false.) then

!< [init_02]
  do ii = 1, nvar
    iwarni(ii) = 1
  enddo

  iwarni(ipr) = 2
  iwarni(iu) = 2
  iwarni(iv) = 2
  iwarni(iw) = 2
!< [init_02]

endif

! --- probes output step

if (.false.) then

!< [init_03]
  nthist = 1
  frhist = -1.d0
!< [init_03]

endif

! --- Number of monitoring points (probes) and their positions
!     (limited to ncaptm=100)

if (.false.) then

!< [init_04]
  ncapt  = 4
  tplfmt = 1 ! time plot format (1: .dat, 2: .csv, 3: both)

  xyzcap(1,1) = 0.30d0
  xyzcap(2,1) = 0.15d0
  xyzcap(3,1) = 0.01d0

  xyzcap(1,2) = 0.30d0
  xyzcap(2,2) = 0.00d0
  xyzcap(3,2) = 0.01d0

  xyzcap(1,3) = 0.30d0
  xyzcap(2,3) =-0.08d0
  xyzcap(3,3) = 0.01d0

  xyzcap(1,4) = 0.60d0
  xyzcap(2,4) =-0.05d0
  xyzcap(3,4) = 0.01d0
!< [init_04]

endif

! --- current variable

!     As for other variables,
!       if we do not assign the following array values,
!       default values will be used

!     ichrvr( ) = chonological output (yes 1/no 0)
!     ilisvr( ) = logging in listing (yes 1/no 0)
!     ihisvr( ) = history output (number of probes and their numbers)
!     if ihisvr(.,1)  = -1, output for all probes

!     Note: Only the fist 8 characters of a name will be used in the most
!           detailed log.

if (.false.) then

!< [init_05]
  ! Current dynamic variables

  ! pressure variable
  ipp = ipprtp(ipr)
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! variable v1x
  ipp = ipprtp(iu)
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! v1y variable
  ipp = ipprtp(iv)
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! v1z variable
  ipp = ipprtp(iw)
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  if (itytur.eq.2) then

    ! turbulent kinetic energy
    ipp = ipprtp(ik)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! turbulent dissipation
    ipp = ipprtp(iep)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  elseif (itytur.eq.3) then

    ! Reynolds stresses
    ipp = ipprtp(ir11)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! Reynolds stresses
    ipp = ipprtp(ir22)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! Reynolds stresses
    ipp = ipprtp(ir33)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! Reynolds stresses
    ipp = ipprtp(ir12)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! Reynolds stresses
    ipp = ipprtp(ir13)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! Reynolds stresses
    ipp = ipprtp(ir23)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! turbulent dissipation
    ipp = ipprtp(iep)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  elseif (iturb.eq.50) then

    ! turbulent kinetic energy
    ipp = ipprtp(ik)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! turbulent dissipation
    ipp = ipprtp(iep)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! phi
    ipp = ipprtp(iphi)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! f_bar
    ipp = ipprtp(ifb)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  elseif (iturb.eq.51) then

    ! turbulent kinetic energy
    ipp = ipprtp(ik)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! turbulent dissipation
    ipp = ipprtp(iep)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! phi
    ipp = ipprtp(iphi)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! alpha
    ipp = ipprtp(ial)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  elseif (iturb.eq.60) then

    ! turbulent kinetic energy
    ipp = ipprtp(ik)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! omega
    ipp = ipprtp(iomg)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  elseif (iturb.eq.70) then

    ! Spalart-Allmaras variable (viscosity-like)
    ipp = ipprtp(inusa)
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

  endif

!< [init_05]
endif

! User scalar variables.

! We may modify here the arrays relative to user scalars, but scalars
!   reserved for specific physics are handled automatically. This explains
!   the tests on 'nscaus', which ensure that the targeted scalars are
!   truly user scalars.
! By specific physics, we mean only those which are handled in specific
!   modules of the code, such as coal, combustion, electric arcs (see usppmo).

if (.false.) then

!< [init_06]
  if (isca(1).gt.0.and.nscaus.ge.1) then
    ipp = ipprtp(isca(1))
    nomvar(ipp)  = 'Scalar 1'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  if (isca(2).gt.0.and.nscaus.ge.2) then
    ipp = ipprtp(isca(2))
    nomvar(ipp)  = 'Scalar 2'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

!< [init_06]
endif

! Other variables

if (.false.) then
!< [init_07]

  ! Density variable (output for post-processing only if variable or
  !                   in the case of specific physics)
  ipp = ipppro(ipproc(irom))
  ichrvr(ipp)   = max(irovar,nmodpp)
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! specific heat
  if (icp .gt. 0) then
    ipp = ipppro(ipproc(icp))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = 0
  endif

  ! laminar viscosity

  ipp = ipppro(ipproc(iviscl))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = 0

  ! turbulent viscosity
  ipp = ipppro(ipproc(ivisct))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Courant number
  ipp = ipppro(ipproc(icour))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

  ! Fourier number
  ipp = ipppro(ipproc(ifour))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

  ! 'csmago' variable for dynamic L.E.S. models
  !    (square of the Samgorinsky "constant")
  if (ismago.gt.0) then
    ipp = ipppro(ipproc(ismago))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! temporal means (example for moment 1)
  if (nbmomt.gt.0) then
    imom = 1
    ipp = ipppro(ipproc(icmome(imom)))
    nomvar(ipp) = 'Time Average 01'
    ichrvr(ipp) = 1
    ilisvr(ipp) = 1
    ihisvr(ipp,1) = -1
  endif

  ! total pressure (not defined in compressible case)
  if (ippmod(icompf).lt.0) then
    ipp = ipppro(ipproc(iprtot))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! local time step
  ipp = ippdt
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! characteristic time of transient velocity/pressure coupling
  ipp = ipptx
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ipp = ippty
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ipp = ipptz
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!< [init_07]
endif

! Specific physics variables

if (.false.) then

!< [init_08]
  ! Transported Variables
  !----------------------

  ! ---- Mass fraction of unburned (or fresh)  gas
  if (iygfm.gt.0) then
    ipp = ipprtp(isca(iygfm))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Mean mixture fraction
  if (ifm.gt.0) then
    ipp = ipprtp(isca(ifm))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Variance of mixture fraction
  if (ifp2m.gt.0) then
    ipp = ipprtp(isca(ifp2m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Fuel Mass fraction
  if (iyfm.gt.0) then
    ipp = ipprtp(isca(iyfm))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Variance of Fuel Mass fraction
  if (iyfp2m.gt.0) then
    ipp = ipprtp(isca(iyfp2m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  if (icoyfp.gt.0) then
    ipp = ipprtp(isca(icoyfp))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Enthalpy
  if (iscalt.gt.0) then
    ipp = ipprtp(isca(iscalt))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! ---- Soot
  if (isoot.eq.1) then
    ipp = ipprtp(isca(ifsm))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(inpm))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
!< [init_08]

  ! --> Variables for coal particles

!< [init_09]
  if (ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0               &
                          .or. ippmod(icolwc).ge.0) then

    do icla = 1, nclacp

      ! - Coal mass fraction (in class ICLA)
      ipp = ipprtp(isca(ixch(icla)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1

      ! - Number of particles for 1 kg mix (from class ICLA)
      ipp = ipprtp(isca(inp(icla)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1

      ! - Enthalpy J/kg (for class ICLA)
      ipp = ipprtp(isca(ih2(icla)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    enddo

  endif

  if (ippmod(iccoal).ge.0) then

      do icla = 1, nclacp

        ! Char mass fraction (in class ICLA)
        ipp = ipprtp(isca(ixck(icla)))
        ichrvr(ipp)  = 1
        ilisvr(ipp)  = 1
        ihisvr(ipp,1)= -1

        ! Coal mass fraction (in class ICLA)
        ipp = ipprtp(isca(ixch(icla)))
        ichrvr(ipp)  = 1
        ilisvr(ipp)  = 1
        ihisvr(ipp,1)= -1

        ! Number of particles for 1 kg mix (from class ICLA)
        ipp = ipprtp(isca(inp(icla)))
        ichrvr(ipp)  = 1
        ilisvr(ipp)  = 1
        ihisvr(ipp,1)= -1

        ! Enthalpy J/kg (for class ICLA)
        ipp = ipprtp(isca(ih2(icla)))
        ichrvr(ipp)  = 1
        ilisvr(ipp)  = 1
        ihisvr(ipp,1)= -1
      enddo

    endif

  ! Coal

  do icha = 1, ncharb

    ! - Mean of 1 mixture fraction
    !   (from light volatiles of char ICHA)
    ipp = ipprtp(isca(if1m(icha)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

    ! - Mean of 2 mixture fraction
    !   (from heavy volatiles of char ICHA)
    ipp = ipprtp(isca(if2m(icha)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

  enddo

  ! - Mean of 3 mixture fraction
  !  (C from heterogeneoux oxidation, of char, by O2)
  if (if3m.gt.0) then
    ipp = ipprtp(isca(if3m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! - Meam of (6 ?) mixture fraction
  !   (C from heterogeneous reaction between char and CO2)
  if (if3mc2.gt.0) then
    ipp = ipprtp(isca(if3mc2))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  if (if4m.gt.0) then
    ipp = ipprtp(isca(if4m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! - Mean of 5 mixture fraction
  !   (water vapor from drying)
  if (if5m.gt.0) then
    ipp = ipprtp(isca(if5m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Variance of 4 mixture fraction (air)
  if (ifvp2m.gt.0) then
    ipp = ipprtp(isca(ifvp2m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Mass fraction of CO2 or CO (relaxation to equilibrium)

  if (iyco2 .gt. 0) then
    ipp = ipprtp(isca(iyco2))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! HCN and NO

  if (ieqnox .ge. 1) then
    ipp = ipprtp(isca(iyhcn))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(iyno))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(ihox))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
!< [init_09]

  ! --> Variables for droplets

!< [init_10]
  do icla = 1, nclafu

    ! Fuel mass fraction
    ipp = ipprtp(isca(iyfol(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

    ! Number of droplets in mix (1/kg)
    ipp = ipprtp(isca(ing(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

    ! Fuel enthalpy (J/kg)
    ipp = ipprtp(isca(ih2(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1

  enddo
!< [init_10]

  ! --> Variables for carrying gas

!< [init_11]
  ! Mean of 1 mixture fraction (fuel vapor)
  if (ifvap.gt.0) then
    ipp = ipprtp(isca(ifvap))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Mean of 3 mixture fraction
  ! (carbon from heterogeneous oxidation of char)
  if (if7m.gt.0) then
    ipp = ipprtp(isca(if7m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Real component of the electrical potential
  if (ipotr.gt.0) then
    ipp = ipprtp(isca(ipotr))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Mass fraction of the different constituants of the phase
  if (ngazg .gt. 1) then
    do iesp = 1, ngazg-1
      ipp = ipprtp(isca(iycoel(iesp)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    enddo
  endif

  ! Specific variables for Joule effect for direct conduction
  ! Imaginary component of electrical potential
  if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
    ipp = ipprtp(isca(ipoti))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! --> Specific variables for electric arc in 3D
  !     vector potential components
  if (ippmod(ielarc) .ge. 2) then
    do idimve = 1, ndimve
      ipp = ipprtp(isca(ipotva(idimve)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    enddo
  endif
!< [init_11]

  ! Variables of State; User defined Variables
  !-------------------------------------------

!< [init_12]
  ! ---- Temperature
  ipp = ipppro(ipproc(itemp))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! ---- Fuel Mass fraction :    YM_Fuel
  ipp = ipppro(ipproc(iym(1)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! ---- Oxydizer Mass fraction : YM_Oxy
  ipp = ipppro(ipproc(iym(2)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! ---- Product Mass fraction : YM_Prod
  ipp = ipppro(ipproc(iym(3)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! --- Source term
  if (itsc.gt.0) then
    ipp = ipppro(ipproc(itsc))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! ---- Gas radiation

  if (iirayo.gt.0) then

    ! ---- Absorption Coefficient
    ipp = ipppro(ipproc(ickabs))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! ---- Term T^4
    ipp = ipppro(ipproc(it4m))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! ---- Term T^3
    ipp = ipppro(ipproc(it3m))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

  endif

  ! Premixed flame, LWC

  do idirac = 1, ndirac

    ipp = ipppro(ipproc(irhol(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'RHOL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iteml(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'TEML', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmel(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'FMEL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmal(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'FMAL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iampl(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'AMPL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(itscl(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'TSCL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(imaml(idirac)))
    write(nomvar(ipp),'(a4,i1)') 'MAML', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo

  !     - Mean Molar Mass
  ipp = ipppro(ipproc(immel))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1
!< [init_12]

  ! --> State variables for coal particles or fuel droplets

!< [init_13]
  do icla = 1, nclacp

    ! - Particles' Temperature K (of class ICLA)
    ipp = ipppro(ipproc(itemp2(icla)))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! - Particles' Density kg/m3 (of class ICLA)
    ipp = ipppro(ipproc(irom2(icla)))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! - Particles' Diameter m (of class ICLA)
    ipp = ipppro(ipproc(idiam2(icla)))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ! - Rate of coal consumption  (s-1) < 0
    !   (for class ICLA)
    ipp = ipppro(ipproc(igmdch(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Rate of light volatiles exhaust (s-1) < 0
    !   (for class ICLA)
    ipp = ipppro(ipproc(igmdv1(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Rate of heavy volatile exhaust (s-1) < 0
    !   (for class ICLA)
    ipp = ipppro(ipproc(igmdv2(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Rate of coal oxidation by O2 (s-1) < 0
    !   (from class ICLA)
    ipp = ipppro(ipproc(igmhet(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Rate of coal gazeification by CO2 (s-1) < 0
    !   (from class ICLA)
    if (ihtco2 .eq. 1) then
      ipp = ipppro(ipproc(ighco2(icla)))
      ichrvr(ipp)   = 0
      ilisvr(ipp)   = 0
      ihisvr(ipp,1) = -1
    endif

    ! - Rate of coal gazeification by H2O (s-1) < 0
    !   (from class ICLA)
    if (ihth2o .eq. 1) then
      ipp = ipppro(ipproc(ighh2o(icla)))
      ichrvr(ipp)   = 0
      ilisvr(ipp)   = 0
      ihisvr(ipp,1) = -1
    endif

    ! - Mass fraction (of class ICLA) in mix
    ipp = ipppro(ipproc(ix2(icla)))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    !- Heat flux (between gases and ICLA class droplets)
    ipp = ipppro(ipproc(ih1hlf(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Evaporation mass flow rate (s-1) < 0
    ipp = ipppro(ipproc(igmeva(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

    ! - Char combsution mass flow rate
    ipp = ipppro(ipproc(igmhtf(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1

  enddo
!< [init_13]

  ! --> State variables for carrier gas phase

!< [init_14]
  ! temperature of gas mixture
  ipp = ipppro(ipproc(itemp1))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of  CHx1m
  ipp = ipppro(ipproc(iym1(1)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of CHx2m
  ipp = ipppro(ipproc(iym1(2)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of CO
  ipp = ipppro(ipproc(iym1(3)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of H2S
  ipp = ipppro(ipproc(iym1(4)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of H2
  ipp = ipppro(ipproc(iym1(5)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of HCN
  ipp = ipppro(ipproc(iym1(6)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of NH3
  ipp = ipppro(ipproc(iym1(7)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of O2
  ipp = ipppro(ipproc(iym1(4)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of CO2
  ipp = ipppro(ipproc(iym1(5)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of H2O
  ipp = ipppro(ipproc(iym1(6)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of SO2
  ipp = ipppro(ipproc(iym1(11)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! mass fraction (among gases) of N2
  ipp = ipppro(ipproc(iym1(7)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Carbon balance
  if (ibcarbone.gt.0) then
    ipp = ipppro(ipproc(ibcarbone))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! Oxygen balance
  if (iboxygen.gt.0) then
    ipp = ipppro(ipproc(iboxygen))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! Hydrogen balance
  if (ibhydrogen.gt.0) then
    ipp = ipppro(ipproc(ibhydrogen))
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  endif

  ! Electric conductivity
  if (ipotr.gt.0) then
    ipp = ipppro(ipproc(ivisls(ipotr)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Joule effect Power
  if (iefjou.gt.0) then
    ipp = ipppro(ipproc(iefjou))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  ! Real component of the current density
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(idjr(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo

  ! Imaginary component of the current density
  if (ippmod(ieljou).eq.4) then
    do idimve = 1, ndimve
      ipp = ipppro(ipproc(idji(idimve)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    enddo
  endif

  if (ippmod(ielarc).ge.1) then

    ! Electromagnetic Forces (Laplace forces)
    do idimve = 1, ndimve
      ipp = ipppro(ipproc(ilapla(idimve)))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    enddo

    ! Absorption oefficient  or Radiative sources term
    if (ixkabe.gt.0) then
      ipp = ipppro(ipproc(idrad))
      ichrvr(ipp)  = 1
      ilisvr(ipp)  = 1
      ihisvr(ipp,1)= -1
    endif
  endif

  ! Electric charge (volumic)
  if (ippmod(ielion).ge.1) then
    ipp = ipppro(ipproc(iqelec))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

!< [init_14]
endif

!----
! Formats
!----


return
end subroutine usipes


!===============================================================================

subroutine user_field_parameters
!===============================

!===============================================================================
! Purpose:
! --------

! Define (redefine) key-value pairs on calculation fields.

! This subroutine is called at the end of the parameters initialization
! stage, after all other routines from this file have been called.

! Note that to determine which fields are defined in a computation, you
! may check the 'config.log' file after a first execution.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use ihmpre
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Local variables

logical       ilved, inoprv
integer       fldid, keyvis, idim1, iflpst, itycat, ityloc

!===============================================================================

! Example: force postprocessing of projection of some variables at boundary
!          with no reconstruction.
!          This is handled automatically if the second bit of a field's
!          'post_vis' key value is set to 1 (which amounts to adding 2
!          to that key value).
!
!          field_get_id returns -1 if field does not exist

!< [example_1]
call field_get_key_id('post_vis', keyvis)

fldid = ivarfl(iu)
call field_get_key_int(fldid, keyvis, iflpst)
if (iand(iflpst, 2) .eq. 0) then
  iflpst = ior(iflpst, 2)
  call field_set_key_int(fldid, keyvis, iflpst)
endif

fldid = ivarfl(ipr)
call field_get_key_int(fldid, keyvis, iflpst)
if (iand(iflpst, 2) .eq. 0) then
  iflpst = ior(iflpst, 2)
  call field_set_key_int(fldid, keyvis, iflpst)
endif
!< [example_1]

!-------------------------------------------------------------------------------

! Example: enforce existence of 'tplus' and 'tstar' fields, so that
!          a boundary temperature or Nusselt number may be computed using the
!          post_boundary_temperature or post_boundary_nusselt subroutines.
!          When postprocessing of these quantities is activated, those fields
!          are present, but if we need to compute them in the
!          cs_user_extra_operations user subroutine without postprocessing them,
!          forcing the definition of these fields to save the values computed
!          for the boundary layer is necessary.

!< [example_2]
itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces
ilved = .true. ! interleaved
inoprv = .false. ! no previous time step values needed

call field_get_id('tplus', fldid)
if (fldid.lt.0) then
  call field_create('tplus', itycat, ityloc, idim1, ilved, inoprv, fldid)
endif

call field_get_id('tstar', fldid)
if (fldid.lt.0) then
  call field_create('tstar', itycat, ityloc, idim1, ilved, inoprv, fldid)
endif
!< [example_2]

return

!===============================================================================

!----
! Formats
!----

return
end subroutine user_field_parameters

!===============================================================================

