!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine uslag1
!================

!===============================================================================
! Purpose:
! -------

!    User subroutine of the Lagrangian particle-tracking module:

!    User subroutine for input of calculation parameters (Fortran commons).
!    This parameters concern physical, numerical and post-processing options.

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
use entsor
use lagdim
use lagpar
use lagran
use ihmpre

!===============================================================================

implicit none

! Local variables

integer          ii , ipv , icha
double precision sio2 , al2o3 , fe2o3 , cao

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
  iilagr = 0
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Particle-tracking mode
!===============================================================================

!     IILAGR = 0 : no particle tracking (default)

!            = 1 : particle-tracking one-way coupling

!            = 2 : particle-tracking two-way coupling

!            = 3 : particle tracking on frozen field
!                  (this option requires a calculation restart ISUITE=1,
!                   all Eulerian fields are frozen (pressure, velocities,
!                   scalars). This option is stronger than ICCVFG)

iilagr = 1

!===============================================================================
! 2. Particle-tracking calculation restart
!===============================================================================

!     ISUILA = 0 : no restart (default)
!            = 1 : restart (this value requires a restart on the continuous
!                  phase too, i.e. ISUITE = 1)

isuila = 0

!     Restart on volume and boundary statistics, and two-way coupling terms
!     (useful if ISUILA = 1)
!     (defaul off: 0 ; on: 1)

if (isuila.eq.1) isuist = 0

!===============================================================================
! 3. Particle tracking: specific models
!===============================================================================

!     IPHYLA = 0 : only transport modeling (default)
!            = 1 : equation on temperature (in Celsius degrees), diameter or mass
!            = 2 : pulverized coal combustion (only available if the continuous
!                   phase is a flame of pulverized coal)

iphyla = 0


! 3.1  equation on temperature, diameter or mass

if (iphyla.eq.1) then

!      equation on diameter
!      (default off: 0 ; on: 1)

  idpvar = 0

!      equation on temperature (in Celsius degrees)
!      (default off: 0 ; on: 1)
!      This option requires a thermal scalar for the continuous phase.

  itpvar = 0

!      equation on mass
!      (default off: 0 ; on: 1)

  impvar = 0

endif

! 3.2 coal fouling

!     Reference internal reports EDF/R&D: HI-81/00/030/A and HI-81/01/033/A

!     Evaluation of the probability for a particle to stick to a wall.
!     This probability is the ratio of a critical viscosity on the
!     viscosity of coal ashes

!              VISREF
!     P(Tp) = --------   for VISCEN >= VISREF
!              VISCEN

!           = 1 otherwise

!     The expression of J.D. Watt and T.Fereday (J.Inst.Fuel-Vol42-p99)
!     is used to evaluate the viscosity of the ashes

!                         ENC1 * 1.0D+7
!     Log  (10*VISCEN) = --------------- + ENC2
!        10                            2
!                        (Tp(C) - 150)

!     In literature, the range of the critical viscosity VISREF is between
!     8 Pa.s and 1.D7 Pa.s  For general purpose 1.0D+4 Pa.s is chosen

if (iphyla.eq.2) then

!       IENCRA = 0 no fouling (default)
!              = 1 fouling

!       * In uslag2.f90, the boundary on which the fouling can occur must be given
!       * Post-processing:  IENSI3 = 1 and IENCBD = 1 (10.2)

  iencra = 0

!     Example of definition of fouling criteria for each coal

! first (and single) coal ICHA = 1

  icha = 1

!       TPRENC : threshold temperature below which no fouling occurs
!                (in Celsius degrees)

  tprenc(icha) = 600.d0

!       VISREF : critical viscosity (Pa.s)

  visref(icha) = 10000.d0

!    > coal composition in mineral matters:
!       (with  SiO2 + Al2O3 + Fe2O3 + CaO + MgO = 100% in mass)

  sio2   =  36.0d0
  al2o3  =  20.8d0
  fe2o3  =   4.9d0
  cao    =  13.3d0

!       ENC1 and ENC2 : coefficients in Watt and Fereday expression

  enc1(icha) = 0.00835d0 * sio2 + 0.00601d0 * al2o3 - 0.109d0

  enc2(icha) = 0.0415d0 * sio2  + 0.0192d0 * al2o3                &
       + 0.0276d0 * fe2o3 + 0.016 * cao - 3.92d0

endif

!===============================================================================
! 4. Number of particles allowed simultaneously inside the computational domain
!===============================================================================

!     * Warning, memory is allocated with NBPMAX

nbpmax = 1000000

!===============================================================================
! 5. Calculation features for the dispersed phases
!===============================================================================


! 5.1 Additional variables
! ------------------------

!     * these additional variables are stored in ETTP and ETTPA arrays
!     * NVLS is the number of additional variables
!     * the upper limit is NUSVAR = 10 (fixed in block common lagpar.h)
!     * one access to additional variables in ETTP ETTPA using the pointer JVLS:

!               current step  -> ETTP(NBPT,JVLS(NVUS))
!               previous step -> ETTPA(NBPT,JVLS(NVUS))

!           NBPT is the number of the considered particle
!             (integer between 1 and NBPART),
!           NVUS is the number of the additional variable
!             (integer between 1 and NVLS),
!     * the integration of the associated differential stochastic equation
!       requires a user intervention in uslaed.f90 subroutine

nvls = 0


! 5.2 Stationary or unsteady continuous phase

!     * if stationary: ISTTIO = 1
!     * if unsteady: ISTTIO = 0
!     * if IILAGR = 3 then ISTTIO = 1

!     Remark: if ISTTIO = 0, then the statistical averages are RESET
!             at each lagrangian iteration

if (iilagr.ne.3) isttio = 0

! 5.3 Two-way coupling: (IILAGR = 2)

if (iilagr.eq.2) then

!     * number of absolute lagrangian iteration (i.e. with restart)
!       from which a time average for two-way coupling source terms is
!       computed (stationary source terms)
!     * if the Lagrangian iteration is lower than NSTITS, source terms are
!       unstationary: they are reset at each lagrangian iteration
!     * useful only if ISTTIO = 1.
!     * the min value for NSTITS is 1

  nstits = 1

!     two-way coupling for dynamic (velocities and turbulent scalars)
!     (default off: 0 ; on: 1)
!     (useful if ICCVFG = 0)

  ltsdyn = 0

!     two-way coupling for mass (if IPHYLA = 1 and IMPVAR = 1)
!     (default off: 0 ; on: 1)

  if(iphyla.eq.1 .and. (impvar.eq.1 .or. idpvar.eq.1)) ltsmas = 0

!     two-way coupling for thermal scalar (if IPHYLA = 1 and IMPVAR = 1, or IPHYLA = 2)
!     or for coal variables (if IPHYLA = 2)
!     (default off: 0 ; on: 1)

  if((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) ltsthe = 0

endif


! 5.4 Volume statistics
! ---------------------

!   5.4.1 Generic parameters
!   ~~~~~~~~~~~~~~~~~~~~~~~~~

!     Calculation of the volume statistics
!     (default off: 0 ; on: 1)

istala = 0


if (istala.eq.1) then

! Threshold for the management of volume statistics
! -------------------------------------------------
!  * the value of the seuil variable is a statistical weight.
!  * each cell of the mesh contains a statistical weight
!    (sum of the statistical weights of all the particles
!     located in the cell); seuil is the minimal value from
!    which the contribution in statistical weight of a particle
!    is not taken into account anymore in the full model
!    of turbulent dispersion, in the resolution of the
!    Poisson equation of correction of the mean velocities, and
!    in the writing of the listing and post-processing.
!

  seuil = 0.d0

! Calculation of the volume statistics from the absolute number
! of Lagrangian iterations
! * idstnt is a  absolute number of Lagrangian iterations
!   (i.e. including calculation restarts)


  idstnt = 1

! Stationary calculation from the absolute Lagrangian iteration nstist
! *  nstist is a  absolute number of Lagrangian iterations
!   (i.e. including calculation restarts) from which the statistics
!    are averaged in time.
! *  useful if the calculation is stationary (isttio=1)
! *  if the number of Lagrangian iterations is lower than nstits, the transmitted
!    source terms are unsteady (i.e. they are reset to zero ar each Lagrangian iteration)
! *  the minimal value acceptable for nstist is 1.

  nstist = idstnt


!   5.4.2 Volume statistical variables
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Activation of the calculation of the particle volume fraction
! Name of the mean : Part_vol_frac

  iactfv = 1

! Activation of the calculation of the particle velocity x-component (mean and variance)
! Name of the mean: Part_velocity_X
! Name of the variance: var_Part_velocity_X

  iactvx = 1

! Activation of the calculation of the particle velocity y-component (average and variance)
! Name of the mean: Part_velocity_Y
! Name of the variance: var_Part_velocity_Y

  iactvy = 1

! Activation of the calculation of the particle velocity z-component (average and variance)
! Name of the mean: Part_velocity_Z
! Name of the variance: var_Part_velocity_Z

  iactvz = 1

! Activation of the calculation of the particle residence time (mean and variance)
! Name of the mean: Part_resid_time
! Name of the variance: var_Part_resid_time

  iactts = 1


!   2) Specific models (iphyla = 1) following the chosen options:
!         Mean and variance of the temperature
!         Mean and variance of the diameter
!         Mean and variance of the mass

  if (iphyla.eq.1) then

    if (itpvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoTempPt'
      NOMLAV(IPV) = 'VaTempPt'
      ihslag(ipv)  = 2
    endif
    if (idpvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoDiamPt'
      NOMLAV(IPV) = 'VaDiamPt'
      ihslag(ipv)  = 2
    endif
    if (impvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoMassPt'
      NOMLAV(IPV) = 'VaMassPt'
      ihslag(ipv)  = 2
    endif

  else if (iphyla.eq.2) then

!   3) Pulverized coal (iphyla = 2) :
!         Mean and variance of the temperature
!         Mean and variance of the mass of reactive coal
!         Mean and variance of the mass of coke
!         Mean and variance of the diameter of the shrinking core

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoTempPt'
    NOMLAV(IPV) = 'VaTempPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoMchPt'
    NOMLAV(IPV) = 'VaMchPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoMcokPt'
    NOMLAV(IPV) = 'VaMcokPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoDckPt'
    NOMLAV(IPV) = 'VaDckPt'
    ihslag(ipv)  = 2

  endif

!   4) Additional volume statistical variables
!      ---------------------------------------
!     * If the user wishes other statistic calculations
!       than the standard ones, he must 1) prescribe
!       their number nvlsts, 2) prescribe their names,
!       3) prescribe ihslag and 4) intervene in the
!       user subroutines uslast and uslaen to implement
!       his new statistics (see the given examples)
!     * Default maximal number of additional statistics: 20.
!       (Otherwise, modify the nussta parameter is the
!       include file lagpar.h)

  nvlsts = 0

  if (nvlsts.gt.0) then
    do ii = 1,nvlsts
      ilvu(ii) = ipv + ii
      WRITE(NOMLAG(ILVU(II)),'(A6,I4.4)') 'MoyLag',II
      WRITE(NOMLAV(ILVU(II)),'(A6,I4.4)') 'VarLag',II
      ihslag(ilvu(ii))  = 1
    enddo
    ipv = ipv + nvlsts
  endif


!   6) Statistics per group:
!     ----------------------
!     * if the user wishes to calculate statistics per group of particles
!       (by default there is no statistics of this kind), he must:
!      1) prescribe nbclst the number of groups (limited to 100)
!      2) assign in uslag2 the group to which belongs each particle
!         through the iuslag array.
!
!     * Be careful, nbclst cannot be modified during a calculation restart
!       (isuila=1); even if the calculation of the statistics is not triggered yet
!       (istala=0).

  nbclst = 0

endif


!===============================================================================
! 6. Option concerning particle inlet
!===============================================================================

!    Continous particle injection during the time step
!    (and not only at the beginning the time step; this option
!     makes it possible to avoid bunches of particles in the vicinity
!     of the inlet zones)
!     (default off: 0 ; on: 1)

injcon = 0


!===============================================================================
! 7. Technic of variance reduction: cloning/fusion of the particles
!===============================================================================

!     Use of the Russian roulette
!                    default off : 0
!                            on  : 1 without Y+ calculation
!                                  2 with Y+ calculation
iroule = 0


!===============================================================================
! 8. Options concerning the numerical treatment of the dispersed phase
!===============================================================================

! Resolution of the Poisson equation for the particle mean velocity
! and correction of the particle instantaneous velocity
!      = 0: not correction of the velocities (default values)
!      = 1: correction of the instantaneous velocities

!     CAUTION: OPTION STRICTLY FOR DEVELOPERS; PLEASE LEAVE THE DEFAULT VALUE FOR A
!     ========= STANDARD USE OF THE CODE.           !

ilapoi = 0


!===============================================================================
! 9. Options concerning the treatment of the dispersed phase
!===============================================================================

!     CAUTION: In this version, the turbulent dispersion works only if
!     -------  the continuous phase is calculated with a k-eps or a Rij-eps model

!-->  Activation of the turbulent dispersion
!     (default on: 1 ; off: 0)

idistu = 1


!-->  Turbulent dispersion imposed to the fluid one.

!     If activated, then particle turbulent dispersion is
!     equal to the fluid-particle one. The crossing-trajectory effects
!     are suppressed ; it is then a case of turbulent diffusion. If the
!     simulated particle density is equal to the fluid density, then
!     we are simulating the displacement of fluid particles.
!     (default off: 0 ; on: 1)

idiffl = 0

!     modcpl :
!          = 0 for the incomplete model (default value)
!          > 0 for the full model, is equal the absolute number
!              of Lagrangian iterations from which the full model is activated
!              modcpl must not be lower than idstnt


modcpl = 0

!     idirla (=1 or 2 or 3) : 1st, 2nd or 3rd direction
!     of the full model. Corresponds to the main direction
!     of the flow. Allow to calculate a non-isotropic Lagrangian timescale
!     (default idirla=1)

if (modcpl.gt.0) idirla = 1


!===============================================================================
! 10. Options concerning the treatment of specific forces
!===============================================================================

!-->  Constants for the van der Waals forces
!     --------------------------------------
!    Hamaker constant for the particle/fluid/substrate system:

cstham = 6.d-20

!-->  Constants for the elecstrostatic forces
!    ----------------------------------------

!    Dielectric constant of the fluid (example: water at 293 K)

epseau = 80.10d0

!    Electrokinetic potential of the first solid (Volt)

phi1 = 50.d-3

!    Electrokinetic potential of the second solid (Volt)

phi2 = -50.d-3

!    Ionic force (mol/l)

fion = 1.d-2


!===============================================================================
! 11. Activation of Brownian motion
!===============================================================================


!--> Activation of Brownian motion:

!     (default off: 0 ; on: 1)

!     CAUTION: OPTION FOR DEVELOPERS ONLY
!     =========

lamvbr = 0

!===============================================================================
! 12. Activation of deposition model (Guingo & Minier, 2008)
!===============================================================================

!--> Activation of the deposition model

!    (default off: 0 ; on: 1)

!     CAUTION: OPTION FOR DEVELOPERS ONLY
!     =========

idepst = 0

!===============================================================================
! 13. Post-processing
!===============================================================================

! 13.1 Post-processing of the trajectories and particle displacements


!     CAUTION: COMPUTATIONALLY EXPENSIVE
!     ^^^^^^^^^^

!   13.1.1 Generic parameters
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~

!     Trajectory mode
!     (default off: 0 ; on: 1)

iensi1 = 0

!     Particle-displacement mode
!     (default off: 0 ; on: 1)

iensi2 = 0


!     Maximal number of particles to visualize = nliste
!     (default nbvis=nliste)
!     Be careful, nbvis must not be greater than nbpmax nor nliste

nbvis = nliste

!     Acquisition period of the data to visualize
!     (default nvisla = 1)

nvisla = 1

!     The liste array contains the numbers of particles that we want to
!     visualize (by default liste(...) = -1, no particle to visualize)
!     Be careful, if the number is negative, there is no visualization)
!
!     Example 1: the user wants to track the first nbvis

do ii = 1, nbvis
  liste(ii) = ii
enddo

!   > Example 2:
!     I want to track 3 5 67 23 1 76 35 36 ...etc..

!     liste(1) = 3
!     liste(2) = 5
!     liste(3) = 67
!     liste(4) = 23
!     liste(5) = 1
!     liste(6) = 76
!     liste(7) = 35
!     liste(8) = 36
!     ...etc...

!   > Remark: the holes, the repetitions of the liste array will be deleted,
!             and the numbers will be sorted in increasing order
!
!

!   13.1.2 Variables to visualize on the trajectories or the particles
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     Variable: velocity of the flow seen
!     (default off: 0 ; on: 1)

ivisv1  = 0

!     Variable: particle velocity
!     (default off: 0 ; on: 1)

ivisv2  = 0

!     Variable: residence time
!     (default off: 0 ; on: 1)

ivistp  = 0

!     Variable: diameter
!     (default off: 0 ; on: 1)

ivisdm  = 0

!     Variable: temperature
!     (default off: 0 ; on: 1)
if (iphyla.eq.1 .and. itpvar.eq.1) iviste  = 0

!     Variable: mass
!     (default off: 0 ; on: 1)

ivismp  = 0


if (iphyla.eq.2) then

!       Variable coal: temperature in Celsius degrees
!       (default off: 0 ; on: 1)

  ivishp = 0

!       Variable coal: diameter of the shrinking core
!       (default off: 0 ; on: 1)

  ivisdk  = 0

!       Variable coal: mass of reactive coal
!       (default off: 0 ; on: 1)

  ivisch  = 0

!       Variable coal: mass of coke
!       (default off: 0 ; on: 1)

  ivisck  = 0

endif


! 13.2 Boundary statistics: visualization of the particle/boundaries interactions
! ------------------------------------------------

!   13.2.1 Generic parameters
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~

!     Particle/boundary interaction mode
!     (default off: 0 ; on: 1)

iensi3 = 0

!    Stationary calculation of the boundary statistics from
!    the absolute Lagrangian iteration nstbor.
!     * nstbor is the absolute number of Lagrangian iterations
!       (i.e. including restarts) from which the statistics are averaged
!       (in time or by number of interactions)
!     * useful if the calculation is stationary (isttio=1)
!     * if the absolute number of Lagrangian iterations is inferior to
!       nstbor, the statistics are unsteady (i.e. they are reset to zero at each
!       Lagrangian iteration)

nstbor = 1

!   seuilf for the management of the boundary statistics
!    * the value of seuilf is a statistical weight
!    * Each boundary face has undergone a number of particle interactions
!      in term of statistical weight (sum of the statistical weights of all
!      the particles that interacted with the boundary face); seuilf is the
!      minimal value from which the contribution of a face (in statistical terms)
!      is not taken into account anymore in the writing of the listing and
!      post-processing.

seuilf = 0.d0


!   13.2.2 Information to be recorded
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!   * To activate them, the user has to set below
!     the corresponding keyword to 1.
!   * The default selection must be validated or modified by the user.
!   * By default the asked information for all the particle/wall interactions
!     are written in the same recording.
!   * The boundary statistic 'number of particle/boundary interactions' must be
!     selected to activate the particle average imoybr(...) = 2

!     Number of particle/boundary interactions
!     (default off: 0 ; on: 1)
inbrbd = 1

!     Particle mass flux associated to particle/boundary interactions
!     (default off: 0 ; on: 1)
iflmbd = 1

!     Angle between particle velocity and the plan of the boundary face
!     (default off: 0 ; on: 1)
iangbd = 0

!     Norm of particle velocity during the interation with the boundary face
!     (default off: 0 ; on: 1)
ivitbd = 0

!     Mass of fouled coal particles
!     (default off: 0 ; on: 1)
 if (iphyla.eq.2 .and. iencra.eq.1) iencbd = 0

!     Additional user information to be recorded
!     ------------------------------------------
!     (for instance, erosion rate, temperature..)
!    * these additional recordings are stored in the parbor array
!    * here we prescribe the nusbor number of additional recordings
!    * the max value of this number is nusbrd=10 (in lagpar.h)

nusbor = 0


!   13.2.3 Name of the recordings for display,
!          Average in time of particle average
!          of the boundary statistics
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!    * A priori the user intervenes only in the additional user information
!      to be recorded: he must prescribe the name of the recording as well as
!      the type of average that he wishes to apply to it for the writing
!      of the listing and the post-processing.


!    * The applied average is prescribed through the imoybr array:
!      - if imoybr(iusb(ii)) = 0 -> no average applied
!      - if imoybr(iusb(ii)) = 1 -> a time average is applied, i.e. the
!       statistic is divided by the last time step in the case of an unsteady
!       calculation with a number of iterations lower than nstbor; or that
!       the statistic is divided by the recording time in the case of a
!       steady calculation.
!      -if imoybr(iusb(ii)) = 2 -> a particle average is applied, i.e. the
!       statistic is divided by the number of recorded particle/boundary
!       interactions (in terms of statistical weight) in parbor(nfabor,inbr)
!       To use this average, inbrbd must be set to 1.
!    * The back-ups in the restart file are performed without applying
!      this average.
!    * The average is applied if the number of interactions (in statistical
!      weight) of the boundary face considered is greater than seuilf;
!      otherwise this average is set to zero.



!===============================================================================
! 14. Lagrangian listing
!===============================================================================
! Lagrangian period for the writing of the Lagrangian listing

ntlal = 1


!===============================================================================

return

end subroutine uslag1
