!-------------------------------------------------------------------------------

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

!> \file lagran.f90
!> \brief Module for Lagrangian model.

module lagran

  !===========================================================================

  !         Trois modules complementaires
  !                            lagran qui porte les non dimensions
  !                            lagdim qui porte les dimensions variables
  !                            lagpar qui porte les parametres

  use lagpar

  implicit none

  !> \defgroup lagran Module for Lagrangian model

  !> \addtogroup lagran
  !> \{

  !=============================================================================

  !> \defgroup base Base

  !> \addtogroup base
  !> \{

  !> activates (>0) or deactivates (=0) the Lagrangian module
  !> the different values correspond to the following modellings:
  !> - = 1 Lagrangian two-phase flow in one-way coupling (no influence of
  !> the particles on the continuous phase)
  !> - = 2 Lagrangian two-phase flow with two-way coupling (influence of
  !> the particles on the dynamics of the continuous phase).
  !> Dynamics, temperature and mass may be coupled independently.
  !> - = 3 Lagrangian two-phase flow on frozen continuous phase. This option can
  !> only be used in case of a calculation restart. All the
  !> Eulerian fields are frozen (including the scalar fields). This option
  !> automatically implies \ref iccvfg = 1
  integer, save ::           iilagr

  !> activation (=1) or not (=0) of a Lagrangian calculation restart.
  !> The calculation restart file read when this option is activated (\ref ficaml)
  !> only contains the data related to the particles (see also \ref isuist)
  !> the global calculation must also be a restart calculation
  integer, save ::           isuila

  !> indicates the steady (=1) or unsteady (=0) state of the
  !> continuous phase flow
  !> in particular, \ref isttio = 1 is needed in order to:
  !> calculate stationary statistics in the volume or at the boundaries
  !> (starting respectively from the Lagrangian iterations \ref nstist and
  !> \ref nstbor)
  !> and calculate time-averaged two-way coupling source terms (from the
  !> Lagrangian iteration \ref nstits)
  !> useful if \ref iilagr=1 or \ref iilagr=2 (if \ref iilagr=3,
  !> then \ref isttio=1 automatically)
  integer, save ::           isttio

  !> \}

  !=============================================================================

  !> \defgroup particle_counter Particle counter
  !> (with and without statistical weight)

  !> \addtogroup particle_counter
  !> \{

  !> number of particles treated during one Lagrangian time step,
  !> must always be lower than \ref nbpmax.
  !> Always useful, but initialised and updated without intervention of the user
  integer, save ::           nbpart

  !> number of new entering particles
  integer, save ::           nbpnew

  !> number of in error deleted particles
  integer, save ::           nbperr

  !> total number of injected particles, since the beginning,
  !> including calculation restarts
  integer, save ::           nbptot

  !> contains normally exited particles
  !> and  in detection error exited particles
  integer, save ::           nbpout

  ! TODO
  integer, save ::           nbpert

  !> number of deposed particles for ever,
  !> but keeped for post-processing
  integer, save ::           ndepot

  !> number of deposed particles
  integer, save ::           nbpdep

  !> number of re-entrained particles
  integer, save ::           nbpres

  !> real value for \ref nbpart
  double precision, save ::  dnbpar

  !> real value for \ref nbpnew
  double precision, save ::  dnbpnw

  !> real value for \ref nbperr
  double precision, save ::  dnbper

  !> real value for \ref nbpout
  double precision, save ::  dnbpou

  !> real value for \ref nbpdep
  double precision, save ::  dnbdep

  !> real value for nbres
  double precision, save ::  dnbres

  !> number of cloned new particles
  integer, save ::           npclon

  !> number of killed by error particles
  integer, save ::           npkill

  !> number of cloned particles
  integer, save ::           npcsup

  !> real value for \ref npclon
  double precision, save ::  dnpclo
  !> real value for \ref npkill
  double precision, save ::  dnpkil
  !> real value for \ref npcsup
  double precision, save ::  dnpcsu

  !> \}

  !=============================================================================

  !> \defgroup specific_physic Specific physics

  !> \addtogroup specific_physic
  !> \{

  !> activates (>0) or deactivates (=0) the physical models associated to the
  !> particles:
  !> - = 1: allows to associate with the particles evolution
  !> equations on their temperature (in degrees Celsius), their diameter and
  !> their mass
  !> - = 2: the particles are pulverised coal particles.
  !> Evolution equations on temperature (in degree Celsius), mass of
  !> reactive coal, mass of char and diameter of the shrinking core are
  !> associated with the particles. This option is available only if the
  !> continuous phase represents a pulverised coal flame.
  integer, save ::           iphyla

  !> activation (=1) or not (=0) of an evolution equation on the particle
  !> temperature (in degrees Celsius)
  !> useful if \ref iphyla=1 and if there is a thermal scalar associated with
  !> the continuous phase
  integer, save ::           itpvar

  !> activation (=1) or not (=0) of an evolution equation on the particle
  !> diameter. useful if \ref iphyla = 1
  integer, save ::           idpvar

  !> activation (=1) or not (=0) of an evolution equation on the particle mass
  !> useful if si \ref iphyla = 1
  integer, save ::           impvar

  !> initialisation temperature (in degree Celsius) for the particles already
  !> present in the calculation domain when an evolution equation on
  !> the particle temperature is activated during a calculation (\ref iphyla =
  !> 1 and \ref itpvar = 1)
  !> useful if \ref isuila = 1 and \ref itpvar = 0 in the previous calculation
  double precision, save ::  tpart

  !> initialisation value for the specific heat (\f$ J.kg^{-1}.K^{-1} \f$)
  !> of the particles already present
  !> in the calculation domain when an evolution equation
  !> on the particle temperature is activated during a calculation
  !> (\ref iphyla = 1 and \ref itpvar = 1)
  !> useful if \ref isuila = 1 and \ref itpvar = 0 in the previous calculation
  double precision, save ::  cppart

  !==================================

  !> \defgroup deposition_model Particle deposition sub-model

  !> \addtogroup deposition_model
  !> \{

  !> - 0: no deposition submodel activated,
  !> - 1: deposition submodel used
  integer, save ::     idepst

  !> - 0: no DLVO conditions
  !> - 1: DLVO conditions
  integer, save ::     idlvo

  !> geometric parameters stored
  integer ngeol
  parameter (ngeol = 13)

  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jimark
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jdiel
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jdfac
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jdifel
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jtraj
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jptdet
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jinjst
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jryplu
  !> Additional pointer in the ITEPA array (contains the particule state)
  integer, save ::   jrinpf

  !> \}

  !=======================

  !> \defgroup reentrained_model Reentrained model

  !> \addtogroup reentrained_model
  !> \{

  !> - 0: no resuspension model
  !> - 1: resuspension model
  integer, save ::   ireent

  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jroll
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jnbasg
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jnbasp
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jdepo
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jfadh
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jmfadh
  !> Additional pointer in ITEPA and TEPA arrays (contains particule state)
  integer, save ::   jndisp

  !> Parameter of the particle resuspension model
  double precision, save :: espasg
  !> Parameter of the particle resuspension model
  double precision, save :: denasp
  !> Parameter of the particle resuspension model
  double precision, save :: modyeq
  !> Parameter of the particle resuspension model
  double precision, save :: rayasp
  !> Parameter of the particle resuspension model
  double precision, save :: rayasg

  !> \}

  !=======================

  !> \defgroup clogging_model Clogging model

  !> \addtogroup clogging_model
  !> \{

  !> - 0: no clogging model
  !> - 1: clogging model
  integer, save ::         iclogst

  !> Parameter of the particle clogging model
   double precision, save :: jamlim

  !> \}
  !> \}

  !=============================================================================

  !> \defgroup lagrangian_step Lagrangian time-step

  !> \addtogroup lagrangian_step
  !> \{

  !> absolute iteration number (including the restarts) in the Lagrangian
  !> module ( i.e. Lagrangian time step number)
  integer, save ::           iplas

  !> relative iteration number (including the restarts) in the Lagrangian
  !> module
  integer, save ::           iplar

  !> duration of a Lagrangian iteration
  double precision, save ::  dtp

  !> physical time of the Lagrangian simulation
  double precision, save :: ttclag

  !> \}

  !=============================================================================

  !> \defgroup error_indicator Error indicator

  !> \addtogroup error_indicator
  !> \{

  !> error indicator
  integer, save :: ierr

  !> \}

  !=============================================================================

  !> \defgroup particles_pointer Particles pointers
  !> for array \ref ettp

  !> \addtogroup particles_pointer
  !> \{

  !> pointer to particle X coordinate for array \ref ettp
  integer, save ::  jxp
  !> pointer to particle Y coordinate for array \ref ettp
  integer, save ::  jyp
  !> pointer to particle Z coordinate for array \ref ettp
  integer, save ::  jzp
  !> pointer to particle X velocity component for array \ref ettp
  integer, save ::  jup
  !> pointer to particle Y velocity component for array \ref ettp
  integer, save ::  jvp
  !> pointer to particle Z velocity component for array \ref ettp
  integer, save ::  jwp
  !> pointer to locally undisturbed X fluid velocity component
  !> for array \ref ettp
  integer, save ::  juf
  !> pointer to locally undisturbed Y fluid velocity component
  !> for array \ref ettp
  integer, save ::  jvf
  !> pointer to locally undisturbed Z fluid velocity component
  !> for array \ref ettp
  integer, save ::  jwf
  !> pointer to particle mass for array \ref ettp
  integer, save ::  jmp
  !> pointer to particle diameter for array \ref ettp
  integer, save ::  jdp
  !> pointer to particle and locally undisturbed fluid flow temperature
  !> (Celsius) for array \ref ettp
  integer, save ::  jtp
  !> pointer to particle and locally undisturbed fluid flow temperature
  !> (Celsius) for array \ref ettp
  integer, save ::  jtf
  !> pointer to particle specific heat for array \ref ettp
  integer, save ::  jcp
  !> pointer to coal particle temperature (\degresC) for array \ref ettp
  integer, save ::  jhp
  !> pointer to water mass (for coal) for array \ref ettp
  integer, save ::  jmwat
  !> pointer to mass of reactive coal of the coal particle for array \ref ettp
  integer, save ::  jmch
  !> pointer to mass of coke of the coal particle for array \ref ettp
  integer, save ::  jmck
  !> pointer to work array for the second order in time for array \ref ettp
  integer, save ::  jtaux
  !> pointer to additional user variable for array \ref ettp
  integer, save ::  jvls(nusvar)

  !> pointer to random number associated with a particle
  !> for array \ref tepa
  integer, save :: jrval
  !> pointer to particle residence time
  !> for array \ref tepa
  integer, save :: jrtsp
  !> pointer to particle statistic weight
  !> for array \ref tepa
  integer, save :: jrpoi
  !> pointer to particle emissivity
  !> for array \ref tepa
  integer, save :: jreps
  !> pointer to coal particle initial diameter
  !> for array \ref tepa
  integer, save :: jrd0p
  !> pointer to coal particle initial density
  !> for array \ref tepa
  integer, save :: jrr0p
  !> pointer to coal particle shrinking core diameter
  !> for array \ref tepa
  integer, save :: jrdck
  !> pointer to coal density
  !> for array \ref tepa
  integer, save :: jrhock

  !> pointer to number of the current cell containing the particle
  !> for \ref itepa array; this
  !> number is re-actualised during the trajectography step
  integer, save :: jisor
  !> pointer to number of the coal particle for \ref itepa array
  integer, save :: jinch
  !> pointer to class of the particle for \ref itepa array
  integer, save :: jclst
  !> pointer to number of additional variables related to the particles
  !> for \ref itepa array.
  !> The additional variables can be accessed in the arrays
  !> \ref ettp and \ret ettpa by means of the
  !> pointer \ref jvls:  \ref ettp(nbpt,jvls(ii)) and
  !> \ref ettpa(nbpt,jvls(ii)) (\ref nbpt is
  !> the index-number of the treated particle, and \ref ii an integer
  !> between 1 and \ref nvls)
  integer, save ::           nvls

  !> \}

  !=============================================================================

  !> \defgroup boundary_conditions Boundary conditions

  !> \addtogroup boundary_conditions
  !> \{

  !> number of boundary zones
  integer, save :: nfrlag
  !> continue injection or not
  integer, save :: injcon
  !> list of number of boundary zones
  integer, save :: ilflag(nflagm)


  !> for all the \ref nfrlag boundary zones previously identified, the number of classes
  !> \ref nbclas (a class is a set of particles sharing the same physical properties
  !> and the same characteristics concerning the injection in the calculation domain)
  !> of entering particles is given: \ref iusncl(izone) = \ref nbclas.
  !> By default, the number of particle classes is zero.
  !> The maximum number of classes is \ref nclagm (parameter stored in lagpar,
  !> whose default value is 20).
  integer, save :: iusncl(nflagm)

  !> for all the \ref nfrlag boundary zones
  !> previously identified, a particle boundary condition type is given.
  !> The categories of particle boundary condition types are marked out by the
  !> key words \ref ientrl, \ref isortl, \ref irebol, \ref idepo1, \ref idepo2, \ref iencrl.
  !> - if \ref iusclb(izone) = \ref ientrl, izone is a particle injection zone.
  !> For each particle class associated with this zone, information must be
  !> provided (see below). If a particle trajectory may cross an injection zone,
  !> then this particle leaves the calculation domain.

  !> - if \ref iusclb(izone) = \ref isortl, the particles interacting with the zone
  !>   \ref izone leave definitely the calculation domain.

  !> - if \ref iusclb(izone) = \ref irebol, the particles undergo an elastic
  !>   rebound on the boundary zone \ref izone.

  !> - if \ref iusclb(izone) = \ref idepo1, the particles settle definitely on the
  !>   boundary zone \ref izone. These particles leave the calculation domain
  !>   and are definitely erased from the calculation

  !> - if \ref iusclb(izone) = \ref idepo2, the particles settle definitevely
  !> on the boundary zone \ref izone and they are kept in the calculation
  !> - if \ref iusclb(izone) = \ref idepo2, the particles settle definitevely
  !> on the boundary zone \ref izone and they are kept in the calculation
  !> domain: the particles do not disappear after touching the boundary zone.
  !> However, using \ref idepo2 type zones necessitates more memory
  !> than using \ref idepo1 type zones.
  !> - if \ref iusclb(izone) = \ref iencrl, the particles which are coal particles
  !> (if \ref iphyla = 2) can become fouled up on the zone \ref izone. The
  !> slagging is a \ref idepo1 type deposit of the coal particle if a certain
  !> criterion is respected. Otherwise, the coal particle rebounds
  !> (\ref irebol type behaviour). This boundary condition type is available
  !> if \ref iencra = 1. A limit temperature \ref tprenc, a
  !> critical viscosity \ref visref and the coal composition
  !> in mineral matters must be given in the subroutine
  !> \ref uslag1.
  integer, save :: iusclb(nflagm)

  !> mean over a zone (if mean per zones is activated)
  integer, save :: iusmoy(nflagm)

  !> Some pieces of information must be given for each particle class
  !> associated with an injection zone.
  !> The first part consists in integers contained in the
  !> array \ref iuslag. There are at the most \ref ndlaim integers. These
  !> pieces of information must be provided for each class \ref iclas and each
  !> particle injection zone \ref izone.
  !> They are marked out by means of ``pointers'':
  !> - \ref iuslag(iclas,izone,ijnbp): number of particles to inject in
  !> the calculation domain per class and per zone.
  !> - \ref iuslag(iclas,izone,ijfre): injection period (expressed in number
  !> of time steps). If the period is null, then there is injection only
  !> at the first absolute Lagrangian time step (including the restart
  !> calculations).
  !> - \ref iuslag(iclas,izone,ijuvw): type of velocity condition:
  !> - if \ref iuslag(iclas,izone,ijuvw) = 1, the particle velocity vector is
  !> imposed, and its components must be given in the array \ref ruslag (see
  !> below).
  !> - if \ref iuslag(iclas,izone,ijuvw) = 0, the particle velocity is imposed
  !> perpendicular to the injection boundary face and with the norm
  !> \ref ruslag(iclas,izone,iuno).
  !> - if \ref iuslag(iclas,izone,ijuvw) = -1, the particle injection velocity
  !> is equal to the fluid velocity at the center of the cell
  !> neighbouring the injection boundary face.
  !> - \ref iuslag(iclas,izone,inuchl): when the particles are coal particles
  !> (\ref iphyla = 2), this part of the array contains the coal index-number,
  !> between 1 and \ref ncharb (defined by the user in the thermochemical
  !> file dp\_FCP, with  \ref ncharb$\leqslant$ncharm = 3).
  integer, allocatable, dimension(:,:,:) :: iuslag

  !> massic flow rate for a boudary zone
  double precision, save ::  deblag(nflagm)

  !> number of particles per class and per boudary zone
  integer, save ::  ijnbp
  !> injection frequency
  !> (if < 0 : particle are introduced only at first iteration
  integer, save ::  ijfre

  !> velocity condition type:
  !> - -1 imposed fluid velocity
  !> -  0 imposed fluid velocity along the normal of
  !> the bondary face, with \ref iuno norm.
  !> -  1 imposed velocity: \ref iupt \ref ivpt \ref iwpt must be given.
  !> -  2 velocity profile given by user.
  integer, save ::  ijuvw

  !> - 1 uniform distribution,
  !> - 2 presence rate profile given by user.
  integer, save ::  ijprpd

  !> - 1 constant temperature profile given in \ref uslag2
  !> - 2 temperature profile given by the user
  integer, save ::  ijprtp

  !> type of user profiles in \ref uslag2:
  !>  - 1: flat profile of diameter given in \ref uslag2
  !>  - 2: user profile to be given
  integer, save ::  ijprdp
  !> coal number of the particle (if \ref iphyla=2)
  integer, save ::  inuchl
  !> number of the statistics group
  integer, save ::  iclst

  !> some pieces of information must be given for each particle class
  !> associated with an injection zone.
  !> The second and last part consists in real numbers
  !> contained in the array \ref ruslag. There are at the most
  !> \ref ndlagm such real numbers. These pieces of information must
  !> be provided for each class \ref iclas and each particle injection zone
  !> \ref izone. They are marked out by means of ``pointers'':
  double precision, allocatable, dimension(:,:,:) :: ruslag

  !> particle velocity magnitude
  integer, save ::  iuno
  !> particle u component by class and zone
  integer, save ::  iupt
  !> particle v component by class and zone
  integer, save ::  ivpt
  !> particle w component by class and zone
  integer, save ::  iwpt
  !> particle temperature
  integer, save ::  itpt
  !> particle diameter
  integer, save ::  idpt
  !> particle diameter variance
  integer, save ::  ivdpt
  !> density
  integer, save ::  iropt
  !> particle specific heat
  integer, save ::  icpt
  !> particle weight
  integer, save ::  ipoit
  !> flow rate
  integer, save ::  idebt
  !> particle emissivity
  integer, save ::  iepsi
  !> particle temperature
  integer, save ::  ihpt
  !> water mass in coal particles
  integer, save ::  imwat
  !> active coal mass
  integer, save ::  imcht
  !> coal mass
  integer, save ::  imckt
  !> diameter of shrinking core
  integer, save ::  idckt

  !> \}

  !=============================================================================

  !> \defgroup statistics Statistics

  !> \addtogroup statistics
  !> \{

  !> mean dispersed phase velocity X component
  integer, save ::  ilvx
  !> mean dispersed phase velocity Y component
  integer, save ::  ilvy
  !> mean dispersed phase velocity Z component
  integer, save ::  ilvz
  !> sum of the statistical weights
  integer, save ::  ilpd
  !> dispersed phase volumetric concentration
  integer, save ::  ilfv
  !> recidence time
  integer, save ::  ilts
  !> phase temperature (\f$ \degresC \f$)
  integer, save ::  iltp
  !> dispersed phase mean diameter
  integer, save ::  ildp
  !> dispersed phase mean mass
  integer, save ::  ilmp
  !> temperature of the coal particle cloud (\f$ \degresC \f$)
  integer, save ::  ilhp
  !> water mass
  integer, save ::  ilmwat
  !> mass of reactive coal of the coal particle cloud
  integer, save ::  ilmch
  !> mass of coke of the coal particle cloud
  integer, save ::  ilmck
  !> shriking core diameter
  integer, save ::  ildck
  !> supplementary user volumetric statistics
  integer, save ::  ilvu(nussta)

  ! TODO : absence de \texttt{ilmdk} ? shrinking core diameter of the coal particle cloud

  ! TODO
  integer, save ::  iactfv
  ! TODO
  integer, save ::  iactvx
  ! TODO
  integer, save ::  iactvy
  ! TODO
  integer, save ::  iactvz
  ! TODO
  integer, save ::  iactts

  !> activation (=1) or not (=0) of the calculation of the volume
  !> statistics related to the dispersed phase.
  !> if \ref istala = 1, the calculation of the statistics is activated
  !> starting from the absolute iteration (including the restarts) \ref idstnt.
  !> by default, the statistics are not stationary (reset to zero at every
  !> Lagrangian iteration). But if \ref isttio=1, since the flow is steady,
  !> the statistics will be averaged overt he different time steps.
  !> the statistics represent the significant results on the particle cloud
  integer, save ::  istala

  !> during a Lagrangian calculation restart, indicates whether the particle
  !> statistics (volume and boundary) and two-way coupling terms are to be read
  !> from a restart file (=1) or reinitialised (=0).
  !> The file to be read is \ref ficmls. Useful if \ref isuila = 1
  integer, save ::  isuist

  !> number of additional user volume statistic
  !> the additional statistics (or their cumulated value in the stationary
  !> case) can be accessed in the array \ref statis by means of the pointer
  !> \ref ilvu: \ref statis (iel,ilvu(ii))
  !> (\ref iel is the cell index-number and \ref ii an integer between
  !> 1 and \ref nvlsts) useful if \ref istala = 1
  integer, save ::  nvlsts

  !> absolute Lagrangian iteration number (includings the restarts) after
  !> which the calculation of the volume statistics is activated
  !> useful if \ref istala = 1
  integer, save ::  idstnt

  !>
  !> absolute Lagrangian iteration number (includings the restarts) after
  !> which the volume statistics are cumulated over time (they are then said
  !> to be stationary).
  !> if the absolute Lagrangian iteration number is lower than \ref nstist,
  !> or if the flow is unsteady (\ref isttio=0), the statistics are reset
  !> to zero at every Lagrangian iteration (the volume statistics are then said
  !> to be non-stationary).
  !> useful if \ref istala=1 and \ref isttio=1
  integer, save ::  nstist

  !> number of iterations during which stationary volume statistics have
  !> been cumulated.
  !> useful if \ref istala=1, \ref isttio=1 and if \ref nstist is
  !> inferior or equal to the current Lagrangian iteration.
  !> \ref npst is initialised and updated automatically by the code, its
  !> value is not to be modified by the user
  integer, save ::  npst

  !> number of iterations during which volume statistics have been
  !> calculated (the potential iterations during which non-stationary
  !> statistics have been calculated are counted in \ref npstt).
  !> useful if \ref istala=1.
  !> \ref npstt is initialised and updated automatically by the code,
  !> its value is not to be modified by the user
  integer, save ::  npstt


  !> if the volume statistics are calculated in a stationary way, \ref tstat
  !> represents the physical time during which the statistics have been cumulated.
  !> if the volume statistics are calculated in a non-stationary way,
  !> then \ref tstat=dtp (it is the Lagrangian time step, because the
  !> statistics are reset to zero at every iteration).
  !> useful if \ref istala=1.
  !> \ref tstat is initialised and updated automatically by the code,
  !> its value is not to be modified by the user
  double precision, save ::  tstat

  !> every cell of the calculation domain contains a certain quantity of
  !> particles, representing a certain statistical weight (sum of the
  !> statistical weights of all the particles present in the cell). \ref seuil
  !> is the limit statistical weight value, below which the contribution of the
  !> cell in term of statistical weight is not taken into account in the volume
  !> statistics (for the complete turbulent dispersion model, in the
  !> Poisson's equation used to correct the mean velocities
  !> or in the listing and post-processing outputs). useful if \ref istala = 1
  double precision, save ::  seuil

  !> name of the volumetric statistics, displayed in the listing
  !> and the post-processing files.
  !> The default value is given above, with ``XXXX''
  !> representing a four digit number (for instance 0001, 0011 ...).
  !> useful if \ref istala = 1.
  !> Warning: this name is also used to reference information in the
  !> restart file (\ref isuist =1).
  !> If the name of a variable is changed between two calculations,
  !> it will not be possible to read its value from the restart file
  character*32, save ::      nomlag(nvplmx)

  ! TODO
  character*32, save ::      nomlav(nvplmx)

  !> historic statistics options
  integer, save ::           ihslag(nvplmx)

  !> statistic per zone and per class
  integer, save ::           nbclst

  !> \}

  !============================================================================

  !> \defgroup source_terms Source terms

  !> \addtogroup source_terms
  !> \{

  !> activation (=1) or not (=0) of the two-way coupling on the dynamics
  !> of the continuous phase.
  !> useful if \ref iilagr = 2 and \ref iccvfg = 0
  integer, save ::  ltsdyn

  !> activation (=1) or not (=0) of the two-way coupling on the mass.
  !> useful if \ref iilagr = 2, \ref iphyla = 1 and \ref impvar = 1
  integer, save ::  ltsmas

  !> if \ref iphyla = 1 and \ref itpvar = 1, \ref ltsthe
  !> activates (=1) or not (=0) the two-way coupling on temperature.
  !> if \ref iphyla = 2, \ref ltsthe activates (=1) or not (=0) the
  !> two-way coupling on the eulerian variables related to pulverised
  !> coal combustion.
  !> useful if \ref iilagr = 2
  integer, save ::  ltsthe

  !> explicit source term for the continuous phase X velocity
  integer, save ::  itsvx

  !> explicit source term for the continuous phase Y velocity
  integer, save ::  itsvy

  !> explicit source term for the continuous phase Z velocity
  integer, save ::  itsvz

  !> implicit source term for the continuous phase velocity and
  !> for the turbulent energy if the \f$k-\varepsilon\f$ model is used
  integer, save ::  itsli

  !> explicit source term for the turbulent dissipation and the
  !> turbulent energy if the \f$k-\varepsilon\f$ turbulence model is used
  !> for the continuous phase
  integer, save ::  itske

  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr11

  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr12
  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr13
  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr22
  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr23
  !> source term for the Reynolds stress
  !> and the turbulent dissipation if the \f$R_{ij}-\varepsilon\f$
  !> turbulence model is used for the continuous phase
  integer, save ::  itsr33

  !> explicit thermal source term for the thermal scalar of the continuous phase
  integer, save ::  itste

  !> implicit thermal source term for the thermal scalar of the continuous phase
  integer, save ::  itsti

  !> mass source term
  integer, save ::  itsmas

  !> source term for the light volatile matters
  integer, save ::  itsmv1(ncharm2)

  !> source term for the heavy volatile matters
  integer, save ::  itsmv2(ncharm2)

  !> source term for the carbon released during heterogeneous combustion
  integer, save ::  itsco

  !> Variance of the air scalar
  integer, save ::  itsfp4

  !> number of absolute Lagrangian iterations (including the restarts)
  !> after which a time-average of the two-way coupling source terms is
  !> calculated.
  !> indeed, if the flow is steady (\ref isttio=1), the average quantities
  !> that appear in the two-way coupling source terms can be calculated over
  !> different time steps, in order to get a better precision.
  !> if the number of absolute Lagrangian iterations is strictly inferior to
  !> \ref nstits, the code considers that the flow has not yet reached its
  !> steady state (transition period) and the averages appearing in the source
  !> terms are reinitialised at each time step, as it is the case for unsteady
  !> flows (\ref isttio=0).
  !> useful if \ref iilagr = 2 and \ref isttio = 1
  integer, save ::  nstits

  !> number of time steps for source terms accumulations
  integer, save ::  npts

  !> nomber of cells, whose vulumetric rate DODO (concentration ?)is greather than 0.8
  integer, save ::  ntxerr


  !> maximum volumetric concentration reached
  double precision, save ::  vmax

  !> maximum massic concentration reached
  double precision, save ::  tmamax

  !> \}

  !=============================================================================

  !> \defgroup fusion_cloning Particles cloning and fusion

  !> \addtogroup fusion_cloning
  !> \{

  ! TODO : INDICATEUR D ACTIVATION DE LA ROULETTE RUSSE
  integer, save ::           iroule

  !> \}

  !=============================================================================
  ! TODO
  !> \defgroup encrustation Encrustation

  !> \addtogroup encrustation
  !> \{

  !> activates (=1) or not (=0) the option of coal particle fouling.
  !> It then is necessary to specify the domain boundaries
  !> on which fouling may take place. useful if \ref iphyla = 2
  integer, save ::  iencra


  !> encrustation data
  integer, save ::  npencr

  !> encrustation data
  double precision, save ::  enc1(ncharm2)
  !> encrustation data
  double precision, save ::  enc2(ncharm2)

  !> limit temperature (in degree Celsius) below which the coal particles do
  !> not cause any fouling (if the fouling model is activated),
  !> useful if \ref iphyla = 2 and \ref iencra = 1
  double precision, save ::  tprenc(ncharm2)

  !>
  !> ash critical viscosity in \f$ kg.m^{-1}.s^{-1} \f$, in the fouling model
  !> cf J.D. Watt et T. Fereday (J.Inst.Fuel, Vol.42-p99).
  !> useful if \ref iphyla = 2 and \ref iencra = 1
  double precision, save ::  visref(ncharm2)

  !> encrustation data
  double precision, save ::  dnpenc

  !> \}

  !=============================================================================

  !> \defgroup physico_chemical Physico-chemical (DLVO) parameters

  !> \addtogroup physico_chemical
  !> \{

  !> Hamaker constant for the particle/fluid/substrate system
  double precision, save ::  cstham

  !> Dielectric constant of the fluid
  double precision, save ::  epseau

  !> Electrokinetic potential of the first solid
  double precision, save ::  phi1

  !> Electrokinetic potential of the second solid
  double precision, save ::  phi2

  !> Ionic force
  double precision, save ::  fion

  !> Faraday constant (C/mol)
  double precision cstfar
  parameter(cstfar = 9.648d4)

  !> Vacuum permittivity (F/m)
  double precision epsvid
  parameter(epsvid = 8.854d-12)

  !> Boltzmann constant (J/K)
  double precision kboltz
  parameter(kboltz = 1.38d-23)

  !> Cut-off distance for adhesion forces (assumed to be the Born distance) (m)
  double precision dcutof
  parameter(dcutof = 1.65d-10)

  !> Characteristic retardation wavelength (m) for Hamaker constant
  double precision lambwl
  parameter(lambwl = 1000.d-9)

  !> \}

  !=============================================================================

  !> \defgroup brownian Brownian motion

  !> \addtogroup brownian
  !> \{

  !> brownnian motion activation
  integer, save :: lamvbr

  !> \}

  !=============================================================================

  !> \defgroup time_scheme Time scheme, turbulent disperion and Poisson's equation

  !> \addtogroup time_scheme
  !> \{

  !> number of lagrangian under step (1 or 2)
  integer, save ::  nor

  !> order of integration for the stochastic differential equations
  !> - = 1 integration using a first-order scheme
  !> - = 2 integration using a second-order scheme
  integer, save ::  nordre

  !> activates (>0) or not (=0) the complete turbulent dispersion model.
  !> When \ref modcpl is strictly positive, its value is interpreted as the
  !> absolute Lagrangian time step number (including restarts) after which the
  !> complete model is applied.
  !> Since the complete model uses volume statistics, \ref modcpl must
  !> either be 0 or be larger than \ref idstnt.
  !> Useful if \ref istala = 1
  integer, save ::  modcpl

  !> direction (1=x, 2=y, 3=z) of the complete model.
  !> it corresponds to the main directions of the flow.
  !> useful if \ref modcpl > 0
  integer, save ::  idirla

  !> activation (=1) or not (=0) of the particle turbulent dispersion.
  !> The turbulent dispersion is compatible only with the RANS turbulent models
  !> (\f$k-\varepsilon\f$, $R_{ij}-\varepsilon$, v2f or $k-\omega\f$).
  !> (\ref iturb=20, 21, 30, 31, 50 or 60).
  integer, save ::  idistu

  !> \ref idiffl=1 suppresses the crossing trajectory effect, making
  !> turbulent dispersion for the particles identical to the turbulent
  !> diffusion of fluid particles.
  !> useful if \ref idistu=1
  integer, save ::  idiffl

  !> activation (=1) or not (=0) of the solution of a Poisson's equation for
  !> the correction of the particle instantaneous velocities
  !> (in order to obtain a null divergence).
  !> this option is not validated and reserved to the development team.
  !> Do not change the default value
  integer, save ::  ilapoi

  !> \}

  !=============================================================================

  !> \defgroup boundary_interactions Particles/boundary interactions statistics

  !> \addtogroup boundary_interactions
  !> \{

  !> number additional user data to record for the calculation
  !> of additional boundary statistics in \ref parbor.
  !> useful if \ref iensi3=1
  integer, save ::  nusbor

  !> number of absolute Lagrangian iterations (including the restarts)
  !> after which the statistics at the boundaries are considered stationary
  !> are veraged (over time or over the number of interactions).
  !> If the number of absolute Lagrangian iterations is lower than \ref nstbor,
  !> or if \ref isttio=0, the statistics are reset to zero at every
  !> Lagrangian iteration (non-stationary statistics).
  !> useful if \ref iensi3=1 and \ref isttio=1
  integer, save ::  nstbor

  !> number of iterations during which stationary boundary statistics have
  !> been cumulated.
  !> useful if \ref iensi3=1, \ref isttio=1 and \ref nstbor inferior
  !> or equal to the current Lagrangian iteration.
  !> \ref npstf is initialised and updated automatically by the code,
  !> its value is not to be modified by the user
  integer, save ::  npstf

  !> number of iterations during which boundary statistics have
  !> been calculated
  !> (the potential iterations during which non-stationary
  !> statistics have been calculated are counted in \ref npstft).
  !> useful if \ref iensi3=1.
  !> \ref npstft is initialised and updated automatically by the code,
  !> its value is not to be modified by the user
  integer, save ::  npstft

  !> activation (=1) or not (=0) of the recording of the number of particle/boundary
  !> interactions, and of the calculation of the associated boundary statistics.
  !> \ref inbrd = 1 is a compulsory condition to use the particulate average
  !> \ref imoybr = 2.
  !> useful if \ref iensi3=1
  integer, save ::  inbrbd

  !> activation (=1) or not (=0) of the recording of the particulate mass flow
  !> related to the particle/boundary interactions, and of the calculation of
  !> the associated boundary statistics.
  !> \ref inbrd = 1 is a compulsory condition to use \ref iflmbd=1.
  !> useful if \ref iensi3=1 and \ref inbrbd=1
  integer, save ::  iflmbd

  !> activation (=1) or not (=0) of the recording of the angle between a
  !> particle trajectory and a boundary face involved in a particle/boundary
  !> interaction, and of the calculation of the associated boundary statistics.
  !> useful if \ref iensi3=1
  integer, save ::  iangbd

  !> activation (=1) or not (=0) of the recording of the velocity of a particle
  !> involved in a particle/boundary interaction, and of the calculation of
  !> the associated boundary statistics.
  !> useful if \ref iensi3=1
  integer, save ::  ivitbd

  ! TODO
  integer, save ::  iencnbbd
  ! TODO
  integer, save ::  iencmabd
  ! TODO
  integer, save ::  iencdibd
  ! TODO
  integer, save ::  iencckbd

  !> number of particle/boundary interactions
  integer, save ::  inbr

  !> \ref iflm: particle mass flow at the boundary faces
  integer, save ::  iflm

  !> \ref iang: mean interaction angle with the boundary faces
  integer, save ::  iang

  !> \ref ivit: mean interaction velocity with the boundary faces
  integer, save ::  ivit

  ! TODO
  integer, save ::  iencnb
  ! TODO
  integer, save ::  iencma
  ! TODO
  integer, save ::  iencdi
  ! TODO
  integer, save ::  iencck

  !> supplementary user boundary statistics
  integer, save ::  iusb(nusbrd)

  !> the recordings in \ref parbor at every particle/boundary interaction are
  !> cumulated values (possibly reset to zero at every iteration in the
  !> non-stationary case). They must therefore be divided by a quantity to
  !> get boundary statistics. The user can choose between two average types:.
  !> - = 0: no average is applied to the recorded cumulated values.
  !> - = 1: a time-average is calculated. The cumulated value
  !> is divided by the physical duration in the case of stationary
  !> averages (\ref isttio=1). The cumulated value is divided by the value of
  !> the last time step in the case of non-stationary averages (\ref isttio=0),
  !> and also in the case of stationary averages while the
  !> absolute Lagrangian iteration number is inferior to \ref nstbor.
  !> - = 2: a particulate average is calculated. The cumulated
  !> value is divided by the number of particle/boundary interactions (in terms
  !> of statistical weight) recorded in \ref parbor(nfabor,inbr). This average
  !> can only be calculated when \ref inbrbd=1. The average is calculated if
  !> the number of interactions (in statistical weight) of the considered
  !> boundary face is strictly higher than \ref seuilf, otherwise the average
  !> at the face is set to zero.
  !> only the cumulated value is recorded in the restart file.
  !> useful if \ref iensi3=1
  integer, save ::  imoybr(nusbrd+10)

  ! TODO
  integer, save ::  inclg
  ! TODO
  integer, save ::  iscovc

  !> if the recording of the boundary statistics is stationary, \ref tstatp
  !> contains the cumulated physical duration of the recording of the boundary
  !> statistics.
  !> if the recording of the boundary statisticss is non-stationary, then
  !> \ref tstat=dtp (it is the Lagrangian time step, because the
  !> statistics are reset to zero at every time step).
  !> useful if \ref iensi3=1
  double precision, save ::  tstatp

  !> every boundary face of the mesh undergoes a certain number of
  !> interactions with particles, expressed in term of statistical weight
  !> (sum of the statistical weights of all the particles which have
  !> interacted with the boundary face). \ref seuilf is
  !> the limit statistical weight value, below which the contribution of the
  !> face is not taken into account in the
  !> statistics at the boundaries for post-processing.
  !> useful if \ref iensi3=1
  double precision, save ::  seuilf

  !> name of the boundary statistics, displayed in the listing
  !> and the post-processing files.
  !> useful if \ref iensi3=1.
  !> Warning: this name is also used to reference information in the restart file
  !> (\ref isuist =1). If the name of a variable is changed between two
  !> calculations, it will not be possible to read its value from the restart file
  character*50, save ::      nombrd(nvplmx)

  !> ifrlag pointer in ia to identify boundary zones and faces
  integer, save ::           iifrla

  !> \}

  !=============================================================================

  !> \defgroup visualization Visualization

  !> \addtogroup visualization
  !> \{

  !> activation (=1) or not (=0) of the recording of the particle/boundary
  !> interactions in  \ref parbor, and of the calculation of the
  !> statistics at the corresponding boundaries, for post-processing
  !> (\textit EnSight6 format).
  !> By default, the statistics are non-stationary (reset to zero at every
  !> Lagrangian iteration). They may be stationary if \ref isttio=1 (i.e.
  !> calculation of a cumulated value over time, and then calculation of an
  !> average over time or over the number of interactions with the boundary).
  integer, save ::  iensi3

  !> associates (=1) or not (=0) the variable ``velocity of the locally
  !> undisturbed fluid flow field'' with the output of particles or trajectories.
  integer, save ::  ivisv1

  !> associates (=1) or not (=0) the variable ``particle velocity''
  !> with the output of particles or trajectories.
  integer, save ::  ivisv2

  !> associates (=1) or not (=0) the variable ``residence time''
  !> with the output of particles or trajectories.
  integer, save ::  ivistp

  !> associates (=1) or not (=0) the variable ``particle diameter''
  !> with the output of particles or trajectories.
  integer, save ::  ivisdm

  !> associates (=1) or not (=0) the variable ``particle temperature''
  !> with the output of particles or trajectories.
  integer, save ::  iviste

  !> associates (=1) or not (=0) the variable ``particle mass''
  !> with the output of particles or trajectories.
  integer, save ::  ivismp

  !> associates (=1) or not (=0) the variable ``shrinking core diameter of
  !> the coal particles'' with the output of particles or trajectories.
  !> useful only if \ref iphyla = 2
  integer, save ::  ivisdk

  !> associates (=1) or not (=0) the variable ``mass of reactive coal of the
  !> coal particles'' with the output of particles or trajectories.
  !> useful only if \ref iphyla = 2
  integer, save ::  ivisch

  !> associates (=1) or not (=0) the variable ``mass of coal of the
  !> coal particles'' with the output of particles or trajectories.
  !> useful only if \ref iphyla = 2
  integer, save ::  ivisck

  ! TODO
  integer, save ::  iviswat

  !> \}

  !=============================================================================

contains

  !=============================================================================

  ! Local function to initialize Lagrangian module parameters for
  ! a given zone ii and class jj

  subroutine init_zone_class_param(ii, jj)

    use cstnum
    implicit none

    ! Arguments

    integer :: ii, jj

    ! Local variables

    integer :: kk

    ! define defaults (impossible values the user should override)

    do kk = 1, ndlaim
      iuslag(ii, jj, kk) = 0
    enddo
    iuslag(ii,jj,ijuvw) = -2
    iuslag(ii,jj,ijprtp) = -2
    iuslag(ii,jj,ijprdp) = -2
    iuslag(ii,jj,ijprpd) = -2
    do kk = 1, ndlagm
      ruslag(ii, jj, kk) = 0.d0
    enddo
    ruslag(ii,jj,iuno)  = -grand
    ruslag(ii,jj,iupt)  = -grand
    ruslag(ii,jj,ivpt)  = -grand
    ruslag(ii,jj,iwpt)  = -grand
    ruslag(ii,jj,ipoit) = -grand
    ruslag(ii,jj,idpt)  = -grand
    ruslag(ii,jj,ivdpt) = -grand
    ruslag(ii,jj,iropt) = -grand
    if (iphyla.eq.1) then
      if (itpvar.eq.1) then
        ruslag(ii,jj,itpt)  = -grand
        ruslag(ii,jj,icpt)  = -grand
        ruslag(ii,jj,iepsi) = -grand
      endif
    else if ( iphyla .eq. 2 ) then
      ruslag(ii,jj,ihpt)   = -grand
      ruslag(ii,jj,imwat)  = -grand
      ruslag(ii,jj,imcht)  = -grand
      ruslag(ii,jj,imckt)  = -grand
      ruslag(ii,jj,icpt)   = -grand
    endif

  end subroutine init_zone_class_param

  !=============================================================================

  !> \brief Initialize Lagrangian module parameters for a given zone and class

  !> \param[in]  i_cz_params  integer parameters for this class and zone
  !> \param[in]  r_cz_params  real parameters for this class and zone

  subroutine lagr_init_zone_class_param(i_cz_params, r_cz_params)  &
    bind(C, name='cs_lagr_init_zone_class_param')

    use, intrinsic :: iso_c_binding
    use cstnum
    implicit none

    ! Arguments

    integer, dimension(ndlaim)          :: i_cz_params
    double precision, dimension(ndlagm) :: r_cz_params

    ! Local variables

    integer :: ii

    ! define defaults (impossible values the user should override)

    do ii = 1, ndlaim
      i_cz_params(ii) = 0
    enddo
    i_cz_params(ijuvw) = -2
    i_cz_params(ijprtp) = -2
    i_cz_params(ijprdp) = -2
    i_cz_params(ijprpd) = -2

    do ii = 1, ndlagm
      r_cz_params(ii) = 0.d0
    enddo
    if (iphyla.eq.1) then
      if (itpvar.eq.1) then
        r_cz_params(itpt)  = -grand
        r_cz_params(icpt)  = -grand
        r_cz_params(iepsi) = -grand
      endif
    else if (iphyla .eq. 2) then
      r_cz_params(ihpt)  = -grand
      r_cz_params(imwat) = -grand
      r_cz_params(imcht) = -grand
      r_cz_params(imckt) = -grand
      r_cz_params(icpt)  = -grand
    endif

  end subroutine lagr_init_zone_class_param

  !=============================================================================

  !> \brief Define Lagrangian module parameters for a given zone and class

  !> \param[in]     class_id     id of given particle class
  !> \param[in]     zone_id      id of given boundary zone
  !> \param[in]     i_cz_params  integer parameters for this class and zone
  !> \param[in]     r_cz_params  real parameters for this class and zone

  subroutine lagr_define_zone_class_param(class_id, zone_id,         &
                                          i_cz_params, r_cz_params)  &
    bind(C, name='cs_lagr_define_zone_class_param')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, value                      :: class_id
    integer, value                      :: zone_id
    integer, dimension(ndlaim)          :: i_cz_params
    double precision, dimension(ndlagm) :: r_cz_params

    ! Local variables

    integer :: ncmax, nzmax, mcmxp, nzmxp, ii, jj, kk
    integer, allocatable, dimension(:,:,:) :: itmp
    integer, dimension(3) :: shpe
    double precision, allocatable, dimension(:,:,:) :: rtmp

    ! Allocate on first pass

    if (.not.allocated(iuslag) .or. .not.allocated(ruslag)) then
      allocate(iuslag(1,1,ndlaim))
      allocate(ruslag(1,1,ndlagm))
      call init_zone_class_param(1, 1)
    endif

    ! Reallocate arrays if required
    ! (use size margin to avoid reallocating too often, though this
    ! should only impact the first time step)

    shpe = shape(iuslag)
    ncmax = shpe(1)
    nzmax = shpe(2)

    if (class_id.gt.ncmax .or. zone_id.gt.nzmax) then

      mcmxp = ncmax
      nzmxp = nzmax

      ncmax = max(class_id, mcmxp+5)
      nzmax = max(zone_id, nzmxp+5)

      ! Save iuslag and ruslag arrays

      allocate(itmp(mcmxp,nzmxp,ndlaim))
      allocate(rtmp(mcmxp,nzmxp,ndlagm))

      do ii = 1, mcmxp
        do jj = 1, nzmxp
          do kk = 1, ndlaim
            itmp(ii,jj,kk) = iuslag(ii,jj,kk)
          enddo
          do kk = 1, ndlagm
            rtmp(ii,jj,kk) = ruslag(ii,jj,kk)
          enddo
        enddo
      enddo

      ! Reallocate iuslag and ruslag arrays

      deallocate(iuslag)
      deallocate(ruslag)
      allocate(iuslag(ncmax,nzmax,ndlaim))
      allocate(ruslag(ncmax,nzmax,ndlagm))

      ! Restore saved values, and initialize new entries

      do ii = 1, mcmxp
        do jj = 1, nzmxp
          do kk = 1, ndlaim
            iuslag(ii,jj,kk) = itmp(ii,jj,kk)
          enddo
          do kk = 1, ndlagm
            ruslag(ii,jj,kk) = rtmp(ii,jj,kk)
          enddo
        enddo
        do jj = nzmxp + 1, nzmax
          call init_zone_class_param(ii, jj)
        enddo
      enddo
      do ii = mcmxp + 1, ncmax
        do jj = 1, nzmxp
          call init_zone_class_param(ii, jj)
        enddo
      enddo

      deallocate(rtmp)
      deallocate(itmp)

    endif

    ! Now copy defined values

    do kk = 1, ndlaim
      iuslag(class_id, zone_id, kk) = i_cz_params(kk)
    enddo
    do kk = 1, ndlagm
      ruslag(class_id, zone_id, kk) = r_cz_params(kk)
    enddo

  end subroutine lagr_define_zone_class_param

  !=============================================================================

  !> \brief Return Lagrangian model status.

  !> \param[out]   model     0 without Lagrangian, 1 or 2 with Lagrangian
  !> \param[out]   restart   1 for Lagrangian restart, 0 otherwise
  !> \param[out]   frozen    1 for frozen Eulerian flow, 0 otherwise

  subroutine lagr_status(model, restart, frozen)  &
    bind(C, name='cs_lagr_status')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer(c_int), intent(out) :: model, restart, frozen

    model = iilagr
    restart = isuila
    frozen = isttio

    return

  end subroutine lagr_status

  !=============================================================================

end module lagran
