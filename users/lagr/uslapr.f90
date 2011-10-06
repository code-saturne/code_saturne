!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine uslapr &
!================

 ( idvar  , iepart , izone  , iclass ,                            &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   xxpart , yypart , zzpart ,                                     &
   tvpart , uupart , vvpart , wwpart , ddpart , ttpart  ,         &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!   User subroutine for the boundary conditions associated to
!   the particles (inlet and treatment of the other boundaries)
!
!   It allows to impose the values of the velocity, the diameter
!   and the temperature for the treated particle.
!
!   if idvar = 0 ==> the volume fraction is retrieved
!   if idvar = 1 ==> the 3 components of the velocity are retrieved
!   if idvar = 2 ==> the diameter is retrieved
!   if idvar = 3 ==> the temperature is retrieved


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idvar            ! i  ! <-- ! type of the value(s) ta calculate              !
! iepart           ! i  ! <-- ! number of the particle cell                    !
! izone            ! i  ! <-- ! number of the particle zone                    !
! iclass           ! i  ! <-- ! number of the particle class                   !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the             !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor)   ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! xxpart           !  r ! <-- ! x-coordinate of the particle                   !
! yypart           !  r ! <-- ! y-coordinate of the particle                   !
! zzpart           !  r ! <-- ! z-coordinate of the particle                   !
! tvpart           !  r ! <-- ! value of the volume fraction                   !
! uupart           !  r ! <-- ! x-component of particle velocity               !
! vvpart           !  r ! <-- ! y-component of particle velocity               !
! wwpart           !  r ! <-- ! z-component of particle velocity               !
! ddpart           !  r ! <-- ! particle diameter                              !
! ttpart           !  r ! <-- ! particle temperature                           !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous          !
! (ncelet,*)       !    !     ! time step                                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments


integer          idvar  , iepart , izone  , iclass

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)

double precision xxpart , yypart , zzpart
double precision tvpart , uupart , vvpart , wwpart
double precision ddpart , ttpart

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)

! Local variables


double precision pis6

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
  !==========
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     MODULE LAGRANGIEN :                                    ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslapr DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Memory management
!===============================================================================


!===============================================================================
! 2. Initialization
!===============================================================================

pis6 = pi / 6.d0

!===============================================================================
! 3. Profile for the volume fraction
!===============================================================================

if (idvar .eq. 0) then

  tvpart = 0.01

endif

!===============================================================================
! 4. Velocity profile
!===============================================================================

if (idvar .eq. 1) then

  uupart = 1.d0
  vvpart = 0.d0
  wwpart = 0.d0

endif

!===============================================================================
! 5. Diameter profile
!===============================================================================

if (idvar .eq. 2) then

  ddpart = 50.d-6

endif


!===============================================================================
! 6. Temperature profile
!===============================================================================

if (idvar .eq. 3) then

  ttpart = 20.d0

endif

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
