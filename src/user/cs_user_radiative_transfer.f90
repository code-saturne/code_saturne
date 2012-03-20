!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine usray3 &
!================

 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                       &
   izfrdp ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   ck     )

!===============================================================================
! Purpose:
! --------

! Absorption coefficient for radiative module
! ----------------------

! It is necessary to define the value of the fluid's absorption
! coefficient Ck.

! For a transparent medium, the coefficient should be set to 0.d0.

! In the case of the P-1 model, we check that the optical length is at
! least of the order of 1.

! Caution:
! ========
!   For specific physics (Combustion, coal, ...),

!   it is Forbidden to define the absorption coefficient here.
!         =========

!   See subroutine ppcabs.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrdp(nfabor    ! ia ! <-- ! zone number for boundary faces                 !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ck(ncelet)       ! ra ! --> ! medium's absorption coefficient                !
!                  !    !     ! (zero if transparent)                          !
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
use cstnum
use parall
use period
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)


double precision ck(ncelet)

! Local variables

integer          iel, ifac, iok
double precision vv, sf, xlc, xkmin, pp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

!===============================================================================
! 0.  This test allows the user to be certain his version if the subroutine
!     is being used, and not that from the library.
!===============================================================================

if (iihmpr.eq.1) then
   return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET RAYONNEMENT                           ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usray3 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0 - Memory management
!===============================================================================


! Stop flag (to determine if faces were forgotten)
iok = 0

!===============================================================================
! Absorption coefficient of the medium (m-1)

! In the case of specific physics (gas/coal/fuel combustion, elec)

! Ck must not be defined here
!============
! (it is determined automatically, possibly from the parametric file)

! In other cases, Ck must be defined (it is zero by default)
!===============================================================================

if (ippmod(iphpar).le.1) then

  do iel = 1, ncel
    ck(iel) = 0.d0
  enddo

  !--> P1 model: standard control of absorption coefficient values.
  !              this coefficient must ensure an optical length
  !              at least of the order of 1.

  if (iirayo.eq.2) then
    sf = 0.d0
    vv = 0.d0

    ! Compute characteristic length of calculation domain

    do ifac = 1,nfabor
      sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
      !==========
    endif

    do iel = 1,ncel
      vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
      !==========
    endif

    xlc = 3.6d0 * vv / sf

    ! Clipping for variable CK

    xkmin = 1.d0 / xlc

    iok = 0

    do iel = 1,ncel
      if (ck(iel).lt.xkmin) then
        iok = iok + 1
      endif
    enddo

    ! Stop at the end of the time step if the optical thickness is too big
    ! (istpp1 = 1 allows stopping cleanly at the end of the current time step).
    pp = xnp1mx/100.0d0
    if (dble(iok).gt.pp*dble(ncel)) then
      write(nfecra,3000) xkmin, dble(iok)/dble(ncel)*100.d0, xnp1mx
      istpp1 = 1
      ! call csexit(1)
      !==========
    endif
  endif

endif

! -------
! Formats
! -------

 3000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:    P1 radiation approximation (usray3)         ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The optical length of the semi-transparent medium       ',/,&
'@      must be at least of the order of one to be in the     ',/,&
'@      domain of validity of the P-1 approximation.          ',/,&
'@    This does not seem to be the case here.                 ',/,&
'@                                                            ',/,&
'@    The minimum absorption coefficient to ensure this       ',/,&
'@      optical length is XKmin = ', e10.4                     ,/,&
'@    This value is not reached for ', e10.4,'%               ',/,&
'@      of the meshe''s cells.                                ',/,&
'@    The percentage of mesh cells for which we allow this    ',/,&
'@      condition not to be rspected is set by default in     ',/,&
'@      usini1 to xnp1mx = ', e10.4,'%                        ',/,&
'@                                                            ',/,&
'@    The calculation is interrupted.                         ',/,&
'@                                                            ',/,&
'@    Check the values of the absorption coefficient Ck       ',/,&
'@      in ppcabs, usray3 or the thermochemistry data file.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine


!===============================================================================


subroutine usray4 &
!================

 ( nvar   , nscal  ,                                              &
   mode   ,                                                       &
   itypfb ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   tparop , hparop , tempk  )

!===============================================================================
! Purpose:
! --------

! User subroutine for input of radiative transfer parameters:

!   Temperature <--> enthalpy convertion
!   Usefull if the thermal scalar is an enthalpy.

!   PRECAUTIONS: ENTHALPY MUST BE CONVERTED IN KELVIN TEMPERATURE

!   Warning: it is forbidden to modify MODE in this subroutine

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! mode             ! i  ! <-- ! convertion mode                                !
!                  !    !     ! mode = 1 enthaly -> temperature                !
!                  !    !     ! mode =-1 temperature -> enthaly                !
! itypfb           ! ia ! <-- ! boundary face types                            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! tparop(nfabor)   ! i  ! <-- ! temperature in kelvin for wall boundary faces  !
! hparop(nfabor)   ! i  ! --> ! enthalpy for wall boundary faces               !
! tempk(ncelet)    ! i  ! --> ! temperature in kelvin                          !
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
use cstnum
use parall
use period
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mode

integer          itypfb(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision tempk(ncelet)
double precision tparop(nfabor), hparop(nfabor)

! Local variables

integer          iel , ifac , iscal

!===============================================================================

!===============================================================================
! 1 - INITIALISATIONS GENERALES
!===============================================================================


iscal = iscalt

!===============================================================================
!  2.1 - Tempearature in kelvin for cells
!===============================================================================

!---> enthalpy -> temperature convertion (MODE =  1)
!     -----------------------------------------------


if (mode.eq.1) then

  do iel = 1,ncel
    call usthht(mode,rtpa(iel,isca(iscal)),tempk(iel))
  enddo

endif


!===============================================================================
!  2.2 - Enthalpy for wall boundary faces
!===============================================================================

!---> Temperature -> enthalpy (MODE = -1)
!     -----------------------------------


if (mode.eq.-1) then

  do ifac = 1,nfabor

    if (itypfb(ifac).eq.iparoi .or.                               &
        itypfb(ifac).eq.iparug )then

      call usthht(mode,hparop(ifac),tparop(ifac))

    endif

  enddo

endif

!----
! END
!----

return
end subroutine


!===============================================================================


subroutine usray5 &
!================

 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                       &
   izfrdp ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   tparoi , qincid , flunet , xlam   , epa    , eps     ,  ck   )

!===============================================================================
!  Purpose:
!  --------

! User subroutine for input of radiative transfer parameters.

! This subroutine is calle twice. The first time is for boundary conditions.
! The second time is for the net radiation flux computation.


!  1. Fisrt call (IAPPEL = 1)
!  ==========================

!    1.1 Boundary conditions fot the radiative intensity (DO model)
!    --------------------------------------------------------------

!       The array COFRUA store the intensity for each boundary faces,
!         depending of the natur of the boundary (Dirichlet condition).
!       The intensity of radiation is defined as the rate of emitted
!         energy from unit surface area through unit solid angle.

!       For example:


! 1/ Gray wall: isotropic radiation field.
!                                    4
!                      eps.sig.tparoi         (1-eps).qincid
!        cofrua   =    --------------    +    --------------
!                            pi                     pi
!  wall intensity     wall emission           reflecting flux.

!     (eps=1: black wall; eps=0: reflecting wall)


! 2/ Free boundary: entering intensity is fixed to zero

!        cofrua   =   0.D0

!    (if the user has more information, he can do something better)



!    1.2 Boundary conditions fior the P-1 model
!    ------------------------------------------



!  2. Second call(IAPPEL = 2)
!  ==========================

!      The density of net radiation flux must calculated
!        consistently with the boundary conditions of the intensity.
!      The density of net flux is the balance between the radiative
!      emiting part of a boudary face (and not the reflecting one)
!      and the radiative absorbing part.

!      The provided example is consistently with the proposed boundary
!      conditions for the intensity.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iappel           ! i  ! <-- ! current subroutine call number                 !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrdp(nfabor)   ! ia ! --> ! boundary faces -> zone number                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! cofrua, cofrub   ! ra ! --> ! boundary conditions for intensity or P-1 model !
!  (nfabor, *)     !    !     !                                                !
! tparoi(nfabor)   ! ra ! <-- ! inside current wall temperature (K)            !
! qincid(nfabor)   ! ra ! <-- ! radiative incident flux  (W/m2)                !
! flunet(nfabor)   ! ra ! --> ! net flux (W/m2)                                !
! xlamp(nfabor)    ! ra ! --> ! conductivity (W/m/K)                           !
! epap(nfabor)     ! ra ! --> ! thickness (m)                                  !
! epsp(nfabor)     ! ra ! --> ! emissivity (>0)                                !
! ck(ncelet)       ! ra ! <-- ! absoprtion coefficient                         !
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
use cstnum
use parall
use period
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision cofrua(nfabor), cofrub(nfabor)

double precision tparoi(nfabor), qincid(nfabor)
double precision xlam(nfabor), epa(nfabor)
double precision eps(nfabor), flunet(nfabor)
double precision ck(ncelet)

! Local variables

integer          ifac, iok
double precision unspi, xit, distbf

!===============================================================================

!===============================================================================
! 0 - Initialization
!===============================================================================


! Stop indicator (forgotten boundary faces)
iok = 0

unspi = 1.d0/pi

!==============================================================================
!  1. First call
!  =============
!===============================================================================

if (iappel.eq.1) then


!===============================================================================
!  1.1 - Boundary conditions:
!        DO model: COFRUA msut be filled with the intensity
!        P-1 model: COFRUA and COFRUB msut be filled
!      The provided examples are sufficient in most of cases.
!===============================================================================




!      A - DO model
!      ^^^^^^^^^^^^




  if (iirayo.eq.1) then

    do ifac = 1,nfabor

!      1.1.1 - Symmetry :
!              ----------
!          Reflecting boundary conditions ( EPS=0 )
!          ----------------------------------------

      if (itypfb(ifac).eq.isymet) then

        cofrua(ifac) = qincid(ifac) * unspi


!      1.1.2 - Inlet/Outlet face: entering intensity fixed to zero
!              (WARNING: the treatment is different from than of P-1 model)
!          -------------------------------------------------

      else if (itypfb(ifac).eq.ientre                             &
          .or. itypfb(ifac).eq.isolib) then

        cofrua(ifac) = epzero


!      1.1.3. - Wall boundary face: calculaed intensity
!               ---------------------------------------

      else if (itypfb(ifac).eq.iparoi                             &
          .or. itypfb(ifac).eq.iparug) then

        cofrua(ifac)  = eps(ifac)*stephn*(tparoi(ifac)**4)*unspi  &
                          + (1.d0-eps(ifac))* qincid(ifac)*unspi

      else

!      1.1.4 - Stop if there are forgotten faces
!              ---------------------------------

!           ==============================================

!             Don't skip this test

!           ==============================================

        write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
        iok = iok + 1
      endif

    enddo




!   B - P-1 model
!   ^^^^^^^^^^^^^





  else if (iirayo.eq.2) then

    do ifac = 1,nfabor

!      1.1.1 - Symmetry or reflecting wall (EPS = 0) :
!              zero flux
!              ----------------------------------------

      if (itypfb(ifac).eq.isymet     .or.                         &
         ((itypfb(ifac).eq.iparoi.or.                             &
           itypfb(ifac).eq.iparug).and.eps(ifac).eq.0d0)) then

        cofrua(ifac) = 0.d0
        cofrub(ifac) = 1.d0


!      1.1.2 - Inlet/Outlet faces: zero flux
!              (WARNING: the treatment is different from than of DO model)
!              ----------------------------------------------------------

      else if (itypfb(ifac).eq.ientre                             &
          .or. itypfb(ifac).eq.isolib) then

        cofrua(ifac) = 0.d0
        cofrub(ifac) = 1.d0


!      1.1.3 - Wall boundary faces
!              -------------------

      else if (itypfb(ifac).eq.iparoi .or.                        &
               itypfb(ifac).eq.iparug ) then

        distbf = distb(ifac)

        xit = 1.5d0 *distbf *ck(ifabor(ifac))                     &
            * (2.d0 /(2.d0-eps(ifac)) -1.d0)

        cofrub(ifac) = 1.d0 / (1.d0 + xit)
        cofrua(ifac) = xit * tparoi(ifac)**4 * cofrub(ifac)

      else

!      1.1.4 - Stop if there are forgotten faces
!              ---------------------------------

!           ==============================================

!             Don't skip this test

!           ==============================================

        write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1
    endif

  enddo

  endif

  if (iok.ne.0) then
    write (nfecra,1100)
    call csexit (1)
    !==========
  endif

!===============================================================================
!  2 - Second call
!  ================
!===============================================================================

else if (iappel.eq.2) then

!===============================================================================
!  2.1 - Net flux dendity for the boundary faces
!      The provided examples are sufficient in most of cases.
!===============================================================================

!    If the boundary conditions given above have been modified
!      it is necessary to change the way in which density is calculated from
!      the net radiative flux consistently.
!    The rule is:
!      the density of net flux is a balance between the emitting energy from a
!      boundary face (and not the reflecting energy) and the absorbing radiative
!      energy. Therefore if a wall heats the fluid by radiative transfer, the
!      net flux is negative


  do ifac = 1,nfabor

    if (itypfb(ifac).eq.iparoi .or.                               &
        itypfb(ifac).eq.iparug) then

!      2.1.1 - Wall faces
!              ----------

      flunet(ifac) =                                              &
      eps(ifac) *(qincid(ifac) - stephn*tparoi(ifac)**4)


!      2.1.2 - Symmetry :
!              ----------

    else if (itypfb(ifac).eq.isymet) then

      flunet(ifac)= zero


!      2.1.3 - Inlet/Outlet
!              ------------

    else if (itypfb(ifac).eq.ientre                               &
        .or. itypfb(ifac).eq.isolib) then

      if (iirayo.eq.1) then

      flunet(ifac)= qincid(ifac) -pi*cofrua(ifac)

      else if (iirayo.eq.2) then

        flunet(ifac)= 0.d0

      endif


!      2.1.4 - Stop if there are forgotten faces
!              ---------------------------------
    else

!           ==============================================

!             Don't skip this test

!           ==============================================

      write (nfecra,2000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1

    endif

  enddo


  if (iok.ne.0) then
    write (nfecra,2100)
    call csexit (1)
    !==========
  endif


endif

! -------
! FORMATS
! -------

 1000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (usray5)                    ',/,&
'@    ========                                                ',/,&
'@              Boundary conditions non inquiries             ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Nature = ',I10         )

 1100 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (usray5)                    ',/,&
'@    ========                                                ',/,&
'@    Boundary conditions are unknown for some faces          ',/,&
'@                                                            ',/,&
'@    The calculation stops.                                  ',/,&
'@                                                            ',/,&
'@    Please verify subroutine usray5.                        ',/, &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (usray5)                    ',/,&
'@    ========                                                ',/,&
'@              Net flux calculation non inquiries            ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Nature = ',I10         )

 2100 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (usray5)                    ',/,&
'@    ========                                                ',/,&
'@    Net radiation flux is unknown for some faces            ',/,&
'@                                                            ',/,&
'@    The calculation stops.                                  ',/,&
'@                                                            ',/,&
'@    Please verify subroutine usray5.                        ',/, &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)



end subroutine
