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


subroutine raycll &
!================

 ( itypfb ,                                                       &
   izfrdp ,                                                       &
   rtp    , rtpa   ,                                              &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   tparoi , qincid , eps    ,  ck   , ckmel  )

!===============================================================================
!  Purpose:
!  --------

!    1. Boundary conditions fot the radiative intensity (DO model)
!    --------------------------------------------------------------

!       The array coefap stores the intensity for each boundary faces,
!         depending of the natur of the boundary (Dirichlet condition).
!       The intensity of radiation is defined as the rate of emitted
!         energy from unit surface area through unit solid angle.

!       For example:


! 1/ Gray wall: isotropic radiation field.
!                                    4
!                      eps.sig.tparoi         (1-eps).qincid
!        coefap   =    --------------    +    --------------
!                            pi                     pi
!  wall intensity     wall emission           reflecting flux.

!     (eps=1: black wall; eps=0: reflecting wall)


! 2/ Free boundary: entering intensity is fixed to zero

!        coefap   =   0.d0

!    (if the user has more information, he can do something better)



!    2. Boundary conditions fior the P-1 model
!    -----------------------------------------
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iappel           ! i  ! <-- ! current subroutine call number                 !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrdp(nfabor)   ! ia ! --> ! boundary faces -> zone number                  !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! coefap, coefbp   ! ra ! --> ! boundary conditions for intensity or P-1 model !
!  cofafp,cofbfp   !    !     !                                                !
! tparoi(nfabor)   ! ra ! <-- ! inside current wall temperature (K)            !
! qincid(nfabor)   ! ra ! <-- ! radiative incident flux  (W/m2)                !
! xlamp(nfabor)    ! ra ! --> ! conductivity (W/m/K)                           !
! epap(nfabor)     ! ra ! --> ! thickness (m)                                  !
! epsp(nfabor)     ! ra ! --> ! emissivity (>0)                                !
! ck(ncelet)       ! ra ! <-- ! absoprtion coefficient                         !
! ckmel(ncelet)    ! tr ! <-- ! coeff d'absorption du melange                  !
!                  !    !     !   gaz-particules de charbon                    !
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

integer          iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision tparoi(nfabor), qincid(nfabor)
double precision eps(nfabor)
double precision ck(ncelet)
double precision ckmel(ncelet)

! Local variables

integer          iel, ifac, iok
double precision unspi, xit, distbf, pimp, hint, qimp, coefmn, xlimit, cfl

!===============================================================================

!===============================================================================
! 0 - Initialization
!===============================================================================

! Stop indicator (forgotten boundary faces)
iok = 0

unspi = 1.d0/pi

!--> Initialization to a non-admissible value for testing after raycll
do ifac = 1, nfabor
  coefap(ifac) = -grand
enddo
!
!===============================================================================
! 1. Boundary conditions for DO model
!    coefap must be filled with the intensity
!===============================================================================

if (iirayo.eq.1) then

  do ifac = 1, nfabor

    ! Dirichlet Boundary Conditions
    !------------------------------

    hint = 0.d0

    ! 1.1 - Symmetry :
    !       ----------
    !   Reflecting boundary conditions ( EPS=0 )
    !   ----------------------------------------

    if (itypfb(ifac).eq.isymet) then

      pimp = qincid(ifac) * unspi

    ! 1.2 - Inlet/Outlet face: entering intensity fixed to zero
    !       (WARNING: the treatment is different from than of P-1 model)
    !   -------------------------------------------------

    else if (itypfb(ifac).eq.ientre .or. itypfb(ifac).eq.isolib) then

      pimp = epzero

    ! 1.3. - Wall boundary face: calculaed intensity
    !        ---------------------------------------

    else if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

      pimp = eps(ifac) * stephn*(tparoi(ifac)**4)*unspi  &
           + (1.d0-eps(ifac))* qincid(ifac)*unspi

    else

      ! 1.4 - Stop if there are forgotten faces
      !       ---------------------------------

      write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1
    endif

    call set_dirichlet_scalar &
         !====================
       ( coefap(ifac), cofafp(ifac),             &
         coefbp(ifac), cofbfp(ifac),             &
         pimp        , hint        , rinfin )

  enddo

!===============================================================================
! 2. Boundary conditions for P-1 model:
!    coefap and coefbp msut be filled
!===============================================================================

else if (iirayo.eq.2) then

  do ifac = 1, nfabor

    iel = ifabor(ifac)
    hint = 1.d0 / (ckmel(iel)*distb(ifac))

    ! 2.1 - Symmetry or reflecting wall (EPS = 0) :
    !       zero flux
    !       ----------------------------------------

    if (itypfb(ifac).eq.isymet.or.                             &
      ((itypfb(ifac).eq.iparoi.or.                             &
        itypfb(ifac).eq.iparug).and.eps(ifac).eq.0d0)) then

      ! Neumann Boundary Condition
      !---------------------------

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )


    ! 2.2 - Inlet/Outlet faces: zero flux
    !       (WARNING: the treatment is different from than of DO model)
    !       ----------------------------------------------------------

    else if (itypfb(ifac).eq.ientre .or. itypfb(ifac).eq.isolib) then

      ! Neumann Boundary Condition
      !---------------------------

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )


    ! 2.3 - Wall boundary faces
    !       -------------------

    else if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug ) then

      distbf = distb(ifac)

      xit = 1.5d0 *distbf *ck(iel) * (2.d0 /(2.d0-eps(ifac)) -1.d0)
      cfl = 1.d0/xit

      ! Convective Boundary Condition
      !------------------------------

      pimp = tparoi(ifac)**4

      call set_convective_outlet_scalar &
           !===========================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           pimp        , cfl         , hint )

    else

    ! 2.4 - Stop if there are forgotten faces
    !       ---------------------------------

      write (nfecra,1000) ifac,izfrdp(ifac),itypfb(ifac)
      iok = iok + 1
    endif

  enddo

endif

if (iok.ne.0) then
  write (nfecra,1100)
  call csexit (1)
endif

!---> Check luminance boundary conditions arrays

!     Attention : dans le cas de l'approx P-1 la valeur de coefap peut
!     etre grande (de l'ordre de tparoi**4), d'ou la valeur de coefmn

iok = 0
xlimit = -grand*0.1d0
! GRAND n'est pas assez grand...
coefmn = rinfin

do ifac = 1, nfabor
  if (coefap(ifac).le.xlimit) then
    iok = iok + 1
    coefmn = min(coefmn, coefap(ifac))
    write(nfecra,3000)ifac,izfrdp(ifac),itypfb(ifac)
  endif
enddo

if (iok.ne.0) then
  write(nfecra,1100)
  call csexit (1)
endif

! -------
! Formats
! -------

 1000 format( &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (raycll)                    ',/,&
'@    ========                                                ',/,&
'@              Boundary conditions non inquiries             ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Nature = ',I10         )

 1100 format( &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: Radiative transfer (raycll)                    ',/,&
'@    ========                                                ',/,&
'@    Boundary conditions are unknown for some faces          ',/,&
'@                                                            ',/,&
'@    The calculation stops.                                  ',/,&
'@                                                            ',/,&
'@    Please verify subroutine raycll.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

3000 format( &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@                CONDITIONS AUX LIMITES MAL RENSEIGNEES      ',/,&
'@                                                            ',/,&
'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )

!----
! End
!----

end subroutine
