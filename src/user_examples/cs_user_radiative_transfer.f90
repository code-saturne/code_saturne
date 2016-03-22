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
! --------
!> \file cs_user_radiative_transfer.f90
!>
!> \brief User subroutine for input of radiative transfer parameters : absorption
!> coefficient and net radiation flux.
!>
!> See \subpage cs_user_radiative_transfer for examples.
!-------------------------------------------------------------------------------

subroutine usray3 &
!================

 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                       &
   izfrdp ,                                                       &
   dt     ,                                                       &
   ck     )

!===============================================================================
! Purpose:
! --------

!> \brief Absorption coefficient for radiative module
!>
!> It is necessary to define the value of the fluid's absorption
!> coefficient Ck.
!>
!> For a transparent medium, the coefficient should be set to 0.d0.
!>
!> In the case of the P-1 model, we check that the optical length is at
!> least of the order of 1.
!>
!> Caution:
!>
!> For specific physics (Combustion, coal, ...),
!> it is Forbidden to define the absorption coefficient here.
!>
!>   See subroutine \ref ppcabs.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     itypfb        boundary face types
!> \param[in]     izfrdp        zone number for boundary faces
!> \param[in]     dt            time step (per cell)
!> \param[out]    ck            medium's absorption coefficient
!>                              (zero if transparent)
!______________________________________________________________________________!

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

!< [arg_1]

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision dt(ncelet)

double precision ck(ncelet)

!< [arg_1]

!< [loc_var_dec_1]

! Local variables

integer          iel, ifac, iok
double precision vv, sf, xlc, xkmin, pp

!< [loc_var_dec_1]

!===============================================================================

!< [abso_coeff]

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

    do ifac = 1, nfabor
      sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
    endif

    do iel = 1, ncel
      vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
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

!< [abso_coeff]

!< [format_1]

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
'@      optical length is XKmin = ', e11.4                     ,/,&
'@    This value is not reached for ', e11.4,'%               ',/,&
'@      of the meshe''s cells.                                ',/,&
'@    The percentage of mesh cells for which we allow this    ',/,&
'@      condition not to be respected is set by default in    ',/,&
'@      cs_user_parameters.f90 to xnp1mx = ', e11.4,'%        ',/,&
'@                                                            ',/,&
'@    The calculation is interrupted.                         ',/,&
'@                                                            ',/,&
'@    Check the values of the absorption coefficient Ck       ',/,&
'@      in ppcabs, usray3 or the thermochemistry data file.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!< [format_1]

end subroutine usray3


!===============================================================================

subroutine usray5 &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   izfrdp ,                                                       &
   dt     ,                                                       &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   tparoi , qincid , flunet , xlam   , epa    , eps     ,  ck   )

!===============================================================================
!  Purpose:
!  --------

!> \brief Compute the net radiation flux:
!>
!> The density of net radiation flux must be calculated
!> consistently with the boundary conditions of the intensity.
!> The density of net flux is the balance between the radiative
!> emiting part of a boudary face (and not the reflecting one)
!> and the radiative absorbing part.
!>
!> The provided example is consistent with the proposed boundary
!> conditions for the intensity.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name             role                                        !
!______________________________________________________________________________!
!> \param[in]     nvar            total number of variables
!> \param[in]     nscal           total number of scalars
!> \param[in]     itypfb          boundary face types
!> \param[out]    izfrdp          boundary faces -> zone number
!> \param[in]     dt              time step (per cell)
!> \param[out]    coefap, coefbp  boundary conditions for intensity or P-1 model
!>                cofafp,cofbfp
!> \param[in]     tparoi          inside current wall temperature (K)
!> \param[in]     qincid          radiative incident flux  (W/m2)
!> \param[out]    flunet          net flux (W/m2)
!> \param[out]    xlamp           conductivity (W/m/K)
!> \param[out]    epap            thickness (m)
!> \param[out]    epsp            emissivity (>0)
!> \param[in]     ck              absoprtion coefficient
!>                                gaz-particules de charbon
!______________________________________________________________________________!

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

!< [arg_2]

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          icodcl(nfabor,nvar)
integer          izfrdp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision tparoi(nfabor), qincid(nfabor)
double precision xlam(nfabor), epa(nfabor)
double precision eps(nfabor), flunet(nfabor)
double precision ck(ncelet)

!< [arg_2]

!< [loc_var_dec_2]

! Local variables

integer          ifac, iok

!< [loc_var_dec_2]

!===============================================================================

!< [init]

!===============================================================================
! 0 - Initialization
!===============================================================================

! Stop indicator (forgotten boundary faces)
iok = 0

!< [init]

!< [net_flux]

!===============================================================================
! 1 - Net flux dendity for the boundary faces
!     The provided examples are sufficient in most of cases.
!===============================================================================

!    If the boundary conditions given above have been modified
!      it is necessary to change the way in which density is calculated from
!      the net radiative flux consistently.
!    The rule is:
!      the density of net flux is a balance between the emitting energy from a
!      boundary face (and not the reflecting energy) and the absorbing radiative
!      energy. Therefore if a wall heats the fluid by radiative transfer, the
!      net flux is negative

do ifac = 1, nfabor

  ! Wall faces
  if (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug) then

    flunet(ifac) = eps(ifac) *(qincid(ifac) - stephn*tparoi(ifac)**4)

  ! Symmetry
  else if (itypfb(ifac).eq.isymet) then

    flunet(ifac) = zero

  ! Inlet/Outlet
  else if (itypfb(ifac).eq.ientre             .or. &
           itypfb(ifac).eq.i_convective_inlet .or. &
           itypfb(ifac).eq.isolib             .or. &
           itypfb(ifac).eq.ifrent                 ) then

    if (iirayo.eq.1) then
      flunet(ifac) = qincid(ifac) -pi*coefap(ifac)
    else if (iirayo.eq.2) then
      flunet(ifac)= 0.d0
    endif

  ! Stop if there are forgotten faces
  else

    write (nfecra,2000) ifac,izfrdp(ifac),itypfb(ifac)
    iok = iok + 1

  endif

enddo


if (iok.ne.0) then
  write (nfecra,2100)
  call csexit (1)
endif

!< [net_flux]

!< [format_2]

! -------
! Formats
! -------

 2000 format( &
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

 2100 format( &
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
'@    Please verify subroutine usray5.                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!< [format_2]

!----
! End
!----

end subroutine usray5
