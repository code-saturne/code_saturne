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

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision dt(ncelet)

double precision ck(ncelet)

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
'@ @@ WARNING: stopped by radiative module                    ',/,&
'@    ========                                                ',/,&
'@     User subroutine usray3 must be completed               ',/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! -------
! Formats
! -------

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

!===============================================================================

! -------
! Formats
! -------

!----
! End
!----

end subroutine usray5
