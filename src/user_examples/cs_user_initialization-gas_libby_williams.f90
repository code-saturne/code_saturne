!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file cs_user_initialization-gas_libby_williams.f90
!> \brief Libby-Williams gas example
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________


subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,nflown:nvar), propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          iel, mode, igg, izone
double precision hinit, coefg(ngazgm)
double precision sommqf, sommqt, sommq, tentm, fmelm


character*80     chaine
integer          iscal, ivar, ii
double precision valmax, valmin

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

! Control output

write(nfecra,9001)

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if ( isuite.eq.0 ) then

  ! a. Preliminary calculations

  sommqf = zero
  sommq  = zero
  sommqt = zero

  !  Deals with multiple inlets
  do izone = 1, nozapm
    sommqf = sommqf + qimp(izone)*fment(izone)
    sommqt = sommqt + qimp(izone)*tkent(izone)
    sommq  = sommq  + qimp(izone)
  enddo

  if (abs(sommq).gt.epzero) then
    fmelm = sommqf / sommq
    tentm = sommqt / sommq
  else
    fmelm = zero
    tentm = t0
  endif

  ! ----- Calculation of the Enthalpy of the gas mixture
  !       (unburned mean gas)

  if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3               &
                           .or. ippmod(icolwc).eq.5 ) then
    coefg(1) = fmelm
    coefg(2) = (1.d0-fmelm)
    coefg(3) = zero
    mode     = -1

    !   Converting the mean temperatur boundary conditions into
    !   enthalpy values
    call cothht                                                   &
    !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot   , th     , ehgazg ,                       &
        hinit  , tentm )
  endif

  do iel = 1, ncel

    ! b. Initialisation

    ! Mass fraction of Unburned (fresh) Gas
    rtp(iel,isca(iyfm))  = 0.0d0*fmelm
    ! Mean Mixture Fraction
    rtp(iel,isca(ifm))   = 0.d0*fmelm
    ! Variance of fuel Mass fraction
    rtp(iel,isca(iyfp2m)) = zero
    ! Variance of Mixture Fraction
    rtp(iel,isca(ifp2m))  = zero

    ! Covariance for NDIRAC >= 3

    if ( ippmod(icolwc).ge. 2 ) then
      rtp(iel,isca(icoyfp))   = zero
    endif

    ! Enthalpy

    if ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3             &
                             .or. ippmod(icolwc).eq.5 ) then
      rtp(iel,isca(iscalt)) = hinit
    endif

  enddo

  ! ---> Control Output of the user defined initialisation values

  write(nfecra,2000)

  do ii  = 1, nscapp
    iscal = iscapp(ii)
    ivar  = isca(iscal)
    valmax = -grand
    valmin =  grand
    do iel = 1, ncel
      valmax = max(valmax,rtp(iel,ivar))
      valmin = min(valmin,rtp(iel,ivar))
    enddo
    if ( irangp.ge.0 ) then
      call parmax(valmax)
      call parmin(valmin)
    endif

    call field_get_label(ivarfl(ivar), chaine)
    write(nfecra,2010)chaine(1:8),valmin,valmax
  enddo
  write(nfecra,2020)

endif
!< [init]

!--------
! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization: variables initialisation by user'   ,/,&
                                                                /)

 2000 format(                                                   /,&
                                                                /,&
' -----------------------------------------------------------' ,/,&
                                                                /,&
' ** INITIALISATION OF VARIABLES FOR Libby-Williams model'     ,/,&
'    --------------------------------------------------------' ,/,&
'           ONLY ONE PASS'                                     ,/,&
' ---------------------------------'                           ,/,&
'  Variable  Valeur min  Valeur max'                           ,/,&
' ---------------------------------'                             )

 2010 format(                                                     &
 2x,     a8,      e12.4,      e12.4                              )

 2020 format(                                                     &
' ---------------------------------'                           ,/)

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization
