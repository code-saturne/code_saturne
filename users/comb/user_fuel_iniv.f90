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

subroutine user_fuel_iniv &
!========================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

!===============================================================================
! PURPOSE  :
! --------

! INITIALISATION OF TRANSPORTED VARIABLES
!    EXTENDED PHYSICS : Heavy Fuel Oil Combustion
!    similar to  USINIV.F
!
! This routine is called at the beginning of every computation
!  (new or continuation) before the time loop
!
! This routine initialize or modify (if continuation)
!  values of transported variables and of the time step
!
! The exemple is ... default value

!
! Physical properties are stored in PROPCE(cell center)
!  PROPFA(inner face) and PROPFB(boundary face)
!  e.g.
!  PROPCE(IEL, IPPROC(IROM  )) is ROM(IEL) mean density kg/m3
!  PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) is FLUMAS(IFAC,IVAR) convective flux
!                                                        of variable IVAR
!  PROPFB(......                      .................................
!
! Physical properties (ROM, VISCL, CP ...) are computed in PPPHYV
!  not to be modified here
!

!   All cells can be identified by using the subroutine 'getcel'.
!    Syntax of getcel:
!     getcel(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.

!       string may contain:
!       - references to colors (ex.: 1, 8, 26, ...
!       - references to groups (ex.: inlet, group1, ...)
!       - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!       These criteria may be combined using logical operators
!       ('and', 'or') and parentheses.
!       Example: '1 and (group2 or group3) and y < 1' will select boundary
!       faces of color 1, belonging to groups 'group2' or 'group3' and
!       with face center coordinate y less than 1.
!
!   All boundary faces may be identified using the 'getfbr' subroutine.
!    Syntax of getfbr:
!     getfbr(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.
!
!     string may contain:
!     - references to colors (ex.: 1, 8, 26, ...
!     - references to groups (ex.: inlet, group1, ...)
!     - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!     These criteria may be combined using logical operators
!     ('and', 'or') and parentheses.
!
!   All internam faces may be identified using the 'getfac' subroutine.
!    Syntax of getfac:
!     getfac(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.
!
!     string may contain:
!     - references to colors (ex.: 1, 8, 26, ...
!     - references to groups (ex.: inlet, group1, ...)
!     - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!     These criteria may be combined using logical operators
!     ('and', 'or') and parentheses.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! LOCAL VARIABLES

integer          iel, ige, mode, icla
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision xkent, xeent, d2s3
double precision dmas , wmco2 , wmh2o , wmn2 , wmo2

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. CONTROL PRINT
!===============================================================================

write(nfecra,9001)

!===============================================================================
! 1.  LOCAL VARIABLES INITIALISATION
!===============================================================================

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION OF TRANSPORTED VARIABLES
!      RONLY IF THE COMPUTATION IS NOT A CONTINUATION
!===============================================================================

if ( isuite.eq.0 ) then

! --> Initialisation of k and epsilon (exemple)

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur.eq.2) then

    do iel = 1, ncel
      rtp(iel,ik)  = xkent
      rtp(iel,iep) = xeent
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      rtp(iel,ir11) = d2s3*xkent
      rtp(iel,ir22) = d2s3*xkent
      rtp(iel,ir33) = d2s3*xkent
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
      rtp(iel,iep)  = xeent
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iep)  = xeent
      rtp(iel,iphi) = d2s3
      rtp(iel,ifb)  = 0.d0
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iomg) = xeent/cmu/xkent
    enddo

  endif

! --> All the computation domain is initialised with air at TINITK
!                   ================================================

! ---- Computation of H1INIT and  H2INIT

  t1init = 1000.d0
  t2init = 1000.d0

! ------ Transported variables for droplets

  h2init = h02fol +  cp2fol*(t2init-trefth)

  do icla = 1, nclafu
    do iel = 1, ncel
      rtp(iel,isca(iyfol(icla))) = zero
      rtp(iel,isca(ing(icla  )))  = zero
      rtp(iel,isca(ih2(icla  )))  = zero
    enddo
  enddo

! ------ Transported variables for the mix (droplets and carrying gases)

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo
!  On considere l'oxydant 1
  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  mode = -1
  call cs_fuel_htconvers1(mode,h1init,coefe,t1init)
 !============================

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = h1init
  enddo

! ------ Transported variables for gaseous mixture
!        (passive scalars, variance, reactive species)

  do iel = 1, ncel
    rtp(iel,isca(ifvap )) = 0.d0
    rtp(iel,isca(if7m  )) = zero
    rtp(iel,isca(ifvp2m)) = zero
    if ( ieqco2 .ge. 1 ) then
      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      rtp(iel,isca(iyco2)) = oxyco2(ioxy)*wmco2/dmas
    endif
    if ( ieqnox .eq. 1 ) then
      rtp(iel,isca(iyhcn)) = zero
      rtp(iel,isca(iyno )) = zero
      rtp(iel,isca(ihox)) = h1init
    endif
  enddo

endif


!----
! FORMATS
!----

 9001 format(                                                     &
'                                                             ',/,&
'  usfuiv : Variables Initialisation for FUel by the USer     ',/,&
'                                                             ',/)

!----
! END
!----
return
end subroutine
