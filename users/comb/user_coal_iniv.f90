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

subroutine user_coal_iniv                               &
!================
 ( nvar   , nscal  ,                                            &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

!===============================================================================
! PURPOSE  :
! --------
! Initialisation of transported variables for extended physics :
!  pulverised coal combustion
!   similar to USINIV
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
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
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
use ppincl
use ppcpfu
use cs_coal_incl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! VARIABLES LOCALES

integer          iel, ige, mode, icla, icha
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3
double precision wmh2o,wmco2,wmn2,wmo2,dmas

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Control Print
!===============================================================================

write(nfecra,9001)

!===============================================================================
! 1.  Local variables initialisation
!===============================================================================

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. Unknows initialisation
!      Only if the computation is'nt a continuation
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

! --> All the domain is filled with the fisrt oxidizer at TINITK
!                   ================================================

! ---- Computation of H1INIT and H2INIT

  t1init = t0
  t2init = t0

! ------ Transported variables for the solid phase
!         initialy lacking

  do icla = 1, nclacp
    icha = ichcor(icla)

    do iel = 1, ncel

      rtp(iel,isca(ixch(icla))) = zero
      rtp(iel,isca(ixck(icla))) = zero
      rtp(iel,isca(inp(icla) )) = zero
      rtp(iel,isca(ih2(icla)))  = zero

      if ( ippmod(icp3pl) .eq. 1 ) then
        rtp(iel,isca(ixwt(icla))) = zero
      endif

    enddo
  enddo

! ------ Transported variables for the mix (solid+carrying gas)²

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

!       Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
!       the composition of the fisrt oxidiser is taken in account

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

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  mode = -1
  call cs_coal_htconvers1(mode,h1init,coefe,f1mc,f2mc,t1init)
 !============================

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = h1init
  enddo

! ------ Transported variables for the mix (passive scalars, variance)

  do icha = 1, ncharb
    do iel = 1, ncel
      rtp(iel,isca(if1m(icha))) = zero
      rtp(iel,isca(if2m(icha))) = zero
    enddo
  enddo

  do iel = 1, ncel
!
    if ( noxyd .ge. 2 ) then
      rtp(iel,isca(if4m)) = zero
    endif
    if ( noxyd .ge. 3 ) then
      rtp(iel,isca(if5m)) = zero
    endif
!
    if ( ippmod(iccoal) .ge. 1 ) then
      rtp(iel,isca(if6m)) = zero
    endif
!
    rtp(iel,isca(if7m)) = zero
!
    if ( ihtco2 .eq. 1 ) then
      rtp(iel,isca(if8m)) = zero
    endif
!
    if ( ihth2o .eq. 1 ) then
      rtp(iel,isca(if9m)) = zero
    endif
!
    rtp(iel,isca(ifvp2m)) = zero
!
    if ( ieqco2.ge.1 ) then

      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas
      rtp(iel,isca(iyco2)) = oxyco2(ioxy)*wmco2/dmas
    endif
!
    if ( ieqnox .eq. 1 ) then
      rtp(iel,isca(iyhcn)) = zero
      rtp(iel,isca(iyno )) = zero
      rtp(iel,isca(ihox )) = h1init
    endif
!
  enddo

endif


!----
! FORMATS
!----

 9001 format(                                                     &
'                                                             ',/,&
'  user_coal_variables _init :                                ',/,&
'                   Variables Initialisation for              ',/,&
'                    Pulverised Coal by USer                  ',/,&
'                                                             ',/)

!----
! END
!----
end subroutine user_coal_iniv
