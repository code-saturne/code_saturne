!-------------------------------------------------------------------------------

!VERS

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usd3pi &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   ra     )

!===============================================================================
! PURPOSE :
! --------
! Initialisation of the transported variables of the 3-point
! chemistry combustion model. The latter ones are used to calculate
! diffusion flame properties by applying detailed thermochemistry data.
!
! FEATURES:
! --------
! - This subroutine is called at the beginning of both a restart-calculation
!   and a new calculation.
! - The user is able to modify the initial conditions of the tranported
!   variables.
! - All variables of state deducible from the transported variables, as well as
!   the detailed thermochemistry data (JANAF-table) are located in the array named
!
!   a) PROPCE (Properties in the cell center).
!   b) PROPFA (Properties at internal faces).
!   c) PROPFB (Properties at boundary faces).
!
!   Examples:
!   PROPCE(IEL,IPPROC(IROM  )) =  IROM of cell IEL
!   PROPCE(IEL,IPPROC(ICP   )) =  ICP  of cell IEL
!
!   PROPFA(IFAC,IPPROF(IFLUMA(IVAR )))=  FLUMAS of IVAR at the internal face IFAC
!
!   PROPFB(IFAC,IPPROB(IROM  ))=  ROMB at the boundary face IFAC
!   PROPFB(IFAC,IPPROB(IFLUMA(IVAR )))=  FLUMAB of IVAR at the boundary face IFAC

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! maxelt           !  e ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor    ! ia ! <-- ! indirection for boundary faces ordering)       !
! itypfb           ! ia ! --> ! boundary face types                            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!                  !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! coefu            ! ra ! --- ! tab de trav                                    !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use mesh

!===============================================================================

implicit none

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          maxelt, lstelt(maxelt)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ra(*)


! LOCAL VARIABLES

integer          idebia, idebra
integer          iel, igg
double precision coefg(ngazgm)


!===============================================================================



!===============================================================================
! 0. CONTROL OUTPUT
!===============================================================================

write(nfecra,9001)

!===============================================================================
! 1.  INITIALISATION OF LOCAL VARIABLES
!===============================================================================

idebia = idbia0
idebra = idbra0

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! 2. INITIALISATION OF TRANSPORTED VARIABLES
!===============================================================================

if ( isuite.eq.0 ) then

  do iel = 1, ncel

! ----- Mean Mixture Fraction
    rtp(iel,isca(ifm))   = fs(1)

! ----- Variance of Mixture Fraction
    rtp(iel,isca(ifp2m)) = zero

! ----- Enthalpy
    if ( ippmod(icod3p).eq.1 ) then
      rtp(iel,isca(ihm)) = hinfue*fs(1)+hinoxy*(1.d0-fs(1))
    endif

  enddo

endif


!----
! FORMATS
!----

 9001 format(                                                           &
'                                                             ',/,&
'  usd3pi : user defined initialisation of variables          ',/,&
'                                                             ',/)


!----
! END
!----


return
end subroutine
