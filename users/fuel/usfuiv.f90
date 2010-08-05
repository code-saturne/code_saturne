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

subroutine usfuiv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   rdevel , rtuser , ra     )

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
!  PROPCE(IEL, IPPROC(IROM  (IPHAS))) is ROM(IEL,IPHAS) mean density kg/m3
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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           !  e ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfac(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor    ! ia ! <-- ! indirection for boundary faces ordering)       !
!  (nfabor, nphas) !    !     !                                                !
! itypfb           ! ia ! --> ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary developpement !
! ituser(nituse    ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfavor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
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
! rdevel(nrdeve)   ! ra ! <-> ! tab reel complementaire developemt             !
! rdevel(nideve)   ! ra ! <-- ! real work array for temporary developpement    !
! rtuser(nituse    ! ra ! <-- ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.f90"
include "pointe.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "entsor.f90"
include "parall.f90"
include "period.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "fuincl.f90"
include "ppincl.f90"
include "ppcpfu.f90"

!===============================================================================

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml), maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! LOCAL VARIABLES

integer          idebia, idebra
integer          iel, ige, mode, iphas, icla

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision xkent, xeent, d2s3

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

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION OF TRANSPORTED VARIABLES
!      RONLY IF THE COMPUTATION IS NOT A CONTINUATION
!===============================================================================

if ( isuite.eq.0 ) then

  iphas = 1

! --> Initialisation of k and epsilon (exemple)

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur(iphas).eq.2) then

    do iel = 1, ncel
      rtp(iel,ik(iphas))  = xkent
      rtp(iel,iep(iphas)) = xeent
    enddo

  elseif (itytur(iphas).eq.3) then

    do iel = 1, ncel
      rtp(iel,ir11(iphas)) = d2s3*xkent
      rtp(iel,ir22(iphas)) = d2s3*xkent
      rtp(iel,ir33(iphas)) = d2s3*xkent
      rtp(iel,ir12(iphas)) = 0.d0
      rtp(iel,ir13(iphas)) = 0.d0
      rtp(iel,ir23(iphas)) = 0.d0
      rtp(iel,iep(iphas))  = xeent
    enddo

  elseif (iturb(iphas).eq.50) then

    do iel = 1, ncel
      rtp(iel,ik(iphas))   = xkent
      rtp(iel,iep(iphas))  = xeent
      rtp(iel,iphi(iphas)) = d2s3
      rtp(iel,ifb(iphas))  = 0.d0
    enddo

  elseif (iturb(iphas).eq.60) then

    do iel = 1, ncel
      rtp(iel,ik(iphas))   = xkent
      rtp(iel,iomg(iphas)) = xeent/cmu/xkent
    enddo

  endif

! --> All the computation domain is initialised with air at TINITK
!                   ================================================

! ---- Computation of H1INIT and  H2INIT

  t1init = t0(iphas)
  t2init = t0(iphas)

! ------ Transported variables for droplets

  h2init = h02fol +  cp2fol*(t2init-trefth)

  do icla = 1, nclafu
    do iel = 1, ncel
      rtp(iel,isca(iyfol(icla))) = zero
      rtp(iel,isca(ing(icla) ))  = zero
      rtp(iel,isca(ihlf(icla)))  = h2init
    enddo
  enddo

! ------ Transported variables for the mix (droplets and carrying gases)

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo
  coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
  coefe(in2) = 1.d0 - coefe(io2)
  mode = -1
  call futhp1                                                     &
  !==========
 ( mode   , h1init , coefe  , t1init )

  do iel = 1, ncel
    rtp(iel,isca(ihm)) = h1init
  enddo

! ------ Transported variables for gaseous mixture
!        (passive scalars, variance, reactive species)

  do iel = 1, ncel
    rtp(iel,isca(ifvap )) = zero
    rtp(iel,isca(ifhtf )) = zero
    rtp(iel,isca(if4p2m)) = zero
    if ( ieqco2 .ge. 1 ) then
      rtp(iel,isca(iyco2)) = zero
    endif
    if ( ieqnox .eq. 1 ) then
      rtp(iel,isca(iyhcn)) = zero
      rtp(iel,isca(iyno )) = zero
      rtp(iel,isca(itaire)) = 20.d0+tkelvi
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
