!-------------------------------------------------------------------------------

!                      Code_Saturne version 2.0.0-rc1
!                      --------------------------

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

subroutine usvpst &
!================

 ( idbia0 , idbra0 , ipart  ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itypps , ifacel , ifabor , ifmfbr , ifmcel , iprfml ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis ,                                     &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Output additional variables on a postprocessing mesh.

! Several "automatic" postprocessing meshes may be defined:
! - The volume mesh (ipart=-1) if 'ichrvl' = 1
! - The boundary mesh (ipart=-2) if 'ichrbo' = 1
! - SYRTHES coupling surface (ipart < -2) if 'ichrsy' = 1
! - Cooling tower exchange zone meshes (ipart < -2) if 'ichrze' = 1
!
! Additional meshes (cells or faces) may also be defined using the
! 'usdpst' user subroutine, (and possibly modified using 'usmpst').

! This subroutine is called once for each post-processing mesh
! (with a different value of 'ipart') for each time step at which output
! on this mesh is active.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ipart            ! i  ! <-- ! number of the post-processing mesh (< 0 or > 0)!
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
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! ncelps           ! i  ! <-- ! number of cells in post-processing mesh        !
! nfacps           ! i  ! <-- ! number of interior faces in post-process. mesh !
! nfbrps           ! i  ! <-- ! number of boundary faces in post-process. mesh !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! itypps(3)        ! ia ! <-- ! global presence flag (0 or 1) for cells (1),   !
!                  !    !     ! interior faces (2), or boundary faces (3) in   !
!                  !    !     ! post-processing mesh                           !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! lstcel(ncelps)   ! ia ! <-- ! list of cells in post-processing mesh          !
! lstfac(nfacps)   ! ia ! <-- ! list of interior faces in post-processing mesh !
! lstfbr(nfbrps)   ! ia ! <-- ! list of boundary faces in post-processing mesh !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statis           ! ra ! <-- ! statistic values (Lagrangian)                  !
!  (ncelet, nvlsta)!    !     !                                                !
! tracel(*)        ! ra ! --- ! work array for post-processed cell values      !
! trafac(*)        ! ra ! --- ! work array for post-processed face values      !
! trafbr(*)        ! ra ! --- ! work array for post-processed boundary face v. !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "optcal.h"
include "numvar.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0
integer          ipart
integer          ndim,   ncelet, ncel,   nfac,   nfabor
integer          nfml,   nprfml
integer          nnod,   lndfac, lndfbr, ncelbr
integer          nvar,   nscal , nphas , nvlsta
integer          ncelps, nfacps, nfbrps
integer          nideve, nrdeve, nituse, nrtuse

integer          itypps(3)
integer          ifacel(2,nfac), ifabor(nfabor)
integer          ifmfbr(nfabor), ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision statis(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

character*32     namevr

integer          ntindp
integer          iel   , ifac  , iloc  , iphas, ivar , iclt
integer          idimt , ii   , jj
integer          idimte, itenso, ientla, ivarpr
integer          imom1, imom2, ipcmo1, ipcmo2, idtcm
double precision pond
double precision rbid(1)

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================


!===============================================================================
! 1. Handle variables to output
!    MUST BE FILLED IN by the user at indicated places
!===============================================================================

! A post-processing id a "part" (using the EnSight vocabulary; the MED and
! CGNS equivalents are "mesh" and "base" respectively).
! The user will have defined post-processing meshes in 'usdpst' ('nbpart'
! post-processing meshes).

! This subroutine is called once for each post-processing mesh
! (with a different value of 'ipart') for each time step at which output
! on this mesh is active. For each mesh and for all variables we wish to
! post-process here, we must define certain parameters and pass them to
! the 'psteva' subroutine, which is in charge of the actual output.
! These parameters are:

! namevr <-- variable name
! idimt  <-- variable dimension
!            (1: scalar, 3: vector, 6: symmetric tensor, 9: tensor)
! ientla <-- when idimt >1, this flag specifies if the array containing the
!            variable values is interlaced when ientla = 1
!            (x1, y1, z1, x2, y2, z2, x3, y3, z3...), or non-interlaced
!            when ientla = 0 (x1,x2,x3,...,y1,y2,y3,...,z1,z2,z3,...).
! ivarpr <-- specifies if the array containing the variable is defined on
!            the "parent" mesh or locally.
!            Even if the 'ipart' post-processing mesh contains all the
!            elements of its parent mesh, their numbering may be different,
!            especially when different element types are present.
!            The 'tracel' array passed as an argument to 'psteva' is built
!            relative to the numbering of the 'ipart' post-processing mesh.
!            To post-process a variable contained for example in the 'rtuser'
!            array, it should first be re-ordered, as shown here:
!              do iloc = 1, ncelps
!                iel = lstcel(iloc)
!                tracel(iloc) = rtuser(iel)
!              enddo
!            An alternative option is provided, to avoid unnecessary copies:
!            an array defined on the parent mesh, such our 'rtuser' example,
!            may be passed directly to 'psteva', specifying that values
!            are defined on the parent mesh instead of the post-processing mesh,
!            by setting the 'ivarpr' argument of 'psteva' to 1.

! Note: be cautious with variable name lengths.

! We allow up to 32 characters here, but names may be truncted depending on the
! output format:

! - 19 characters for EnSight
! - 32 characters for MED

! The nam length is not limited internally, so in case of 2 variables whoses
! names differ only after the 19th character, the corresponding names will
! both appear in the ".case" file; simply renaming one of the field descriptors
! in this text file will correct the output.

! Whitespace at the beginning or the end of a line is truncated automatically.
! Depending on the format used, prohibited characters (under EnSight, characters
! (  ) ] [ + - @           ! # * ^ $ / as well as white spaces and tabulations
! are automatically replaced by the _ character.

! Examples:

! For post-processing mesh 2, we output the velocity, pressure, and prescribed
! temperature at boundary faces (as well as 0 on possible interior faces)

! For post-processing mesh 1, we output all the variables usually
! post-processed, using a more compact coding.

! Examples given here correspond to the meshes defined in usdpst.f90 and
! modified in usmpst.f90.


!===============================================================================
! 1.1. Examples of volume variables on the main volume mesh (ipart = -1)
!===============================================================================

if (ipart.eq.1) then

  ! Initialisation
  ! no user intervention required --------------

  ! Initialize variable name
  do ii = 1, 32
    namevr(ii:ii) = ' '
  enddo

  ! Variable name
    namevr = 'interpolated Temperature'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 1

  ! Values are not interlaced (dimension 1 here, so no effect).
  ientla = 0

  ! Computation of variable values on internal faces
  ! To simplify the example, we use here a simple linear interpolation
  ! In muti-processor calculations, if neighouring cells are used, it is 
  ! necessary to exchange informations as usual. In calculations with periodicity
  ! and several processors, it is in addition necessary to call the routine 'parcom'
  ! before calling the routine 'percom'.

  ! user intervention:

  if (irangp.ge.0) then
    call  parcom (rtp(1,isca(1)))
    !===========
  endif

  if (iperio.eq.1) then
    ivar = isca(1)
    idimte = 0
    itenso = 0
    call percom &
    !==========
  ( idimte, itenso,   &
    rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),  &
    rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),  &
    rtp(1,ivar), rtp(1,ivar), rtp(1,ivar))
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(iloc) =         pond  * rtp(ii,isca(1))  &
                 + (1.d0 - pond) * rtp(jj,isca(1))
  enddo

  ! Writing of calculated values

  ivarpr = 0

  call psteva(ipart , namevr, idimt, ientla, ivarpr,  &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)

elseif (ipart .eq. 2) then

  ! Post-processing of velocity
  ! ---------------------------

  ! Initialisation
  ! no user intervention required --------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Temperature'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 1

  ! Values are not interlaced (dimension 1 here, so no effect).
  ientla = 0

  do iloc = 1, ncelps
    iel = lstcel(iloc)
    tracel(iloc)= rtp(iel,isca(1))
  enddo

  ivarpr = 0

  call psteva(ipart , namevr, idimt, ientla, ivarpr,   &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)

endif ! end of test on post-processing mesh number

return

end subroutine
