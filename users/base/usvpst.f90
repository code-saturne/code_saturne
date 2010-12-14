!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
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
! lstcel(ncelps)   ! ia ! <-- ! list of cells in post-processing mesh          !
! lstfac(nfacps)   ! ia ! <-- ! list of interior faces in post-processing mesh !
! lstfbr(nfbrps)   ! ia ! <-- ! list of boundary faces in post-processing mesh !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
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

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use pointe
use entsor
use optcal
use numvar
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0, idbra0
integer          ipart
integer          nvar,   nscal , nphas , nvlsta
integer          ncelps, nfacps, nfbrps
integer          nideve, nrdeve, nituse, nrtuse

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

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
integer          iel, ifac, iloc, iphas, ivar, iclt
integer          idimt, ii , jj
integer          ientla, ivarpr
integer          imom1, imom2, ipcmo1, ipcmo2, idtcm
double precision pond
double precision rvoid(1)

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

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

if (ipart .eq. -1) then

  ! 1.1.1 Output of k=1/2(R11+R22+R33) for the Rij-epsilon model
  !       ------------------------------------------------------

  iphas = 1   ! We only consider phase 1 in this example

  if (itytur(iphas) .eq. 3) then

    ! Initialize variable name
    do ii = 1, 32
      namevr(ii:ii) = ' '
    enddo

    ! Variable name
    namevr = 'Turb energy'

    ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
    idimt = 1

    do iloc = 1, ncelps
      iel = lstcel(iloc)
      tracel(iloc) = 0.5d0*(  rtp(iel,ir11(iphas)) &
                            + rtp(iel,ir22(iphas)) &
                            + rtp(iel,ir33(iphas)) )
    enddo

    ! Values are not interlaced (dimension 1 here, so no effect).
    ientla = 0

    ! Values are defined on the work array, not on the parent.
    ivarpr = 0

    ! Output values; as we have no face values, we can pass a
    ! trivial array rvoid instead of trafac and trafbr.
    call psteva(ipart, namevr, idimt, ientla, ivarpr,  &
    !==========
                ntcabs, ttcabs, tracel, rvoid, rvoid)

  endif


  ! 1.1.2 Output of a combination of moments
  !       ----------------------------------
  ! We assume in this example that we have 2 temporal means (moments):
  !   <u>  for imom=1
  !   <uu> for imom=2
  ! We seek to plot <u'u'>=<uu>-<U>**2

  iphas = 1   ! We only consider phase 1 in this example

  if (nbmomt .ge. 2) then

    ! Initialize variable name
    do ii = 1, 32
      namevr (ii:ii) = ' '
    enddo

    ! Moment numbers:
    imom1 = 1
    imom2 = 2

    ! Position in 'propce' of the array of temporal accumulation for moments,
    ! propce(iel,ipcmom)
    ipcmo1 = ipproc(icmome(imom1))
    ipcmo2 = ipproc(icmome(imom2))

    ! Variable name
    namevr = '<upup>'

    ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
    idimt = 1

    ! The temporal accumulation for moments must be divided by the accumulated
    ! time, which id an array of size ncel or a single real number:
    ! - array of size ncel if idtmom(imom) > 0 : propce(iel, idtcm)
    ! - or simple real     if idtmom(imom) < 0 : dtcmom(idtcm)

    ! To improve this example's readability, we assume moments imom1 and imom2
    ! have been computed on the same time window.

    if(idtmom(imom1).gt.0) then
      idtcm = ipproc(icdtmo(idtmom(imom1)))
      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(iloc) =    propce(iel,ipcmo2)/max(propce(iel,idtcm),epzero)      &
                       - (propce(iel,ipcmo1)/max(propce(iel,idtcm),epzero))**2
      enddo
    elseif(idtmom(imom1).lt.0) then
      idtcm = -idtmom(imom1)
      do iloc = 1, ncelps
        iel = lstcel(iloc)
        tracel(iloc) =    propce(iel,ipcmo2)/max(dtcmom(idtcm),epzero)      &
                       - (propce(iel,ipcmo1)/max(dtcmom(idtcm),epzero))**2
      enddo
    endif

    ! Values are not interlaced (dimension 1 here, so no effect).
    ientla = 0

    ! Values are defined on the work array, not on the parent.
    ivarpr = 0

    ! Output values; as we have no face values, we can pass a
    ! trivial array rvoid instead of trafac and trafbr.
    call psteva(ipart, namevr, idimt, ientla, ivarpr,   &
    !==========
                ntcabs, ttcabs, tracel, rvoid, rvoid)

  endif

!===============================================================================
! 1.2. Examples of volume variables on the boundary mesh (ipart = -2)
!===============================================================================

else if  (ipart .eq. -2) then

  ! 1.2.1 Output of the density at the boundary
  !       -------------------------------------

  iphas = 1   ! We only consider phase 1 in this example

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Density at boundary'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 1

  ! Values are not interlaced (dimension 1 here, so no effect).
  ientla = 0

  ! We directly use the propfb array defined on the parent mesh.
  ivarpr = 1

  ! Output values; as we have only boundary face values, we can pass a
  ! trivial array rvoid instead of tracel and trafac.
  call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rvoid, rvoid,            &
              propfb(1,ipprob(irom(iphas))))


  ! 1.2.2 Output of the domain number in parallel
  !       ---------------------------------------

  ! This variable is independent of time, so we output it once
  ! only (see 'ntindp' below)

  if (ipass.eq.0 .and. irangp.ge.0) then

    ipass = ipass + 1

    ! Initialize variable name
    do ii = 1, 32
      namevr (ii:ii) = ' '
    enddo

    ! Variable name
    namevr = 'domain number'

    ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
    idimt = 1

    do iloc = 1, nfbrps
      ! ifac = lstfbr(iloc)
      trafbr(iloc) = irangp + 1
    enddo

    ! Values are not interlaced (dimension 1 here, so no effect).
    ientla = 0

    ! Values are defined on the work array, not on the parent.
    ivarpr = 0

    ! This variable is time-invariant;
    ! We assign a negative time to it and output it once only to avoid
    ! duplicating it at each output time.
    ntindp = -1

    ! Output values; as we have only boundary face values, we can pass a
    ! trivial array rvoid instead of trafac and trafbr.
    call psteva(ipart, namevr, idimt, ientla, ivarpr,  &
    !==========
                ntindp, ttcabs, rvoid, rvoid, trafbr)

  endif


!===============================================================================
! 1.3. Examples of volume variables on user meshes 1 or 2
!===============================================================================

else if  (ipart.eq.1 .or. ipart.eq.2) then

  ! 1.3.1 Output of the velocity
  !       ----------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Interpol velocity'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 3

  ! Values are interlaced.
  ientla = 1

  iphas = 1   ! We only consider phase 1 in this example

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(rtp(1,iu(iphas)), rtp(1,iv(iphas)), rtp(1,iw(iphas)))
    !==========
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(1 + (iloc-1)*idimt)  =            pond  * rtp(ii,iu(iphas))   &
                                   + (1.d0 - pond) * rtp(jj,iu(iphas))
    trafac(2 + (iloc-1)*idimt)  =            pond  * rtp(ii,iv(iphas))   &
                                   + (1.d0 - pond) * rtp(jj,iv(iphas))
    trafac(3 + (iloc-1)*idimt)  =            pond  * rtp(ii,iw(iphas))   &
                                   + (1.d0 - pond) * rtp(jj,iw(iphas))

  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    trafbr(1 + (iloc-1)*idimt) = rtp(ii, iu(iphas))
    trafbr(2 + (iloc-1)*idimt) = rtp(ii, iv(iphas))
    trafbr(3 + (iloc-1)*idimt) = rtp(ii, iw(iphas))

  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)


  ! 1.3.2 Output of the pressure
  !       ----------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Interpol pressure'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 1

  ! Values are not interlaced (dimension 1 here, so no effect).
  ientla = 0

  iphas = 1   ! We only consider phase 1 in this example

  ! Variable number
  ivar = ipr(iphas)

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ivar))
    !==========
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(iloc) =           pond  * rtp(ii, ivar)  &
                   + (1.d0 - pond) * rtp(jj, ivar)
  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    trafbr(iloc) = rtp(ii, ivar)
  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)


  ! 1.3.3 Output of the boundary temperature
  !       ----------------------------------

  iphas = 1   ! We only consider phase 1 in this example

  if (iscalt(iphas) .gt. 0) then

    ! Initialize variable name
    do ii = 1, 32
      namevr (ii:ii) = ' '
    enddo

    ! Variable name
    namevr = 'Boundary temperature'

    ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
    idimt = 1

    ! Values are not interlaced (dimension 1 here, so no effect).
    ientla = 0

    ! Set value to 0 for interior faces

    do iloc = 1, nfacps
      trafac(iloc) = 0.d0
    enddo

    ! Compute variable values on boundary faces.

    ivar = isca(iscalt(iphas))
    iclt = iclrtp(ivar,icoef)

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      ii = ifabor(ifac)
      trafbr(iloc) = coefa(ifac,iclt)+coefb(ifac,iclt)*rtp(ii, ivar)
    enddo

    ! Values are defined on the work array, not on the parent.
    ivarpr = 0

    ! Output values
    call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
    !==========
                ntcabs, ttcabs, tracel, trafac, trafbr)

  endif

  ! The examples below illustrate how to output a same variable in different
  ! ways (interlaced or not, using an indirection or not).


  ! 1.3.4 Output of the centers of gravity, interlaced
  !       --------------------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'face cog (interlaced)'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 3

  ! Values are interlaced
  ientla = 1

  do iloc = 1, nfacps

    ifac = lstfac(iloc)

    trafac(1 + (iloc-1)*idimt ) = cdgfac(1, ifac)
    trafac(2 + (iloc-1)*idimt ) = cdgfac(2, ifac)
    trafac(3 + (iloc-1)*idimt ) = cdgfac(3, ifac)
  enddo

  ! Compute variable values on boundary faces

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)

    trafbr(1 + (iloc-1)*idimt ) = cdgfbo(1, ifac)
    trafbr(2 + (iloc-1)*idimt ) = cdgfbo(2, ifac)
    trafbr(3 + (iloc-1)*idimt ) = cdgfbo(3, ifac)
  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)


  ! 1.3.5 Output of the centers of gravity, non-interlaced
  !       --------------------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'face cog (non-interlaced)'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 3

  ! Values are not interlaced
  ientla = 0

  do iloc = 1, nfacps

    ifac = lstfac(iloc)

    trafac(iloc)            = cdgfac(1, ifac)
    trafac(iloc + nfacps)   = cdgfac(2, ifac)
    trafac(iloc + 2*nfacps) = cdgfac(3, ifac)

  enddo

  ! Compute variable values on boundary faces

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)

    trafbr(iloc)            = cdgfbo(1, ifac)
    trafbr(iloc + nfbrps)   = cdgfbo(2, ifac)
    trafbr(iloc + 2*nfbrps) = cdgfbo(3, ifac)
  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)


  ! 1.3.6 Output of the centers of gravity, with indirection (parent-based)
  !       --------------------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'face cog (parent)'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 3

  ! Values are interlaced
  ientla = 1

  ! Values are defined on the parent.
  ivarpr = 1

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,   &
  !==========
              ntcabs, ttcabs, rvoid, cdgfac, cdgfbo)


!===============================================================================
! 1.4. Examples of volume variables on user meshes 3 or 4
!===============================================================================

else if  (ipart.ge.3 .and. ipart.le.4) then

  ! 1.4.1 Output of the velocity
  !       ----------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Velocity'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 3

  ! Values are not interlaced.
  ientla = 0

  iphas = 1   ! We only consider phase 1 in this example

  ! Variable number
  ivar = iu(iphas)

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(rtp(1,iu(iphas)), rtp(1,iv(iphas)), rtp(1,iw(iphas)))
    !==========
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(iloc) =                       pond  * rtp(ii, iu(iphas))   &
                               + (1.d0 - pond) * rtp(jj, iu(iphas))
    trafac(iloc + nfacps)    =           pond  * rtp(ii, iv(iphas))   &
                               + (1.d0 - pond) * rtp(jj, iv(iphas))
    trafac(iloc + 2*nfacps)  =           pond  * rtp(ii, iw(iphas))   &
                               + (1.d0 - pond) * rtp(jj, iw(iphas))
  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    trafbr(iloc )           = rtp(ii, iu(iphas))
    trafbr(iloc + nfbrps)   = rtp(ii, iv(iphas))
    trafbr(iloc + 2*nfbrps) = rtp(ii, iw(iphas))
  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,         &
  !==========
              ntcabs, ttcabs, rtp(1,ivar), trafac, trafbr)


  ! 1.5.2 Output of the pressure
  !       ----------------------

  ! Initialize variable name
  do ii = 1, 32
    namevr (ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Pressure'

  ! Variable dimension (1: scalar, 3: vector, 6/9: symm/non-symm tensor)
  idimt = 1

  ! Values are not interlaced (dimension 1 here, so no effect).
  ientla = 0

  iphas = 1   ! We only consider phase 1 in this example

  ! Variable number
  ivar = ipr(iphas)

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ivar))
    !==========
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(iloc)  =           pond  * rtp(ii, ivar)   &
                    + (1.d0 - pond) * rtp(jj, ivar)
  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    trafbr(iloc) = rtp(ii, ivar)
  enddo

  ! Values are defined on the work array, not on the parent.
  ivarpr = 0

  ! Output values
  call psteva(ipart, namevr, idimt, ientla, ivarpr,         &
  !==========
              ntcabs, ttcabs, rtp(1,ivar), trafac, trafbr)


endif ! end of test on post-processing mesh number

return

end subroutine
