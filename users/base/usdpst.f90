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

subroutine usdpst &
!=================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

! Define additional post-processing writers and meshes.
!
! Post-processing writers allow outputs in different formats or with
! different format options and output frequancy than the default writer.
!
! Post-processing meshes are defined as a subset of the main meshe's
! cells or faces (interior and boundary).

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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
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
! lstcel(ncelet)   ! ia ! --- ! work array (list of cells)                     !
! lstfac(nfac)     ! ia ! --- ! work array (list of interior faces)            !
! lstfbr(nfabor)   ! ia ! --- ! work array (list of boundary faces)            !
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
use optcal
use entsor
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelet), lstfac(nfac), lstfbr(nfabor)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          indmod, icas, nbcas, ipart, nbpart, ipref, icat
integer          ntchrl, indgrp

integer          nlcel, nlfac , nlfbr, nlfam
integer          iel, ifac, ifam, ii
integer          idebia, idebra
integer          iflag1, iflag2, iel1, iel2
character*32     nomcas, nomfmt, nommai
character*96     nomrep, optfmt

double precision frchrl
double precision xfac  , yfac  , zfac

integer, allocatable, dimension(:) :: fam_list, fam_mask

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

nbcas  = 0
nbpart = 0

! "pointeurs" to the first free positions in 'ia' and 'ra'

idebia = idbia0
idebra = idbra0

!===============================================================================
! Create output writers for post-processing
! (one per case and per format, to be adapted by the user)
!===============================================================================

! Number of writers (case in the EnSight sense, study in the MED sense,
!                    or root of a CGNS tree)

nbcas = 4

do icas = 1, nbcas

  ! Miscellaneous initializations

  do ii = 1, len(nomcas)
    nomcas (II:II) = ' '
  enddo
  do ii = 1, len(nomrep)
    nomrep (ii:ii) = ' '
  enddo
  do ii = 1, len(nomfmt)
    nomfmt (ii:ii) = ' '
  enddo
  do ii = 1, len(optfmt)
    optfmt (ii:ii) = ' '
  enddo

  ! User definition:

  ! 'nomcas' and 'nomrep' respectively define the file names prefix and
  ! the corresponding directory path.
  ! If 'nomrep' is a local name of the "xxxx.ensight" or "xxxx.med" form,
  ! the script will automatically retreive the results to the 'RESU'
  ! directory, under a name such as XXXX.ENSIGHT.$DATE or XXXX.MED.$DATE.
  ! If 'nomrep' is of another form, it will have to be defined as a
  ! generic user output dire or directory so as to be copied.

  ! A user may also defined 'nomrep' as an absolute path, outside of the
  ! execution directory, in which case the results are output directly
  ! to that directory, and not managed by the script.

  ! 'nomfmt' allows choosing the output format ("EnSight Gold",
  ! "MED_fichier", or "CGNS").

  ! 'optfmt' allows the addition of a list of comma-separated
  ! format-specific output options:
  ! - EnSight:
  !      "text" ou "binary" (default),
  ! - EnSight, MED, or CGNS:
  !     "discard_polygons" to ignore polygons in output.
  !     "discard_polyhedra" to ignore polyhedra in output.
  ! - EnSight or MED :
  !     "divide_polygons" to divide polygons into triangles
  !     "divide_polyhedra" to divide polyhedra into tetrahedra and pyramids

  ! 'indmod' indicates if the meshes output using this writer will be:
  !     0: fixed,
  !     1: deformables with constant topology constante,
  !     2 : modifyable (may be redefined during the calculation through
  !         the 'usmpst' user subroutine).
  !     10: as indmod = 0, with a vertex displacement field
  !     11: as indmod = 1, with a vertex displacement field
  !     12: as indmod = 2, with a vertex displacement field

  ! 'ntchrl' defines the default output frequency (in time-steps)
  ! 'frchrl' defines the default output frequency (in seconds)
  ! Output at a specific time may still be forced or inhibited using the
  ! 'usnpst' user subroutine).

  frchrl = -1.d0

  if (icas .eq. 1) then

    nomcas = 'chr'
    nomrep = 'EnSight'
    nomfmt = 'EnSight Gold'
    optfmt = 'binary, discard_polygons'
    indmod = 0
    ntchrl = 4

  else if (icas .eq. 2) then

    nomcas = 'chr'
    nomrep = 'EnSight_text'
    nomfmt = 'ensight'
    optfmt = 'text, discard_polyhedra'
    indmod = 1
    ntchrl = ntchr

  else if (icas .eq. 3) then

    nomcas = 'modif'
    nomrep = 'EnSight'
    nomfmt = 'ensight'
    optfmt = 'discard_polyhedra'
    indmod = 2
    ntchrl = ntchr
    ntchrl = 2

  else if (icas .eq. 4) then

    nomcas = 'CHR'
    nomrep = ' '
    nomfmt = 'MED'
    optfmt = ' '
    indmod = 1
    ntchrl = ntchr

  endif

  ! Create writer

  call pstcwr (icas  , nomcas, nomrep, nomfmt, optfmt, indmod, ntchrl, frchrl)
  !==========

enddo

! Define number of additional postprocessing output meshes
!=========================================================

! 'nbpart' is the number of parts which will be generated (in the EnSight
! sense; the MED and CGNS equivalent terms are mesh and base respectively).

! A "part" may be any volume or surface defined through a selection of the
! main meshe's cells of faces.

! Example:
!
! 4 "parts", corresponding respectively to a mixed "interior faces"
! / "exterior faces" extraction, an extraction containing only
! interior faces, and 2 time-varying mesh pieces.

! We will later add a 5th "part", which is an alias of the second.

nbpart = 4

! Start of loop on user-defined parts
!====================================

do ipart = 1, nbpart

  ! Miscellaneous initializations
  !==============================

  indgrp = 1
  nlcel = 0
  nlfac = 0
  nlfbr = 0
  do iel = 1, ncelet
    lstcel(iel) = 0
  enddo
  do ifac = 1, nfac
    lstfac(ifac) = 0
  enddo
  do ifac = 1, nfabor
    lstfbr(ifac) = 0
  enddo

  do ii = 1, len(nommai)
    nommai(ii:ii) = ' '
  enddo

  ! Mark cells or faces included in the mesh (to be adapted by the user)
  !=====================================================================

  ! Note that this subroutine is called before boundary conditions
  ! are defined.

  ! Part 1:
  !   We select interior faces separating cells of color 2 from those
  !   of colors 3, (assuming no cell has both colors), as well as
  !   boundary faces of color 4.

  if (ipart .eq. 1) then

    nommai = 'Cut 1'

    ! Interior faces

    allocate(fam_list(nfml))
    allocate(fam_mask(nfml))

    fam_mask = 0

    ! Build mask on families matching colors 2 (1), 3 (2)

    call getfam('2', nlfam, fam_list)
    !==========

    do ifam = 1, nlfam
      fam_mask(fam_list(ifam)) = 1
    enddo

    call getfam('3', nlfam, fam_list)
    !==========

    do ifam = 1, nlfam
      fam_mask(fam_list(ifam)) = 2
    enddo

    deallocate(fam_list)

    ! Now that mask is built, test for adjacency

    do ifac = 1, nfac

      ! Adjacent cells
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)

      ! Adjacent cell flags
      iflag1 = fam_mask(ifmcel(iel1))
      iflag2 = fam_mask(ifmcel(iel2))

      ! Should the face belong to the extracted mesh ?
      if ((iflag1.eq.1.and.iflag2.eq.2).or.(iflag1.eq.2.and.iflag2.eq.1)) then
        nlfac = nlfac+1
        lstfac(nlfac)= ifac
      endif

    enddo

    deallocate(fam_mask)

    ! Boundary faces

    call getfbr('4', nlfbr, lstfbr)
    !==========

  ! Part 2:
  !   We select interior faces with y = 0.5

  else if (ipart .eq. 2) then

    nommai = 'Cut 2'

    ! Interior faces

    call getfac('plane[0, -1, 0, 0.5, epsilon = 0.0001]', nlfac, lstfac)
    !==========


  ! Part 3:
  !   We select all cells, and will modify the selection in 'usmpst'.

  else if (ipart .eq. 3) then

    nommai = 'Volume v > 0.5'

    nlcel = ncel
    nlfac = 0
    nlfbr = 0

  ! Part 4:
  !   We select all boundary faces, and will modify the selection in 'usmpst'.

  else if (ipart .eq. 4) then

    nommai = 'Surface "iso" v'

    nlcel = 0
    nlfac = 0
    nlfbr = nfabor

  endif

  ! Create post-processing mesh
  !============================

  call pstcma (ipart, nommai, indgrp, &
  !==========
               nlcel, nlfac, nlfbr, lstcel, lstfac, lstfbr)

  ! Associate extracted mesh and writer (to be adapted by the user)
  !================================================================

  if ((ipart .eq. 1) .or. (ipart .eq. 2)) then

    ! Associate post-processing meshes 1 and 2 with writers  1 and 2.
    icas = 1
    call pstass(ipart, icas)
    !==========
    icas = 2
    call pstass(ipart, icas)
    !==========

  else if ((ipart .eq. 3) .or. (ipart .eq. 4)) then

    ! Associate post-processing meshes 1 and 2 with writer 3.
    icas = 3
    call pstass(ipart, icas)
    !==========

  endif

  ! End of loop on user-defined parts
  !==================================

enddo

!===============================================================================
! Handle possible aliases; an alias is useful when we seek to assign an
! additional writer to an already defined mesh, with which we choose
! to output specific variables without duplicating the mesh in memory.
!===============================================================================

! Part 5: alias for part 2 (for independent varibles output)

ipart = 5
ipref = 2

call pstalm(ipart, ipref)
!==========

! Associate part 5 (alias of 2) with writer 3
icas = 3
call pstass(ipart, icas)
!==========

!===============================================================================
! Assign optional categories to user meshes
!===============================================================================

! Available categories are -1 (volume mesh) and -2 (boundary mesh).
! When a category is assigned to a user mesh, all standard
! (i.e. non-user) outputs which apply to the main category mesh
! also apply to the user mesh.

! In this example, variables output on the main volume mesh will
! also be output on the subset defined by post-processing mesh 3.
ipart = 3
icat = -1
call pstcat(ipart, icat)
!==========

return

!===============================================================================
! Formats
!===============================================================================

end subroutine
