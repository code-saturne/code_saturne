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

subroutine usetcl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb , izfppp ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose  :
! --------

!    User routine for extended physic
!    Electric module
!    Allocation of boundary conditions for varaibles unknown during usclim



!    This user subroutine is compulsory
!    =============================================


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Boundary condition types
! ========================

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in the 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.

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

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use elincl

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac  , ii     , iphas , iel
integer          i     , ntf    , nb    , id , itrouv
integer          izone
integer          nborne(nbtrmx)
integer          ilelt , nlelt

double precision rnbs2,capaeq
double precision sir(nelemx)   ,sii(nelemx)
double precision sirb(nbtrmx,6),siib(nbtrmx,6)
double precision ur(nbtrmx,6)  ,ui(nbtrmx,6)
double precision sirt(nbtrmx)  ,siit(nbtrmx)
character*200    chain

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  Allocation of Boundary Conditions
!       Loop on boundary faces
!       Allocation of family tyoe and properties
!       Allocation of boundary conditions
!===============================================================================

! 2.1 - Computation of intensity (A/m2) for each electrode
!       -----------------------------------------------------

!   Pre - initialisation

do i= 1,nbelec
  sir(i) = 0.d0
  sii(i) = 0.d0
enddo

do ntf= 1,nbtrf
  sirt(ntf) = 0.d0
  siit(ntf) = 0.d0
enddo

if(ntcabs.lt.(ntpabs+2)) then
  do ntf = 1,nbtrf
    uroff(ntf) = 0.d0
    uioff(ntf) = 0.d0
  enddo
endif

!     Loop on selected boundary faces

iphas = 1

do i=1,nbelec

  CHAIN = ' '
  write(chain,3000) ielecc(i)

  if ( ielect(i).ne. 0 ) then
    call getfbr(chain,nlelt,lstelt)
    !==========

    do ilelt = 1, nlelt

      ifac = lstelt(ilelt)

      iel = ifabor(ifac)

      do id=1,ndimve
        sir(i) = sir(i)                                           &
             +propce(iel,ipproc(idjr(id)))*surfbo(id,ifac)
      enddo

      if ( ippmod(ieljou) .eq. 4 ) then
        do id=1,ndimve
          sii(i) = sii(i)                                         &
               +propce(iel,ipproc(idji(id)))*surfbo(id,ifac)
        enddo
      endif

    enddo

  endif

enddo


! 2.2 - Definition of Voltage on each termin of transformers
!       ----------------------------------------------------------------

!  2.2.1 Computation of Intensity on each termin of transformers

do i=1,nbelec
  sirb(ielect(i),ielecb(i)) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siib(ielect(i),ielecb(i)) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirb(ielect(i),ielecb(i)) = sirb(ielect(i),ielecb(i))         &
                               +sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
       siib(ielect(i),ielecb(i)) = siib(ielect(i),ielecb(i))      &
                                  +sii(i)
    endif
  endif
enddo

!  2.2.2 RVoltage on each termin

do ntf=1,nbtrf

!      Primary and Secondary in Triangle

  if (ibrpr(ntf) .eq. 0 .and. ibrsec(ntf) .eq. 0 ) then

    nborne(ntf) = 3

    rnbs2 = 3.d0*rnbs(ntf)*rnbs(ntf)
    ur(ntf,1)=  1.154675d0*tenspr(ntf)/rnbs(ntf)                  &
      + (zr(ntf)*sirb(ntf,1)-zi(ntf)*siib(ntf,1))/rnbs2

    ur(ntf,2)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,2)-zi(ntf)*siib(ntf,2))/rnbs2
    ur(ntf,3)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,3)-zi(ntf)*siib(ntf,3))/rnbs2

    ui(ntf,1)=  0.d0                                              &
      + (zi(ntf)*sirb(ntf,1)+zr(ntf)*siib(ntf,1))/rnbs2
    ui(ntf,2)= -1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,2)+zr(ntf)*siib(ntf,2))/rnbs2
    ui(ntf,3)=  1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,3)+zr(ntf)*siib(ntf,3))/rnbs2

  else

    WRITE(NFECRA,*) 'Matrice sur le Transfo a ecrire'
    call csexit(1)

  endif
enddo

!  2.2.3 Total intensity for a transformer
!         (zero valued WHEN Offset established)

do ntf=1,nbtrf
  sirt(ntf) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siit(ntf) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirt(ielect(i)) = sirt(ielect(i)) + sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
      siit(ielect(i)) = siit(ielect(i)) + sii(i)
    endif
  endif
enddo

!  2.2.4 Take in account of Offset

capaeq = 3.d0

do ntf=1,nbtrf
  uroff(ntf) = uroff(ntf) + sirt(ntf)/capaeq
  if ( ippmod(ieljou) .eq. 4 ) then
    uioff(ntf) = uioff(ntf) + siit(ntf)/capaeq
  endif
enddo

! A reference transformer is assumed to have an Offset zero valued

if ( ntfref .gt. 0 ) then
  uroff(ntfref) = 0.d0
  uioff(ntfref) = 0.d0
endif

do ntf=1,nbtrf
  do nb=1,nborne(ntf)
    ur(ntf,nb) = ur(ntf,nb) + uroff(ntf)
    if ( ippmod(ieljou) .eq. 4 ) then
      ui(ntf,nb) = ui(ntf,nb) + uioff(ntf)
    endif
  enddo
enddo

! Print of UROFF (real part of offset potential)

write(nfecra,1500)
do ntf=1,nbtrf
  write(nfecra,2000) ntf,uroff(ntf)
enddo
write(nfecra,1501)

!  2.2.5 Take in account of Boundary Conditions


!     Loop on selected Boundary Faces

iphas = 1

do i=1,nbelec

  CHAIN = ' '
  write(chain,3000) ielecc(i)

  call getfbr(chain,nlelt,lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    iel = ifabor(ifac)

    itypfb(ifac,iphas) = iparoi

!     - Zone number
    izone = i

!      - Allocation of zone number
    izfppp(ifac) = izone

    if ( ielect(i) .ne. 0 ) then
      icodcl(ifac,isca(ipotr))   = 1
      rcodcl(ifac,isca(ipotr),1) = ur(ielect(i),ielecb(i))

      if ( ippmod(ieljou).eq.4  ) then
        icodcl(ifac,isca(ipoti))   = 1
        rcodcl(ifac,isca(ipoti),1) = ui(ielect(i),ielecb(i))
      endif

    else

      ii = ipotr
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0

      if ( ippmod(ieljou).eq. 4   ) then
        ii = ipoti
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif

  enddo

enddo

! 3 - Test, if not any reference transformer
!      a piece of wall may be at ground.

if ( ntfref .eq. 0 ) then

  itrouv = 0
  do ifac = 1, nfabor

    iphas = 1

    if ( itypfb(ifac,iphas) .eq. iparoi ) then

      if (icodcl(ifac,isca(ipotr)) .eq. 1 ) then

        if ( ippmod(ieljou).eq.3 ) then

          if ( abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20 ) then
            itrouv = 1
          endif

        else if ( ippmod(ieljou).eq.4 ) then

          if (icodcl(ifac,isca(ipoti)) .eq. 1 ) then

            if (abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20         &
           .and.abs(rcodcl(ifac,isca(ipoti),1)).lt.1.e-20 ) then
              itrouv = 1
            endif
          endif
        endif
      endif
    endif
  enddo

  if ( itrouv .eq. 0 ) then
    write(nfecra,1000)
    call csexit (1)
  endif

endif

!----
! FORMATS
!----

 1000 format(1X,' ERROR in JOULE : '                          ,/, &
       1X,' ====================   '                          ,/, &
      10X,' Lack of reference : choose a transformer for wich',/, &
      10X,' offset is assumed zero or a face at ground on the',/, &
          ' boundary')
 1500 format(/,2X,' ** INFORMATIONS ON TRANSFOMERS           ',/, &
         2X,'    ---------------------------------------'/    ,/, &
         1X,'      ---------------------------------'         ,/, &
         1X,'      Number of Transfo        UROFF    '        ,/, &
         1X,'      ---------------------------------')
 1501 format(1X,'      ---------------------------------')
 2000 format(10x,i6,12x,e12.5)
 3000 format(i7)

!----
! END
!----

end subroutine
