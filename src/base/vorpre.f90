!-------------------------------------------------------------------------------

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

subroutine vorpre &
!================

 ( idbia0 , idbra0 , ifinia , ifinra ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nvar   , nscal  , nphas  , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   irepvo ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   propce , propfa , propfb ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE DE PREPATATION DE LA METHODE DES VORTEX
!    Gestion memoire, connectivites, ...
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
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
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! irepvo           ! te ! <-- ! tab entier pour reperage des faces de          !
!                  !    !     ! bord pour la methode des vortex                !
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
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "vortex.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , ifinia , ifinra
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          irepvo(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, ii, iphas
integer          ient, ipcvis, ipcrom
integer          iappel
integer          isurf(nentmx)
double precision xx, yy, zz
double precision xxv, yyv, zzv

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

nvomax = 0
do ient = 1, nnent
  nvomax = max(nvort(ient),nvomax)
enddo

! NVOMAX = nombre max de vortex (utilise plus tard)

do ient = 1, nnent
  icvor2(ient) = 0
enddo

do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    icvor2(ient) = icvor2(ient) + 1
  endif
enddo

! ICVOR2 = compteur du nombre local de faces
!   utilisant des vortex a l'entree IENT

icvmax = 0
if(irangp.ge.0) then
  do ient = 1, nnent
    icvor(ient) = icvor2(ient)
    call parcpt(icvor(ient))
    !==========
    icvmax = max(icvmax,icvor(ient))
  enddo
else
  do ient = 1, nnent
    icvor(ient) = icvor2(ient)
    icvmax = max(icvmax,icvor(ient))
  enddo
endif

! ICVOR = nombre global de faces utilisant des vortex a l'entree IENT

! ICVMAX = max du nombre global de faces utilisant des vortex
! (toutes entrees confondues).

iappel = 2
call memvor                                                       &
!==========
 ( idebia , idebra , iappel , nfabor , ifinia , ifinra )

idebia = ifinia
idebra = ifinra

!===============================================================================
! 2. CONSTRUCTION DE LA " GEOMETRIE GOBALE "
!===============================================================================

do ient = 1, nnent
  icvor2(ient) = 0
  xsurfv(ient) = 0.d0
  isurf(ient)  = 0
enddo

! Chaque processeur stocke dans les tableaux RA(IW1X),...
! les coordonnees des faces ou il doit ensuite utiliser des vortex

iphas  = 1
ipcvis = ipproc(iviscl(iphas))
ipcrom = ipproc(irom(iphas))
do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    iel = ifabor(ifac)
    icvor2(ient) = icvor2(ient) + 1
    ra(iw1x+(ient-1)*icvmax+icvor2(ient)-1)= cdgfbo(1,ifac)
    ra(iw1y+(ient-1)*icvmax+icvor2(ient)-1)= cdgfbo(2,ifac)
    ra(iw1z+(ient-1)*icvmax+icvor2(ient)-1)= cdgfbo(3,ifac)
    ra(iw1v+(ient-1)*icvmax+icvor2(ient)-1) =                     &
      propce(iel,ipcvis)/propce(iel,ipcrom)
    xsurfv(ient) = xsurfv(ient) + sqrt(surfbo(1,ifac)**2          &
      + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
!         Vecteur surface d'une face de l'entree
    if (isurf(ient).eq.0) then
      surf(1,ient) = surfbo(1,ifac)
      surf(2,ient) = surfbo(2,ifac)
      surf(3,ient) = surfbo(3,ifac)
      isurf(ient)  = 1
    endif
  endif
enddo

if(irangp.ge.0) then
  do ient = 1, nnent
    call parsom(xsurfv(ient))
    !==========
  enddo
endif

! -------------
! En parallele
! -------------
if(irangp.ge.0) then
  do ient = 1, nnent
    call paragv                                                   &
    !==========
 ( icvor2(ient), icvor(ient),                                     &
   ra(iw1x  + (ient-1)*icvmax)   ,                                &
   ra(ixyzv + (ient-1)*3*icvmax) )
    call paragv                                                   &
    !==========
 ( icvor2(ient), icvor(ient),                                     &
   ra(iw1y  + (ient-1)*icvmax)              ,                     &
   ra(ixyzv + (ient-1)*3*icvmax +   icvmax) )
    call paragv                                                   &
    !==========
 ( icvor2(ient), icvor(ient),                                     &
   ra(iw1z  + (ient-1)*icvmax)              ,                     &
   ra(ixyzv + (ient-1)*3*icvmax + 2*icvmax) )
    call paragv                                                   &
    !==========
 ( icvor2(ient), icvor(ient),                                     &
   ra(iw1v  + (ient-1)*icvmax)   ,                                &
   ra(ivisv + (ient-1)*3*icvmax) )
  enddo

!  -> A la fin de cette etape, tous les processeurs connaissent
!     les coordonees des faces d'entree

else
! ----------------------
! Sur 1 seul processeur
! ----------------------
  do ient = 1,nnent
    do ii = 1, icvor(ient)
      ra(ixyzv+(ient-1)*3*icvmax+ii-1)=                           &
           ra(iw1x+(ient-1)*icvmax + ii -1)
      ra(ixyzv+(ient-1)*3*icvmax+icvmax + ii-1) =                 &
           ra(iw1y+(ient-1)*icvmax + ii -1)
      ra(ixyzv+(ient-1)*3*icvmax+2*icvmax + ii-1) =               &
           ra(iw1z+(ient-1)*icvmax + ii -1)
      ra(ivisv+(ient-1)*icvmax+ii-1) =                            &
           ra(iw1v+(ient-1)*icvmax+ii-1)
    enddo
  enddo
endif

!===============================================================================
! 3. CONSTRUCTION DE LA CONNECTIVITE
!===============================================================================

do ient = 1, nnent
  icvor2(ient) = 0
  do ifac = 1, icvmax
    ia(iifagl+(ient-1)*icvmax+ifac-1) = 0
  enddo
enddo

! On cherche ensuite le numero de la ligne du tableau RA(IXYZV) qui est
! associe a la Ieme face d'entree utilisant des vortex (dans la
! numerotation chronologique que suit ICVOR2).

do ifac = 1, nfabor
  ient = irepvo(ifac)
  if(ient.ne.0) then
    icvor2(ient) = icvor2(ient) + 1
    do ii = 1, icvor(ient)
      xx = cdgfbo(1,ifac)
      yy = cdgfbo(2,ifac)
      zz = cdgfbo(3,ifac)
      xxv = ra(ixyzv+(ient-1)*3*icvmax+ii-1)
      yyv = ra(ixyzv+(ient-1)*3*icvmax+icvmax+ii-1)
      zzv = ra(ixyzv+(ient-1)*3*icvmax+2*icvmax+ii-1)
      if(abs(xxv-xx).lt.epzero.and.abs(yyv-yy).lt.epzero.and.     &
           abs(zzv-zz).lt.epzero) then
        ia(iifagl+(ient-1)*icvmax+icvor2(ient)-1) = ii
      endif
    enddo
  endif
enddo

! La methode de vortex va generer un tableau de vitesse RA(IUVOR)
! qui aura la meme structure que RA(IXYZV).
! Le tableau RA(IXYZV) sera envoyee a tous les processeurs
! la vitesse a imposer à la Ieme face se trouvera à la ligne IA(IIFAGL+I)

iappel = 3
call memvor                                                       &
!==========
 ( idebia , idebra , iappel , nfabor , ifinia , ifinra )

! ---
! FIN
! ---

return
end subroutine
