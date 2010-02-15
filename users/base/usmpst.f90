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

subroutine usmpst &
!================

 ( idbia0 , idbra0 , ipart  ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse , imodif ,                   &
   itypps , ifacel , ifabor , ifmfbr , ifmcel , iprfml ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis ,                                     &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )
!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR POUR LA MODIFICATION DES LISTES DE CELLULES
! OU FACES INTERNES ET DE BORD DEFINISSANT UN MAILLAGE DE POST
! TRAITEMENT EXISTANT ; CETTE ROUTINE EST APPELEE AUX PAS DE
! TEMPS AUQUEL CE MAILLAGE EST ACTIF, ET UNIQUEMENT POUR LES
! MAILLAGES POST UTILISATEUR PRINCIPAUX (NON ALIAS), SI TOUS LES
! "WRITERS" ASSOCIES A CE MAILLAGE OU SES ALIAS PERMETTENT
! CETTE MODIFICATION
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ipart            ! e  ! <-- ! numero du maillage post                        !
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
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! imodif           ! e  ! <-- ! 0 si maillage non modifie par cette            !
!                  !    !     ! fonction, 1 si modifie                         !
! itypps(3)        ! te ! <-- ! indicateur de presence (0 ou 1) de             !
!                  !    !     ! cellules (1), faces (2), ou faces de           !
!                  !    !     ! de bord (3) dans le maillage post              !
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
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! te ! --- ! macro tableau entier                           !
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
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet)         !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statis           ! tr ! <-- ! statistiques (lagrangien)                      !
!ncelet,nvlsta)    !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
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
include "pointe.h"
include "entsor.h"
include "optcal.h"
include "numvar.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ipart
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nvlsta
integer          ncelps , nfacps , nfbrps
integer          nideve , nrdeve , nituse , nrtuse, imodif

integer          itypps(3)
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
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

integer          ifac  , iphas
integer          ii   , jj
double precision vmin2, v2, w2


!===============================================================================

!     Remarque : le tableau ITYPPS permet de savoir si le maillage post
!                contient a l'origine des cellules, des faces internes,
!                ou des faces de bord (sur l'ensemble des processeurs).

!                Ceci permet d'avoir un traitement "generique" qui
!                peut fonctionner pour tous les numeros de maillage,
!                mais si le maillage post est vide a un instant de
!                post traitement donne, on ne saura plus s'il contenait
!                des cellules ou faces. Dans ce cas, il est preferable
!                d'utiliser explicitement le numero du maillage post
!                pour bien determiner s'il doit contenir des cellules
!                ou des faces.

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
!     1. TRAITEMENT DES MAILLAGES POST A REDEFINIR
!         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
!===============================================================================

!     Exemple :
!               pour les maillage post utilisateur, on ne conserve que
!               les mailles auxquelles la vitesse est superieure à
!               un seuil donne.


if (ipart.eq.3) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

!       SI LE MAILLAGE POST CONTIENT DES CELLULES
!       -----------------------------------------

  if (itypps(1) .eq. 1) then

    do ii = 1, ncel

      iphas = 1

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2
      if (v2 .ge. vmin2) then
        ncelps = ncelps + 1
        lstcel(ncelps) = ii
      endif

    enddo

!       SI LE MAILLAGE POST CONTIENT DES FACES INTERNES
!       -----------------------------------------------

  else if (itypps(2) .eq. 1) then

    do ifac = 1, nfac

      iphas = 1

      ii = ifacel(1, ifac)
      jj = ifacel(2, ifac)

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2
      w2 =   rtp(jj, iu(iphas))**2 + rtp(jj, iv(iphas))**2        &
           + rtp(jj, iw(iphas))**2

      if (v2 .ge. vmin2 .or. w2 .ge. vmin2) then
        nfacps = nfacps + 1
        lstfac(nfacps) = ifac
      endif

    enddo

!       SI LE MAILLAGE POST CONTIENT DES FACES DE BORD
!       ----------------------------------------------

  else if (itypps(3) .eq. 1) then

    do ifac = 1, nfabor

      iphas = 1

      ii = ifabor(ifac)

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2

      if (v2 .ge. vmin2) then
        nfbrps = nfbrps + 1
        lstfbr(nfbrps) = ifac
      endif

    enddo

  endif

!       Fin du test sur le type de mailles deja existantes

else if (ipart.eq.4) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

!       SELECTION DES FACES INTERNES
!       ----------------------------

  do ifac = 1, nfac

    iphas = 1

    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)

    v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2          &
         + rtp(ii, iw(iphas))**2
    w2 =   rtp(jj, iu(iphas))**2 + rtp(jj, iv(iphas))**2          &
         + rtp(jj, iw(iphas))**2

    if (     (v2 .ge. vmin2 .and. w2 .lt. vmin2)                  &
        .or. (v2 .lt. vmin2 .and. w2 .ge. vmin2)) then
      nfacps = nfacps + 1
      lstfac(nfacps) = ifac
    endif

  enddo

!       SELECTION DES FACES DE BORD
!       ---------------------------

  do ifac = 1, nfabor

    iphas = 1

    ii = ifabor(ifac)

    v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2          &
         + rtp(ii, iw(iphas))**2

    if (v2 .ge. vmin2) then
      nfbrps = nfbrps + 1
      lstfbr(nfbrps) = ifac
    endif

  enddo

endif
!     Fin du test sur le numero de maillage post.


return

end subroutine
