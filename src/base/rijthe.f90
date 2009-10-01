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

subroutine rijthe &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , grarox , graroy , graroz , smbr   ,          &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! TERMES DE GRAVITE
!   POUR Rij et EPSILON
! VAR  = R11 R22 R33 R12 R13 R23 EP

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de phase                                !
! ivar             ! e  ! <-- ! numero de variable                             !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! grarox,y,z       ! tr ! <-- ! tableau de travail pour grad rom               !
!  (ncelet)        !    !     !                                                !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , ivar   , isou   , ipp

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision grarox(ncelet), graroy(ncelet), graroz(ncelet)
double precision smbr(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iel
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph

double precision uns3, const, kseps
double precision prdtur, r1t, r2t, r3t
double precision g11, g22, g33, g12, g13, g23, gkks3
double precision g11p, g22p, g33p
double precision phit11, phit22, phit33, phit12, phit13, phit23
double precision aa, bb


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ir11ip = ir11(iphas)
ir22ip = ir22(iphas)
ir33ip = ir33(iphas)
ir12ip = ir12(iphas)
ir13ip = ir13(iphas)
ir23ip = ir23(iphas)
ieiph  = iep (iphas)

if(iscalt(iphas).gt.0.and.nscal.ge.iscalt(iphas)) then
  prdtur = sigmas(iscalt(iphas))
else
  prdtur = 1.d0
endif

const = -1.5d0*cmu/prdtur
uns3  = 1.d0/3.d0

!===============================================================================
! 2. TERMES POUR RIJ :
!      ROM*VOLUME*dRij/dt =
!                     ... + (Gij - CRIJ3*(Gij-Delta ij Gkk/3))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!===============================================================================


if     (ivar.eq.ir11ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit11 = -crij3*(g11-gkks3)

    smbr(iel) = smbr(iel) + (g11+phit11)*volume(iel)

  enddo

elseif (ivar.eq.ir22ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit22 = -crij3*(g22-gkks3)

    smbr(iel) = smbr(iel) + (g22+phit22)*volume(iel)

  enddo

elseif (ivar.eq.ir33ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit33 = -crij3*(g33-gkks3)

    smbr(iel) = smbr(iel) + (g33+phit33)*volume(iel)

  enddo

elseif (ivar.eq.ir12ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g12 = const*kseps*     (r1t*gy+r2t*gx)

    phit12 = -crij3* g12

    smbr(iel) = smbr(iel) + (g12+phit12)*volume(iel)

  enddo

elseif (ivar.eq.ir13ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g13 = const*kseps*     (r1t*gz+r3t*gx)

    phit13 = -crij3* g13

    smbr(iel) = smbr(iel) + (g13+phit13)*volume(iel)

  enddo

elseif (ivar.eq.ir23ip) then

  do iel = 1, ncel

    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g23 = const*kseps*     (r2t*gz+r3t*gy)

    phit23 = -crij3* g23

    smbr(iel) = smbr(iel) + (g23+phit23)*volume(iel)

  enddo

!===============================================================================
! 3. TERMES POUR EPSILON :
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*(EPS/K)*MAX(0,(Gkk/2))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!            On simplifie (EPS/K) en notant
!                GijP = -(1.5 CMU/PRDTUR)         (Rit Gj + Rjt Gi)
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*        MAX(0,(GkkP/2))*VOLUME
!===============================================================================


elseif (ivar.eq.ieiph ) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    g11p = const*      2.d0*(r1t*gx       )
    g22p = const*      2.d0*(r2t*gy       )
    g33p = const*      2.d0*(r3t*gz       )

    aa = 0.d0
    bb = 0.5d0*(g11p+g22p+g33p)
    smbr(iel) = smbr(iel) + ce1*max(aa,bb)*volume(iel)

  enddo

endif

return

end subroutine
