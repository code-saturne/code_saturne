!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine memnav &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviscf , iviscb , ivisfi , ivisbi ,                            &
   idam   , ixam   ,                                              &
   idrtp  , igrdp  , ismbr  , irovsd ,                            &
   iw1    , iw2    , iw3    , iw4    , iw5    , iw6    , iw7    , &
   iw8    , iw9    , iw10   , idfrcx , ifrchy , idfrhy ,          &
   icoefu , iesflm , iesflb , itrava , iximpa ,                   &
   iuvwk  ,                                                       &
   ifinia , ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE NAVIER-STOKES

!-------------------------------------------------------------------------------
! Arguments
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
! iviscf, b        ! e  ! --> ! "pointeur" sur viscf, viscb                    !
! ivisfi, bi       ! e  ! --> ! "pointeur" sur viscfi, viscbi                  !
! idam, ixam       ! e  ! --> ! "pointeur" sur dam, xam                        !
! idrtp            ! e  ! --> ! "pointeur" sur drtp                            !
! igrdp            ! e  ! --> ! "pointeur" sur grdp                            !
! ismbr            ! e  ! --> ! "pointeur" sur smbr                            !
! irovsd           ! e  ! --> ! "pointeur" sur rovsdt                          !
! iw1,2,...,10     ! e  ! --> ! "pointeur" sur w1 a w9                         !
! idfrcx           ! e  ! --> ! "pointeur" sur dfrcxt                          !
! ifrchy           ! e  ! --> ! "pointeur" sur frchy                           !
! idfrhy           ! e  ! --> ! "pointeur" sur dfrchy                          !
!iesflm, iesflb    ! e  ! --> ! "pointeur" sur esflum et esflub                !
! itrava           ! e  ! --> ! "pointeur" sur trava                           !
! iximpa           ! e  ! --> ! "pointeur" sur ixmpa                           !
! iuvwk            ! e  ! --> ! "pointeur" sur uvwk                            !
! ifinia           ! e  ! --> ! pointeur de la premiere case libre             !
!                  !    !     !  dans ia en sortie                             !
! ifinra           ! e  ! --> ! pointeur de la premiere case libre             !
!                  !    !     !  dans ia en sortie                             !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          iviscf , iviscb , ivisfi , ivisbi
integer          idam   , ixam
integer          idrtp  , igrdp  , ismbr  , irovsd
integer          iw1    , iw2    , iw3    , iw4    , iw5    , iw6
integer          iw7    , iw8    , iw9    , iw10
integer          idfrcx , ifrchy , idfrhy
integer          icoefu , iesflm , iesflb
integer          itrava , iximpa , iuvwk
integer          ifinia , ifinra

integer          idebia , idebra , iphas, irij

integer          iescat
integer          iirnpn

!===============================================================================

!---> INITIALISATION

idebia = idbia0
idebra = idbra0


!---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA

irij = 0
do iphas = 1, nphas
  if(itytur(iphas).eq.3.and.irijnu(iphas).eq.1) then
    irij = 1
  endif
enddo

iescat = 0
do iphas = 1, nphas
  if(iescal(iestot,iphas).gt.0) then
    iescat = 1
  endif
enddo

! Attention :
!   ci-dessous, on fait pointer VISCFI et VISCBI
!   sur VISCF et VISCB respectivement
!   dans les cas ou ils sont identiques.

!     Pour la norme de resolp
iirnpn = 0
if(irnpnw.eq.1) iirnpn = 1

ivisfi =       idebra
iviscf =       ivisfi + nfac  *irij
ivisbi =       iviscf + nfac
iviscb =       ivisbi + nfabor*irij
idam   =       iviscb + nfabor
ixam   =       idam   + ncelet
idrtp  =       ixam   + nfac*2
igrdp  =       idrtp  + ncelet
ismbr  =       igrdp  + ncelet*3
irovsd =       ismbr  + ncelet
iw1    =       irovsd + ncelet
iw2    =       iw1    + ncelet
iw3    =       iw2    + ncelet
iw4    =       iw3    + ncelet
iw5    =       iw4    + ncelet
iw6    =       iw5    + ncelet
iw7    =       iw6    + ncelet
iw8    =       iw7    + ncelet
iw9    =       iw8    + ncelet
iw10   =       iw9    + ncelet
idfrcx =       iw10   + ncelet*iirnpn
ifrchy =       idfrcx + ncelet*3*nphas*iphydr
idfrhy =       ifrchy + ncelet*ndim*icalhy
icoefu =       idfrhy + ncelet*ndim*icalhy
iesflm =       icoefu + nfabor*ndim
iesflb =       iesflm + nfac*iescat
iximpa =       iesflb + nfabor*iescat
if(nterup.gt.1) then
  iuvwk  =       iximpa + ncelet*ndim*nphas
  itrava =       iuvwk  + ncelet*ndim*nphas
  ifinra =       itrava + ncelet*ndim*nphas
else
  iuvwk  =       iximpa
  itrava =       iuvwk
  ifinra =       itrava
endif

! Dans une phase d'optimisation memoire, on pourra faire pointer
!  esflum et esflub sur fluint et flubrd dans le cas (rare)
!  avec estimateur total et reconstruction par moindres carres

!---> VERIFICATION

CALL IASIZE('MEMNAV',IFINIA)
!==========

CALL RASIZE('MEMNAV',IFINRA)
!==========

return
end
