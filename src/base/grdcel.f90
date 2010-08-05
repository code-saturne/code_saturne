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

subroutine grdcel &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp ,                                     &
   dpdx   , dpdy   , dpdz   ,                                     &
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! APPEL DES DIFFERENTES ROUTINES DE CALCUL DE GRADIENT CELLULE

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
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ivar             ! e  ! <-- ! numero de la variable                          !
!                  !    !     !   destine a etre utilise pour la               !
!                  !    !     !   periodicite uniquement (pering)              !
!                  !    !     !   on pourra donner ivar=0 si la                !
!                  !    !     !   variable n'est ni une composante de          !
!                  !    !     !   la vitesse, ni une composante du             !
!                  !    !     !   tenseur des contraintes rij                  !
! imrgra           ! e  ! <-- ! methode de reconstruction du gradient          !
!                  !    !     !  0 reconstruction 97                           !
!                  !    !     !  1 moindres carres                             !
!                  !    !     !  2 moindres carres support etendu              !
!                  !    !     !    complet                                     !
!                  !    !     !  3 moindres carres avec selection du           !
!                  !    !     !    support etendu                              !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
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
! pvar  (ncelet    ! tr ! <-- ! variable (pression)                            !
! coefap,coefbp    ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! fextx,y,z        ! tr ! <-- ! force exterieure generant la pression          !
!   (ncelet)       !    !     !  hydrostatique                                 !
! dpdx,dpdy        ! tr ! --> ! gradient de pvar                               !
! dpdz (ncelet     !    !     !                                                !
! dpdxa (ncelet    ! tr ! --- ! tableau de travail pour le grad de p           !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "pointe.f90"
include "period.f90"
include "parall.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          ivar   , imrgra , inc    , iccocg , nswrgp
integer          imligp ,iwarnp  , iphydp , nfecra
double precision epsrgp , climgp , extrap

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision dpdxa(ncelet),dpdya(ncelet),dpdza(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          idimte , itenso
integer          iiu(nphsmx),iiv(nphsmx),iiw(nphsmx)
integer          iitytu(nphsmx)
integer          iir11(nphsmx),iir22(nphsmx),iir33(nphsmx)
integer          iir12(nphsmx),iir13(nphsmx),iir23(nphsmx)
integer          imlini

double precision climin

!===============================================================================


idebia = idbia0
idebra = idbra0

!===============================================================================
! 0. PREPARATION POUR PERIODICITE DE ROTATION
!===============================================================================

! Par defaut, on traitera le gradient comme un vecteur ...
!   (i.e. on suppose que c'est le gradient d'une grandeurs scalaire)

! S'il n'y a pas de rotation, les echanges d'informations seront
!   faits par percom (implicite)

! S'il y a une ou des periodicites de rotation,
!   on determine si la variables est un vecteur (vitesse)
!   ou un tenseur (de Reynolds)
!   pour lui appliquer dans percom le traitement adequat.
!   On positionne IDIMTE et ITENSO
!   et on recupere le gradient qui convient.
! Notons que si on n'a pas, auparavant, calcule et stocke les gradients
!   du halo on ne peut pas les recuperer ici (...).
!   Aussi ce sous programme est-il appele dans phyvar (dans perinu perinr)
!   pour calculer les gradients au debut du pas de temps et les stocker
!   dans DUDXYZ et DRDXYZ

! Il est necessaire que ITENSO soit toujours initialise, meme hors
!   periodicite, donc on l'initialise au prealable a sa valeur par defaut.

idimte = 1
itenso = 0

if(iperio.eq.1) then

!       On recupere d'abord certains pointeurs necessaires a PERING

    call pergra                                                   &
    !==========
  ( nphsmx , nphas  ,                                             &
    iiu    , iiv    , iiw    ,                                    &
    iitytu ,                                                      &
    iir11  , iir22  , iir33  , iir12  , iir13  , iir23  )

  call pering                                                     &
  !==========
  ( nphas  , ivar   ,                                             &
    idimte , itenso , iperot , iguper , igrper ,                  &
    iiu    , iiv    , iiw    , iitytu ,                           &
    iir11  , iir22  , iir33  , iir12  , iir13  , iir23  ,         &
    dpdx   , dpdy   , dpdz   ,                                    &
    ra(idudxy) , ra(idrdxy)  )
endif

!===============================================================================
! 1. CALCUL DU GRADIENT
!===============================================================================


call cgdcel                                                       &
!==========
 ( ncelet , ncel   , nfac   , nfabor , ncelbr , ivar   ,          &
   imrgra , inc    , iccocg , nswrgp , idimte , itenso , iphydp , &
   iwarnp , nfecra , imligp , epsrgp , extrap , climgp ,          &
   ifacel , ifabor , ia(iicelb) , ia(iisymp) ,                    &
   volume , surfac , surfbo , ra(isrfbn) , ra(ipond),             &
   ra(idist) , ra(idistb) , ra(idijpf) , ra(idiipb) , ra(idofij) ,&
   fextx  , fexty  , fextz  ,                                     &
   xyzcen , cdgfac , cdgfbo, coefap , coefbp , pvar   ,           &
   ra(icocgb) , ra(icocg)   ,                                     &
   ra(icocib) , ra(icoci)   ,                                     &
   dpdx   , dpdy   , dpdz   ,                                     &
   dpdxa  , dpdya  , dpdza  )


return
end subroutine
