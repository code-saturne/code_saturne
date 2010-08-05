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

subroutine ppprcl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfppp ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   rcodcl , coefu  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    PREPARATION DU REMPLISSAGE DES CONDITIONS AUX LIMITES

!           AIGUILLAGE SPECIFIQUE AUX PHYSIQUES PARTICULIERES


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
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! coefu            ! tr ! --- ! tab de trav                                    !
!  nfabor,3        !    !     !  (vitesse en i'                 )              !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
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

include "paramx.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "cstnum.f90"
include "entsor.f90"
include "pointe.f90"
include "ppppar.f90"
include "ppthch.f90"
include "coincl.f90"
include "cpincl.f90"
include "fuincl.f90"
include "ppincl.f90"
include "cfpoin.f90"
include "atincl.f90"

!===============================================================================

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
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
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
double precision coefu(nfabor,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, izone, icha, iclapc
integer          iphas, ivar

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

! ---> Combustion gaz USEBUC
!      Flamme de diffusion : chimie 3 points

if ( ippmod(icod3p).ge.0 ) then

  do izone = 1, nozppm
    qimp(izone)   = zero
    iqimp(izone)  = 0
    ientox(izone) = 0
    ientfu(izone) = 0
  enddo

  tinoxy        = zero
  tinfue        = zero

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion gaz USEBUC
!      Flamme de premelange : modele EBU

elseif ( ippmod(icoebu).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    qimp(izone)   = zero
    icalke(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    fment(izone)  = zero
    tkent(izone)  = zero
    ientgf(izone) = 0
    ientgb(izone) = 0
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo


! ---> Combustion charbon pulverise USCPCL

elseif ( ippmod(icp3pl).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientcp(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    do icha = 1, ncharm
      qimpcp(izone,icha) = zero
      timpcp(izone,icha) = zero
      do iclapc = 1, ncpcmx
        distch(izone,icha,iclapc) = zero
      enddo
    enddo
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion charbon pulverise couple Lagrangien USCPLC

elseif ( ippmod(icpl3c).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    do icha = 1, ncharm
      qimpcp(izone,icha) = zero
    enddo
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion fuel  USFUCL

elseif ( ippmod(icfuel).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    ientfl(izone) = 0
    qimpfl(izone) = zero
    timpfl(izone) = zero
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Compressible

elseif ( ippmod(icompf).ge.0 ) then

!     Zones
  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

!     Marqueur d'utilisation de Rusanov au bord (0 = non)
!     Marqueur de flux conductif impose au bord (0 = non)
  do iphas = 1, nphas
    do ifac = 1, nfabor
      ia(iifbru+ifac-1+(iphas-1)*nfabor) = 0
      ia(iifbet+ifac-1+(iphas-1)*nfabor) = 0
    enddo
  enddo

!     Flux de Rusanov au bord pour Qdm et E
  do iphas = 1, nphas
    do ifac = 1, nfabor
      propfb(ifac,ipprob(ifbrhu(iphas))) = 0.d0
      propfb(ifac,ipprob(ifbrhv(iphas))) = 0.d0
      propfb(ifac,ipprob(ifbrhw(iphas))) = 0.d0
      propfb(ifac,ipprob(ifbene(iphas))) = 0.d0
    enddo
  enddo

!     Initialisation des RCODCL(IFAC,.,1) à -RINFIN
!       pour savoir si l'utilisateur les a modifies (ils sont
!       initialises par defaut à 0)
  do ivar = 1, nvar
    do ifac = 1, nfabor
      rcodcl(ifac,ivar,1) =-rinfin
    enddo
  enddo

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

elseif ( ippmod(ieljou).ge.1 .or.                                 &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Version ecoulements atmospheriques

elseif ( ippmod(iatmos).ge.0  ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo
  do izone = 1, nozppm
    iprofm(izone) = 0
  enddo

!     Initialisation des RCODCL(IFAC,.,1) à RINFIN
!       pour savoir si l'utilisateur les a modifies (ils sont
!       initialises par defaut à 0)
  do ivar = 1, nvar
    do ifac = 1, nfabor
      rcodcl(ifac,ivar,1) = rinfin
    enddo
  enddo

! ---> Version aerorefrigerants

elseif ( ippmod(iaeros).ge.0 ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

endif

!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
