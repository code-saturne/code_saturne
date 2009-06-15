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

subroutine raypun &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , iphas  ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   cofrua , cofrub ,                                              &
   flurds , flurdb ,                                              &
   dtr    , viscf  , viscb  ,                                     &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   theta4 , thetaa , sa     ,                                     &
   qx     , qy     , qz     ,                                     &
   qincid , eps    , tparoi ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , ckmel  ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS
!   AVEC L'APPROXIMATION P-1

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
! iphas            ! e  ! <-- ! numero de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
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
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
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
! cofrua,cofrub    ! tr ! --- ! conditions aux limites aux                     !
!(nfabor)          !    !     !    faces de bord pour la luminance             !
! flurds,flurdb    ! tr ! --- ! pseudo flux de masse (faces internes           !
!(nfac)(nfabor)    !    !     !    et faces de bord )                          !
! dtr(ncelet)      ! tr ! --- ! dt*cdtvar                                      !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! theta4(ncelet    ! tr ! --- ! pseudo temperature radiative                   !
! thetaa(ncelet    ! tr ! --- ! pseudo temp rar pdt precedent (nulle)          !
! sa (ncelet)      ! tr ! --> ! part d'absorption du terme source rad          !
! qxqyqz(ncelet    ! tr ! --> ! composante du vecteur densite de flux          !
!                  !    !     ! radiatif explicite                             !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
! eps (nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
! tparoi(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ckmel(ncelet)    ! tr ! <-- ! coeff d'absorption du melange                  !
!                  !    !     !   gaz-particules de charbon                    !
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
include "cstnum.h"
include "pointe.h"
include "radiat.h"
include "parall.h"
include "period.h"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor,nphas)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision cofrua(nfabor), cofrub(nfabor)
double precision flurds(nfac), flurdb(nfabor)

double precision dtr(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)

double precision theta4(ncelet), thetaa(ncelet)
double precision sa(ncelet)
double precision qx(ncelet), qy(ncelet), qz(ncelet)
double precision qincid(nfabor), tparoi(nfabor), eps(nfabor)

double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ckmel(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

character*80     cnom

integer          idebia, idebra
integer          ifac  , iel
integer          iconv1, idiff1, ndirc1, ireso1
integer          nitmap, nswrsp, nswrgp, iwarnp
integer          imgr1 , imligp, ircflp, ischcp, isstpp, iescap
integer          ncymap, nitmgp
integer          inum
integer          idtva0, ivar0
integer          inc, iccocg, iphydp
integer          idimte , itenso
double precision epsrgp, blencp, climgp, epsilp, extrap, epsrsp
double precision aa, aaa, aaaa, relaxp, thetap


!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. PARAMETRAGE DU SOLVEUR ET INITIALISATION
!===============================================================================

!--> Gradient Conjugue

ireso1 = 0

!--> Parametrage de CODITS

! IVAR0= 0  LA VARIABLE N'EST ICI NI RIJ NI VITESSE
ivar0   = 0
nitmap  = 1000
!     IMRGRA  = 0
nswrsp  = 1
nswrgp  = 100
imligp  = -1
ircflp  = 1
ischcp  = 1
isstpp  = 0
iescap  = 0
imgr1   = 1
ncymap  = 100
nitmgp  = 10
iwarnp  = iimlum
blencp  = zero
epsilp  = 1.d-8
epsrsp  = 1.d-8
epsrgp  = 1.d-5
climgp  = 1.5d0
extrap  = zero
relaxp  = 1.d0

!--> Il y a des dirichlets

ndirc1 = 1

!--> Pas de convection pour le modele P1

iconv1 = 0

!--> Equation de diffusion

idiff1 = 1

!--> Remise a zero des tableaux avant resolution

do iel = 1,ncel
  drtp(iel)   = zero
  theta4(iel) = zero
  thetaa(iel) = zero
enddo

do ifac = 1,nfac
  flurds(ifac) = zero
enddo

do ifac = 1,nfabor
  flurdb(ifac) = zero
enddo

!===============================================================================
! 2. COEFFICIENT DE DIFFUSION AUX FACES
!===============================================================================

do iel = 1,ncel
  ckmel(iel) = 1.d0 / ckmel(iel)
enddo

call viscfa                                                       &
!==========
   ( idebia , idebra ,                                            &
     ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,        &
     nprfml , nnod   , lndfac , lndfbr , ncelbr ,                 &
     nideve , nrdeve , nituse , nrtuse , imvisf ,                 &
     ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                 &
     ipnfac , nodfac , ipnfbr , nodfbr ,                          &
     idevel , ituser , ia     ,                                   &
     xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,        &
     volume , ckmel  , viscf  , viscb  ,                          &
     rdevel , rtuser , ra     )

!===============================================================================
! 3.  RESOLUTION
!===============================================================================

!     Parametre pour schemas en temps et stationnaire
thetap = 1.d0
idtva0 = 0

CNOM = ' '
WRITE(CNOM,'(A)') 'Rayon P1'
inum = 1
nomvar(inum) = cnom

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml , nprfml ,   &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtva0 , ivar0  , iconv1 , idiff1 , ireso1 , ndirc1 , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgr1  , ncymap , nitmgp , inum   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   volume ,                                                       &
   thetaa , thetaa , cofrua , cofrub , cofrua , cofrub ,          &
   flurds , flurdb ,                                              &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , theta4 ,                                     &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! 4. Vecteur densite de flux radiatif
!===============================================================================

!    En periodique et parallele, echange avant calcul du gradient

!    Parallele
if (irangp.ge.0) then
  call parcom (theta4)
  !==========
endif

!    Periodique
if (iperio.eq.1) then
  idimte = 0
  itenso = 0
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    theta4 , theta4 , theta4 ,                                    &
    theta4 , theta4 , theta4 ,                                    &
    theta4 , theta4 , theta4)
endif

!     Calcul de la densite du flux radiatif QX, QY, QZ

inc     = 1
iccocg  = 1
imligp  = -1
iwarnp  = iimlum
epsrgp  = 1.d-8
climgp  = 1.5d0
extrap  = 0.d0
nswrgp  = 100
ivar0   = 0
iphydp  = 0

call grdcel                                                       &
!==========
   ( idebia , idebra ,                                            &
     ndim   , ncelet , ncel   , nfac   , nfabor , nfml,           &
     nprfml ,                                                     &
     nnod   , lndfac , lndfbr , ncelbr , nphas  ,                 &
     nideve , nrdeve , nituse , nrtuse ,                          &
     ivar0  , imrgra , inc    , iccocg , nswrgp , imligp,         &
     iphydp ,                                                     &
     iwarnp , nfecra , epsrgp , climgp , extrap ,                 &
     ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                 &
     ipnfac , nodfac , ipnfbr , nodfbr ,                          &
     idevel , ituser , ia     ,                                   &
     xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,        &
     volume ,                                                     &
     w7     , w7     , w7     ,                                   &
     theta4 , cofrua , cofrub ,                                   &
     w1     , w2     , w3     ,                                   &
     w4     , w5     , w6     ,                                   &
     rdevel , rtuser , ra     )

aa = - stephn * 4.d0 / 3.d0

do iel = 1,ncel
  aaa = aa * ckmel(iel)
  qx(iel) = w1(iel) * aaa
  qy(iel) = w2(iel) * aaa
  qz(iel) = w3(iel) * aaa
enddo

!===============================================================================
! 5. Terme Source Radiatif d'absorption et densite de flux incident
!===============================================================================

!     Calcul de la part d'absorption du terme Source Radiatif

aa = 4.d0 * stephn
do iel = 1,ncel
  sa(iel) = aa * theta4(iel)
enddo

!     Calcul du flux incident Qincid

do ifac = 1, nfabor

  if (itypfb(ifac,iphas).eq.iparoi .or.                           &
      itypfb(ifac,iphas).eq.iparug ) then

!--> Premiere version plus chere et legerement plus precise

    aaaa = tparoi(ifac)**4

    aaa  = 1.5d0 * ra(idistb-1+ifac) / ckmel(ifabor(ifac))        &
           * ( 2.d0 /(2.d0-eps(ifac)) -1.d0 )
    aa   = ( aaa * aaaa + theta4(ifabor(ifac)) )                  &
         / (1.d0 + aaa)

    qincid(ifac) = stephn * (2.d0 * aa - eps(ifac) * aaaa)        &
                       / (2.d0 - eps(ifac))

!--> Deuxieme version plus cheap mais moins precise

!         QINCID(IFAC) = STEPHN *
!    &    (2.D0 * THETA4(IFABOR(IFAC)) - EPS(IFAC) * TPAROI(IFAC)**4)
!    &  / (2.D0 - EPS(IFAC))

  else
    qincid(ifac) = stephn * theta4(ifabor(ifac))                  &
               + ( qx(iel) * surfbo(1,ifac) +                     &
                   qy(iel) * surfbo(2,ifac) +                     &
                   qz(iel) * surfbo(3,ifac) ) /                   &
                   (0.5d0 * ra(isrfbn-1+ifac) )
  endif

enddo

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end
